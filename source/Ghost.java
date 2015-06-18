
/**
 * Write a description of class Ghost here.
 * UWA CITS3001 Project
 * @author Mitchell Poole and Sally Sorrell 
 * @version Sun 31/5/2015 Final version 
 */

import java.awt.Point;
import java.util.*;
public class Ghost
{

    //Data Structures for A* and Expectimax
    private HashMap<Point, Double> gScore;
    private HashMap<Point, Double> fScore;
    private HashMap<Point, Point> DecisionMap;
    private HashMap<Double, Point> ScoreMap;

    private Mode            GhostMode;
    private Mode            PreviousMode;
    private Orientation     GhostOrientation;
    private Orientation     PreviousOrientation;
    private int[]           GhostScatterTarget;
    private Target          GhostTarget;
    private boolean         dead;
    
    private int             GhostSpeed;

    private java.awt.Color  GhostColour;
    private Point           GhostPosition;
    private Point           OriginalPosition; //to go back home when dead
    
    //given
    public static enum Mode{CHASE,SCATTER,PANIC};
    public static enum Orientation{UP,DOWN,LEFT,RIGHT};
    public static enum Target {PACMAN,OFFSET,AIMLESS,SCATTER};
    public static final int[] SPEED={0,1,2,4,8,16,32};
    public static final int OFFSET = 4;
    
    public static final java.awt.Color PANIC_COLOUR = new java.awt.Color(60,100,175);
    
    /*
     * The SPEED indicates the legal speed, representing how many pixels the center of the Pac-Man or the ghost moves per frame. 
     * Note that these speed are proper divisors of the MazeViewer's CELL_SIZE.
     */

    public Ghost(java.awt.Point pos, java.awt.Color colour, int[] scatterTarget){
        /* creates a ghost with a specified position
         * colour
         * a two element integer array specifying the tile position of the scatter target
         * 
         * and sets its initial orientation to Orientation.UP
         * mode to Mode.SCATTER
         * initial travelling speed to CELL_SIZE. 
         */
        
        GhostPosition = pos;
        GhostColour = colour;
        GhostTarget = Target.SCATTER;
        GhostScatterTarget = scatterTarget;
        GhostOrientation = Orientation.UP;
        PreviousOrientation = GhostOrientation;
        setSpeed(MazeViewer.CELL_SIZE);
        GhostMode = Mode.SCATTER;
        
        dead = false;
        OriginalPosition = GhostPosition;
        PreviousMode = GhostMode;
				
	fScore = new HashMap<Point, Double>();
	gScore = new HashMap<Point, Double>();
    }
    
    public Ghost(java.awt.Point pos, java.awt.Color colour, Orientation ori, int speed, Mode m, Target t, int[] scatterTarget) {
       
        /* 
         * creates a ghost with client specified position, colour, orientation, speed, mode, 
         * target scheme and a two element integer array specifying the tile position of the scatter target. 
         * 
         * Your code should ensure the consistency between the specified mode and target scheme. 
         * Note, only Target.SCATTER can be used in Mode.SCATTER, only Target.AIMLESS can be used in Mode.PANIC.
         */
        
        GhostPosition = pos;
        GhostColour = colour;
        GhostOrientation = ori;
        setSpeed(speed);
        setMode(m, t);
        GhostScatterTarget = scatterTarget;
        dead = false;
        OriginalPosition = GhostPosition;
        PreviousOrientation = GhostOrientation;
        PreviousMode = GhostMode;

	fScore = new HashMap<Point, Double>();
	gScore = new HashMap<Point, Double>();
	DecisionMap = new HashMap<Point, Point>();
	ScoreMap  = new HashMap<Double, Point>();
    }
    
    public Mode getMode(){
        //returns the current mode of the ghost
        return GhostMode;
    }

    public int getSpeed() {
        //returns the travelling speed of the ghost.
        return GhostSpeed;
    }
    
    public void setSpeed(int newSpeed){
        //changes the travelling speed to newSpeed.
        
        //if in speed array
        boolean found = false;
        for(int i = 0; i < SPEED.length; i++){
            if (SPEED[i] == newSpeed) {
                found = true;
                break;
            }
        }
        
        if (found) {
            GhostSpeed = newSpeed;
        } else {
            throw new IllegalArgumentException("Error - Newspeed not in speed array");
        }
            
    }
    
    public java.awt.Color getColour(){
        //returns the colour of the ghost. 
        return GhostColour;
    }
    
    public void setColour(java.awt.Color colour){
        //sets the colour of the ghost. 
        GhostColour = colour;
    }
    
    public void setMode(Mode m, Target t){
        //changes the mode to m and target scheme to t. Again your code should ensure consistencies between the input.
        
       if (
            (t == Target.SCATTER && m == Mode.SCATTER) ||
            (t == Target.AIMLESS && m == Mode.PANIC)||
            ((t == Target.PACMAN || t == Target.OFFSET) && m == Mode.CHASE)
        ){
            GhostMode = m;  
            GhostTarget = t;
        } else {
            throw new IllegalArgumentException("Inconsitency between mode & target scheme");
        }        
    }
    
    public Orientation getOrientation(){
        //returns the current orientation of the ghost.
        return GhostOrientation;
    }

    public Orientation getPrevOrientation(){
        //returns the current orientation of the ghost.
        return PreviousOrientation;
    }		
    
    public Point getPosition(){
        //returns the current position of the ghost.
        return GhostPosition;
    }
    
    public boolean isDead(){
        return dead;
    }
    
    public boolean isPanic(){
        //returns true if the ghost is frightened and in PANIC mode.
        return GhostMode == Mode.PANIC;
    }
    
    public static void setPanic(Ghost[] ghosts, boolean panic){
        //sets a set of ghosts in panic mode or release them from panic mode. Remember setting the ghost in PANIC means the target scheme will be AIMLESS. 
        for (int i = 0; i < ghosts.length; i++){
            if(panic) {
                if(ghosts[i].GhostMode != Mode.PANIC)
                    ghosts[i].PreviousMode = ghosts[i].GhostMode;
                ghosts[i].GhostMode = Mode.PANIC;
                ghosts[i].setSpeed(MazeViewer.CELL_SIZE/2);
            } else {
                //return to previous !!!
                ghosts[i].GhostMode = ghosts[i].PreviousMode;
                ghosts[i].setSpeed(MazeViewer.CELL_SIZE);
                
                ghosts[i].GhostPosition.x = ((int) (ghosts[i].GhostPosition.x/2))*2; //round to nearest 2 pixels
                ghosts[i].GhostPosition.y = ((int) (ghosts[i].GhostPosition.y/2))*2; //round to nearest 2 pixels
                
            }
        }
    }
    
    public boolean atGrid(){
        boolean canx = false;
        if ((GhostPosition.x-MazeViewer.CELL_SIZE/2)%MazeViewer.CELL_SIZE == 0) canx = true;
        boolean cany = false;
        if ((GhostPosition.y-MazeViewer.CELL_SIZE/2)%MazeViewer.CELL_SIZE == 0) cany = true;

        return canx && cany;
    }
    
    public void move(Maze maze){
        //automatically moves the ghost in the user specified maze to the next position with proper orientation.
				
				/*
				Use alpha beta pruning
				
				*/
        boolean doMove = false;
        
        PacMan pacman = maze.getPacMan();
        
        //conditions to collide
        
        if((GhostPosition.x == pacman.getPosition().x && GhostPosition.y == pacman.getPosition().y || 
        
        Math.abs(GhostPosition.x - pacman.getPosition().x)== 1 && GhostPosition.y == pacman.getPosition().y ||
        GhostPosition.x == pacman.getPosition().x && Math.abs(GhostPosition.y - pacman.getPosition().y) == 1 || 
        
        Math.abs(GhostPosition.x - pacman.getPosition().x)== 2 && GhostPosition.y == pacman.getPosition().y ||
        GhostPosition.x == pacman.getPosition().x && Math.abs(GhostPosition.y - pacman.getPosition().y) == 2 ) && !dead){
            
            maze.doCollide(isPanic(), GhostPosition);
            if(isPanic()){
                dead = true;
            }
        }	//makes it dead
        
        if(atGrid()){
            
            Orientation[] Ori = getOrientations(maze);
            
           if(dead){
                GhostOrientation = targetDirection(OriginalPosition,maze);
                
                if(Math.abs(GhostPosition.x - OriginalPosition.x) <= 4 && Math.abs(GhostPosition.y - OriginalPosition.y) <= 4){
                    dead = false;
                }
            } else {
                switch(GhostMode){
                    case PANIC:
                    //random from possible directions 
                    GhostOrientation = Ori[(int) (Math.random()*Ori.length)];
                    break;
                    
                    case SCATTER: 
                    //follow scatter target
                    Point target = new Point (GhostScatterTarget[0], GhostScatterTarget[1]);
	            	//GhostOrientation = targetDirection(target, maze, true);
			
			//TEST LINES 
			GhostOrientation = AStar(maze);
			//GhostOrientation = ExpectiMax(maze, GhostPosition, 4 );
                    break;
                    
                    case CHASE:
                    if(GhostTarget == Target.PACMAN){
                        //follow pacman
                        //GhostOrientation = targetDirection(maze.getPacMan().getPosition(), maze);
			
			//Test Lines 
			GhostOrientation = AStar(maze);
			//GhostOrientation = ExpectiMax(maze, GhostPosition, 4 );		
                    } else {
                        //offset
                        //GhostOrientation = targetDirection(offset(maze.getPacMan().getPosition(),maze.getPacMan().getOrientation()), maze);
			//More test lines 
			GhostOrientation = AStar(maze);
			//GhostOrientation = ExpectiMax(maze, GhostPosition, 4 );	
                    }
			break;
                }
            }

            switch (maze.locationStatus(nextPos(GhostPosition,GhostOrientation))){
            case INVALID:
            //wrap around to other side of screen
                if (GhostPosition.x == MazeViewer.CELL_SIZE/2){ //left side
                    GhostPosition.x = maze.getMap().length*MazeViewer.CELL_SIZE-MazeViewer.CELL_SIZE/2;
                    
                } else if (GhostPosition.x == maze.getMap().length*MazeViewer.CELL_SIZE-MazeViewer.CELL_SIZE/2){ //right side
                    GhostPosition.x = MazeViewer.CELL_SIZE/2;
                }
            case DEAD:
                //won't happen (literally what)
                break;
            case ENERGISER:
            case LEGAL:
            case DOT:
            doMove = true;
            }
        } else {
            doMove = true;
        }
        
        if(doMove){
             PreviousOrientation = GhostOrientation;
             GhostPosition = nextPos2(GhostPosition, GhostOrientation);
        }
    }

    private Orientation targetDirection(Point target, Maze maze){
        return targetDirection(target, maze, false);
    }
    
/* A cleaner form of A* Search Algorithm

*/
	/* Comparator for the F(x) = g(x) + h(x) */
    public class fComp implements Comparator<Point>{
    	@Override
    	public int compare(Point a, Point b) {
            if (fScore.get(a) > fScore.get(b) ) return 1;
	    else if(fScore.get(a) < fScore.get(b)) return -1;
	    else return 0;
	}	
    }	

	
    /*
    * Uses the possible Orientations to get the neighbors
    */
    private Point[] getNeighbors(Maze maze, Point Pos) {
    	Orientation[] posOri = getOrientations(maze);
    	int size = posOri.length;
    	Point[] points = new Point[size];
		
    	for(int i = 0; i < posOri.length; i++) { //ads points into array
            if(posOri[i] == Orientation.LEFT ) {
	        if(nextPos(Pos, Orientation.LEFT).x ==  MazeViewer.CELL_SIZE/2) { //wraps around special case
		    points[i] = nextPos(Pos, Orientation.LEFT);
		    points[i].x = maze.getMap().length*MazeViewer.CELL_SIZE-MazeViewer.CELL_SIZE/2;
		}
		else
		    points[i] = nextPos(Pos, Orientation.LEFT);	
		}
		if(posOri[i] == Orientation.RIGHT) {
		    if ( nextPos(Pos, Orientation.RIGHT).x == maze.getMap().length*MazeViewer.CELL_SIZE-MazeViewer.CELL_SIZE/2) {
			points[i] = nextPos(Pos, Orientation.RIGHT);
		        points[i].x =  MazeViewer.CELL_SIZE/2;
		    }
		    else
		        points[i] = nextPos(Pos, Orientation.RIGHT);
		}
		if(posOri[i] == Orientation.UP)
		    points[i] = nextPos(Pos, Orientation.UP);
		if(posOri[i] == Orientation.DOWN)
		    points[i] = nextPos(Pos, Orientation.DOWN);
	}
	return points;
    }

	/*
	* Returns the Orientation B is in relation to A
	*/
    private Orientation whichOrientation(Point a, Point b) { //Which orientation is b in relation to a
    //	     UP
    // 	     .
    // LEFT .... a .... RIGHT 
    //	     .
    //	    DOWN
    	if(a.y == b.y && a.x < b.x) return Orientation.RIGHT;
    	else if( a.y == b.y && a.x > b.x) return Orientation.LEFT; 
	else if( a.x == b.x && a.y < b.y) return Orientation.DOWN;
	else return Orientation.UP;
    }
    /*
    * Calculates the gScore between two points
    */
    private double gScore(Point a, Point b) {
    	return Math.sqrt(Math.pow(Math.abs(a.x-b.x),2) + Math.pow(Math.abs(a.x-b.x),2));	}
    /*
    * Calculates the heuristic between two Points
    */
    private double heuristic(Point a, Point b) {
    	return Math.abs(a.x - b.x) + Math.abs(a.y - b.y);
    }

    /*
    * Uses the AStar algorithm to calculate the optimal path, then returns to orientation of the first two nodes
    * A PriorityQueue OpenList is used to store the unchecked nodes of the graph travesal
    * A HashSet Closed List is to store the checked or finished nodes of the graph
    * A Hashmap is used to store the F and G values
    * Implemented psudeocode from: http://en.wikipedia.org/wiki/A*_search_algorithm 
    */
    private Orientation AStar(Maze maze) {
	fScore.clear();
	gScore.clear();
		
	//To make sure the get neighors will equal Pacman eventually
	Point Pacman = maze.getPacMan().getPosition();
	Pacman.x = Pacman.x - Pacman.x%MazeViewer.CELL_SIZE + MazeViewer.CELL_SIZE/2;
	Pacman.y = Pacman.y - Pacman.x%MazeViewer.CELL_SIZE + MazeViewer.CELL_SIZE/2;
	Point start = GhostPosition;
	start.x = start.x - start.x%MazeViewer.CELL_SIZE + MazeViewer.CELL_SIZE/2; 
	start.y = start.y - start.x%MazeViewer.CELL_SIZE + MazeViewer.CELL_SIZE/2;

	PriorityQueue<Point> openList = new PriorityQueue<Point>(1, new fComp());
	HashSet<Point> closedList = new HashSet<Point>();
	HashMap<Point, Point> cameFrom = new HashMap<Point, Point>();
		
	Point currentPos = start;
		
	openList.add(currentPos);
	
	gScore.put(currentPos, 0.0);
	fScore.put(currentPos, heuristic(currentPos, Pacman));	
	
	while(!openList.isEmpty()) {
	    if(currentPos.equals(Pacman)) {
	        openList.clear();
	    }
	    else {
	    	currentPos = openList.poll();
	    	closedList.add(currentPos);
	    	Point[] neighbors = getNeighbors(maze, currentPos);
		for(int i = 0; i < neighbors.length; i++) { //The Status check because getNeighbors returns invalid Points
		    if(closedList.contains(neighbors[i]) || maze.locationStatus(neighbors[i]) == Maze.Status.INVALID ) //in closed set 
		    { continue; }
		    double tentativeG = gScore.get(currentPos) +  heuristic(currentPos, neighbors[i]);
		    if(!openList.contains(neighbors[i]) || tentativeG < gScore.get(neighbors[i])) {

			cameFrom.put(neighbors[i], currentPos);
			gScore.put(neighbors[i], tentativeG);
			fScore.put(neighbors[i], gScore.get(neighbors[i]) + heuristic(neighbors[i], Pacman));

			if(!openList.contains(neighbors[i]))
			    openList.add(neighbors[i]);
			}

		}
	    }
			 //some reason will go into an infinte loop 
	}
		// need 1st and second positions
	Point secondPos = new Point(0,0);
		
	while(!currentPos.equals(start)) {
	    if(cameFrom.get(currentPos).equals(start))
	        secondPos = currentPos;
	    currentPos = cameFrom.get(currentPos);
	}
	return whichOrientation(start, secondPos);
    }		

	/** END OF A STAR WORKING

	**/
	
	/**
	* Expectimax
	*/


    private final int MAX_DEPTH = 6;
		
    private final double PAC_DIST_SCALE = 0.4; //DistP Evaluations 
    private final double ENERGISE_SCALE =  0.1; //EnerRem Evaluations
    private final double REM_DOT_SCALE = 0.15; //DotRem Evaluations 
    private final double CORNER_SCALE = 0.35; //TraP Evaluations

    private final double NEG_INFINITY = Double.NEGATIVE_INFINITY;
    private final double INFINITY = Double.POSITIVE_INFINITY;

    private int ENERGISER_COUNT = 0;

    private double evaluationFunction(Maze maze , Point node) {
    	double pacDist = PAC_DIST_SCALE*100/heuristic(node, maze.getPacMan().getPosition()); //better when dist is smaller 
	double trapping = CORNER_SCALE*trapping(maze, node); // both close to pacman but farther away from each other
	double enerDist = ENERGISE_SCALE*ENERGISER_COUNT*50; //the more around the better
	double remDots = REM_DOT_SCALE*dotsRem(maze); //same with energisers 
	return pacDist + trapping + enerDist + remDots;
    }


    /*
    * Gets the number of both dots and energizers 
    */
    private int dotsRem(Maze maze) {
    	Maze.Status[][] dots = maze.getMap();
	int dotCount=0;
	for(int i = 0; i < dots.length; i++) {
	    for( int j = 0; j< dots[0].length; j++) {
	        if(dots[i][j] == Maze.Status.DOT)
		    dotCount++;
	        else if (dots[i][j] == Maze.Status.ENERGISER)
		    ENERGISER_COUNT++;
	    }
	}
	return dotCount;	
    }



    /*
    * Compares the distances of the ghosts relative to pacman and pacman's orientations
    *
    */
    private int trapping(Maze maze, Point node) {
	Ghost[] ghosts = maze.getGhosts();
	int trapCount = 0;
	for(int i = 0; i < ghosts.length; i++) {
	    double pacDist = heuristic(ghosts[i].getPosition(), maze.getPacMan().getPosition());
	    for(int j=0;j<ghosts.length;j++) {
	    	if(ghosts[i].equals(ghosts[j])) { continue; }
		if(pacDist < heuristic(ghosts[j].getPosition(), maze.getPacMan().getPosition()) ) {
		    if(maze.getPacMan().getOrientations(maze).length <2) //Pacman can only go straight
			trapCount += 100; //this is arbitrary
		    }
	    }
	}
	return trapCount;
    }
    /*
    * Uses the returned value of the Expectimax to Recurse through the hashmap storing the path
    */
    private Orientation ExpectiMax(Maze maze, Point start, int depth ) { //to return a Orientation/Playstate
	DecisionMap = new HashMap<Point,Point>();
	ScoreMap = new HashMap<Double,Point>();

	Point backTrack = ScoreMap.get( expectiMax(maze, start, depth, true));
	Point second = new Point(0,0);
	while(!backTrack.equals(start)) {
	    if(DecisionMap.get(backTrack).equals(start)) 
		second = backTrack;
	    backTrack = DecisionMap.get(backTrack); 
	    
	}
	DecisionMap.clear();
	ScoreMap.clear();
	return whichOrientation(backTrack, second);
    }
    /*
    *
    */
    private double expectiMax(Maze maze, Point node, int depth, boolean maxPlay) throws IllegalArgumentException {
	if(depth > MAX_DEPTH) {
	    throw new IllegalArgumentException();
	}
	double alpha = 0;
	if(depth == 0) { //no terminal nodes since Ghosts and Pacman can keep moving
	    return evaluationFunction(maze, node);
	}	
	if (maxPlay) { //
	    alpha = INFINITY;
	    Point[] children = getNeighbors(maze, node);
	    for(int i = 0; i < children.length; i++) {
	        if(maze.locationStatus(children[i]) == Maze.Status.INVALID) { 
	        }
		else 
		    alpha = Math.min(alpha,  expectiMax(maze, children[i], depth-1,false));
		    ScoreMap.put(alpha,children[i]);
		    DecisionMap.put(children[i], node);
		}
	}
	else { //
	    Point[] children = getNeighbors(maze, node);
	    for(int i = 0; i < children.length; i++) {
	    	if(maze.locationStatus(children[i]) == Maze.Status.INVALID) { 
		}
		else 
		    alpha = (double)(1/children.length)*expectiMax(maze, children[i], depth-1,true);
		    ScoreMap.put(alpha, children[i]);
		    DecisionMap.put(children[i],node);
		}
	}
	return alpha;
    }
		
		
    /**
    * End of Expecitmax function
    * 
    *
    **/


    private Orientation targetDirection(Point target, Maze maze, boolean doScale){
        
        Orientation[] possibleOrientations = getOrientations(maze);

        if(possibleOrientations.length == 1){
            return possibleOrientations[0];
        }

        Orientation thisway = possibleOrientations[0];        
        int dist;
        int tempDist;
        
        if(doScale){
            dist = Math.abs(target.x*MazeViewer.CELL_SIZE - nextPos(GhostPosition,possibleOrientations[0]).x) + Math.abs(target.y*MazeViewer.CELL_SIZE - nextPos(GhostPosition,possibleOrientations[0]).y);
        } else {
            dist = Math.abs(target.x - nextPos(GhostPosition,possibleOrientations[0]).x) + Math.abs(target.y - nextPos(GhostPosition,possibleOrientations[0]).y);
        }

        //find shortest route
        for(int i = 0; i < possibleOrientations.length; i++){
            
            if(doScale){ 
                tempDist = Math.abs(target.x*MazeViewer.CELL_SIZE - nextPos(GhostPosition,possibleOrientations[i]).x) + Math.abs(target.y*MazeViewer.CELL_SIZE - nextPos(GhostPosition,possibleOrientations[i]).y);
            }
	    else {
                tempDist = Math.abs(target.x - nextPos(GhostPosition,possibleOrientations[i]).x) + Math.abs(target.y - nextPos(GhostPosition,possibleOrientations[i]).y);
            }
            
            if(dist > tempDist){
                dist = tempDist;
                thisway = possibleOrientations[i];
            }
            /*
            if(tempDist <= MazeViewer.CELL_SIZE){
                if(Math.random()<0.9){
                    return possibleOrientations[i];
                }
            }
            */
        }

        return thisway;
    }
    
    private Point nextPos(Point positioncopy, Orientation ori){
        Point position = new Point(positioncopy); //translate copy of point
        int dist = MazeViewer.CELL_SIZE; //+-16
        switch (ori){
            case UP: 
                position.translate(0,-dist); break;
            case DOWN: 
                position.translate(0,dist); break;
            case LEFT: 
                position.translate(-dist,0); break;
            case RIGHT: 
                position.translate(dist,0); break;
        }
        return position;
    }
    
    private Point nextPos2(Point positioncopy, Orientation ori){
        Point position = new Point(positioncopy); //translate copy of point
        int dist = GhostSpeed/(MazeViewer.CELL_SIZE/2); //ghost speed
        switch (ori){
            case UP: 
                position.translate(0,-dist); break;
            case DOWN: 
                position.translate(0,dist); break;
            case LEFT: 
                position.translate(-dist,0); break;
            case RIGHT: 
                position.translate(dist,0); break;
        }
        return position;
    }
    
    private Point offset(Point positioncopy, Orientation ori){
        Point position = new Point(positioncopy); //translate copy of point
        int dist = MazeViewer.CELL_SIZE*4; //for offset
        switch (ori){
            case UP: 
                position.translate(0,-dist); break;
            case DOWN: 
                position.translate(0,dist); break;
            case LEFT: 
                position.translate(-dist,0); break;
            case RIGHT: 
                position.translate(dist,0); break;
        }
        return position;
    }
    
    private Orientation[] getOrientations(Maze maze){
        //up down left right
        
        boolean[] temp = new boolean[4];
        
        for(int i = 0; i < 4; i++){
            temp[i] = false; //set all directions not possible
        }
        
        if(!(maze.locationStatus(nextPos(GhostPosition, Orientation.UP)) == Maze.Status.DEAD)){ //if maze is not dead    
            if(GhostOrientation != Orientation.DOWN){ //and haven't just come that way
                temp[0] = true; //direction is possible
            }
        } 
        
        if(!(maze.locationStatus(nextPos(GhostPosition, Orientation.DOWN)) == Maze.Status.DEAD)){ 
            if(GhostOrientation != Orientation.UP){
                temp[1] = true;
            }
        } 
        
        if(!(maze.locationStatus(nextPos(GhostPosition, Orientation.LEFT)) == Maze.Status.DEAD)){
            if(GhostOrientation != Orientation.RIGHT && PreviousOrientation != Orientation.RIGHT){
                temp[2] = true;
            }
        } 
        
        if(!(maze.locationStatus(nextPos(GhostPosition, Orientation.RIGHT)) == Maze.Status.DEAD)){
            if(GhostOrientation != Orientation.LEFT && PreviousOrientation != Orientation.LEFT){
                temp[3] = true;
            }
        }
        
        int count = 0;
        
        for(int i = 0; i < 4; i++){
            if(temp[i]){
                count++; //count possible directions
            }
        }
        
        Orientation[] tempOri = new Orientation[count]; //new array of possible directions
        
        int i = 0;
        
        if(temp[0]){
            tempOri[i]=Orientation.UP; //populate array
            i++;
        }
        
        if(temp[1]){
            tempOri[i]=Orientation.DOWN;
            i++;
        }
        
        if(temp[2]){
            tempOri[i]=Orientation.LEFT;
            i++;
        }
        
        if(temp[3]){
            tempOri[i]=Orientation.RIGHT;
        }
        
        return tempOri;
    }
    
    static public void setToOriginal(Ghost[] ghosts){
        
        for(int i = 0; i < ghosts.length; i++){
            
            ghosts[i].GhostPosition = ghosts[i].OriginalPosition;
            
        }
    }
}
