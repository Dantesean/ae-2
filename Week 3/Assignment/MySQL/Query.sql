SELECT 
	collisions.*,
	(SELECT temp FROM 
		(SELECT
		  temp, (
		    3959 * acos (
		      cos ( radians(collisions.latitude) )
		      * cos( radians(temps.lat ) )
		      * cos( radians(temps.lon ) - radians(collisions.longitude) )
		      + sin ( radians(collisions.latitude) )
		      * sin( radians(temps.lat ) )
		    )
		  ) AS distance
		FROM temps
		WHERE date = collision_date
		HAVING distance <= 10
		ORDER BY distance
		LIMIT 1) t1
	) t2 
FROM collisions
WHERE collision_date IS NOT NULL
AND longitude IS NOT NULL
AND latitude IS NOT NULL
AND collision_date <= '2017-12-31'  -- 2018 is incomplete
AND collision_date >= '2010-01-01'
LIMIT 10