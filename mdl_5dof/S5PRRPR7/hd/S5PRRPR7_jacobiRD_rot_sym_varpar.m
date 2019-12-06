% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRPR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:32
	% EndTime: 2019-12-05 16:38:32
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(5));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(5));
	t60 = cos(pkin(9));
	t58 = sin(pkin(9));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0; 0, -t62 * t64, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0; 0, -t63 * t64, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:33
	% EndTime: 2019-12-05 16:38:33
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (36->22), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t211 = sin(pkin(5));
	t214 = sin(qJ(3));
	t227 = t211 * t214;
	t216 = cos(qJ(3));
	t226 = t211 * t216;
	t213 = cos(pkin(5));
	t215 = sin(qJ(2));
	t225 = t213 * t215;
	t217 = cos(qJ(2));
	t224 = t213 * t217;
	t223 = qJD(2) * t215;
	t222 = qJD(3) * t214;
	t221 = qJD(3) * t216;
	t220 = qJD(3) * t217;
	t219 = t211 * qJD(2) * t217;
	t210 = sin(pkin(9));
	t212 = cos(pkin(9));
	t206 = -t210 * t215 + t212 * t224;
	t207 = t210 * t217 + t212 * t225;
	t208 = -t210 * t224 - t212 * t215;
	t218 = t210 * t225 - t212 * t217;
	t205 = t218 * qJD(2);
	t204 = t208 * qJD(2);
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t1 = [0, t205 * t216 - t208 * t222, -t204 * t214 + (-t210 * t227 + t216 * t218) * qJD(3), 0, 0; 0, -t203 * t216 - t206 * t222, -t202 * t214 + (-t207 * t216 + t212 * t227) * qJD(3), 0, 0; 0, (-t214 * t220 - t216 * t223) * t211, -t214 * t219 + (-t213 * t214 - t215 * t226) * qJD(3), 0, 0; 0, -t205 * t214 - t208 * t221, -t204 * t216 + (-t210 * t226 - t214 * t218) * qJD(3), 0, 0; 0, t203 * t214 - t206 * t221, -t202 * t216 + (t207 * t214 + t212 * t226) * qJD(3), 0, 0; 0, (t214 * t223 - t216 * t220) * t211, -t216 * t219 + (-t213 * t216 + t215 * t227) * qJD(3), 0, 0; 0, t204, 0, 0, 0; 0, t202, 0, 0, 0; 0, t219, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:33
	% EndTime: 2019-12-05 16:38:33
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (64->32), mult. (240->82), div. (0->0), fcn. (252->10), ass. (0->35)
	t299 = sin(pkin(5));
	t303 = sin(qJ(3));
	t320 = t299 * t303;
	t305 = cos(qJ(3));
	t319 = t299 * t305;
	t302 = cos(pkin(5));
	t304 = sin(qJ(2));
	t318 = t302 * t304;
	t306 = cos(qJ(2));
	t317 = t302 * t306;
	t316 = t303 * t304;
	t315 = t304 * t305;
	t314 = qJD(3) * t303;
	t313 = qJD(3) * t305;
	t312 = qJD(3) * t306;
	t311 = qJD(2) * t299 * t306;
	t310 = t303 * t312;
	t298 = sin(pkin(9));
	t301 = cos(pkin(9));
	t293 = -t298 * t304 + t301 * t317;
	t294 = t298 * t306 + t301 * t318;
	t295 = -t298 * t317 - t301 * t304;
	t309 = t298 * t318 - t301 * t306;
	t290 = t294 * qJD(2);
	t308 = t290 * t305 + t293 * t314;
	t292 = t309 * qJD(2);
	t307 = -t292 * t305 + t295 * t314;
	t300 = cos(pkin(10));
	t297 = sin(pkin(10));
	t291 = t295 * qJD(2);
	t289 = t293 * qJD(2);
	t288 = -t303 * t311 + (-t299 * t315 - t302 * t303) * qJD(3);
	t287 = -t291 * t303 + (-t298 * t320 + t305 * t309) * qJD(3);
	t286 = -t289 * t303 + (-t294 * t305 + t301 * t320) * qJD(3);
	t1 = [0, t291 * t297 - t307 * t300, t287 * t300, 0, 0; 0, t289 * t297 - t308 * t300, t286 * t300, 0, 0; 0, (-t300 * t310 + (t297 * t306 - t300 * t315) * qJD(2)) * t299, t288 * t300, 0, 0; 0, t291 * t300 + t307 * t297, -t287 * t297, 0, 0; 0, t289 * t300 + t308 * t297, -t286 * t297, 0, 0; 0, (t297 * t310 + (t297 * t315 + t300 * t306) * qJD(2)) * t299, -t288 * t297, 0, 0; 0, t292 * t303 + t295 * t313, t291 * t305 + (t298 * t319 + t303 * t309) * qJD(3), 0, 0; 0, -t290 * t303 + t293 * t313, t289 * t305 + (-t294 * t303 - t301 * t319) * qJD(3), 0, 0; 0, (-qJD(2) * t316 + t305 * t312) * t299, t305 * t311 + (-t299 * t316 + t302 * t305) * qJD(3), 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:34
	% EndTime: 2019-12-05 16:38:34
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (272->94), mult. (904->192), div. (0->0), fcn. (1000->12), ass. (0->69)
	t464 = sin(pkin(10));
	t475 = cos(qJ(2));
	t506 = t464 * t475;
	t466 = sin(pkin(5));
	t471 = sin(qJ(3));
	t505 = t466 * t471;
	t474 = cos(qJ(3));
	t504 = t466 * t474;
	t467 = cos(pkin(10));
	t470 = sin(qJ(5));
	t503 = t467 * t470;
	t473 = cos(qJ(5));
	t502 = t467 * t473;
	t501 = t467 * t474;
	t500 = t467 * t475;
	t469 = cos(pkin(5));
	t472 = sin(qJ(2));
	t499 = t469 * t472;
	t498 = t469 * t475;
	t497 = t472 * t474;
	t496 = qJD(2) * t472;
	t495 = qJD(2) * t475;
	t494 = qJD(3) * t471;
	t493 = qJD(3) * t474;
	t492 = qJD(3) * t475;
	t491 = qJD(5) * (t464 * t472 + t474 * t500) * t466;
	t490 = qJD(5) * t471;
	t489 = qJD(5) * t475;
	t465 = sin(pkin(9));
	t488 = t465 * t499;
	t468 = cos(pkin(9));
	t487 = t468 * t505;
	t486 = t466 * t495;
	t485 = t471 * t492;
	t484 = t474 * t492;
	t480 = -t465 * t472 + t468 * t498;
	t449 = t480 * qJD(2);
	t454 = t465 * t475 + t468 * t499;
	t450 = t454 * qJD(2);
	t479 = -t450 * t474 - t480 * t494;
	t483 = -t449 * t464 - t479 * t467 - t480 * t490;
	t455 = t465 * t498 + t468 * t472;
	t451 = t455 * qJD(2);
	t452 = -qJD(2) * t488 + t468 * t495;
	t478 = -t452 * t474 + t455 * t494;
	t482 = t451 * t464 + t455 * t490 - t478 * t467;
	t442 = t454 * t471 + t468 * t504;
	t456 = t468 * t475 - t488;
	t481 = -t456 * t471 + t465 * t504;
	t445 = t456 * t474 + t465 * t505;
	t458 = t466 * t497 + t469 * t471;
	t457 = -t469 * t474 + t472 * t505;
	t477 = -qJD(5) * (t454 * t464 + t480 * t501) - t450 * t471 + t480 * t493;
	t476 = -qJD(5) * (-t455 * t501 + t456 * t464) - t452 * t471 - t455 * t493;
	t447 = -t457 * qJD(3) + t474 * t486;
	t446 = t458 * qJD(3) + t471 * t486;
	t443 = t454 * t474 - t487;
	t441 = t458 * t467 - t466 * t506;
	t440 = (-t467 * t485 + (-t467 * t497 + t506) * qJD(2)) * t466;
	t437 = t466 * t464 * t496 + t447 * t467;
	t436 = t445 * t467 + t455 * t464;
	t435 = t443 * t467 - t464 * t480;
	t434 = t481 * qJD(3) - t451 * t474;
	t433 = t445 * qJD(3) - t451 * t471;
	t432 = -t442 * qJD(3) + t449 * t474;
	t431 = -qJD(3) * t487 + t449 * t471 + t454 * t493;
	t428 = t434 * t467 + t452 * t464;
	t427 = t432 * t467 + t450 * t464;
	t1 = [0, t476 * t470 - t482 * t473, -t433 * t502 + t434 * t470 + (t445 * t473 - t481 * t503) * qJD(5), 0, -t428 * t470 + t433 * t473 + (-t436 * t473 + t470 * t481) * qJD(5); 0, t477 * t470 - t483 * t473, -t431 * t502 + t432 * t470 + (t442 * t503 + t443 * t473) * qJD(5), 0, -t427 * t470 + t431 * t473 + (-t435 * t473 - t442 * t470) * qJD(5); 0, -t470 * t491 + t440 * t473 + (t471 * t473 * t489 + (-t471 * t496 + t484) * t470) * t466, -t446 * t502 + t447 * t470 + (t457 * t503 + t458 * t473) * qJD(5), 0, -t437 * t470 + t446 * t473 + (-t441 * t473 - t457 * t470) * qJD(5); 0, t482 * t470 + t476 * t473, t433 * t503 + t434 * t473 + (-t445 * t470 - t481 * t502) * qJD(5), 0, -t428 * t473 - t433 * t470 + (t436 * t470 + t473 * t481) * qJD(5); 0, t483 * t470 + t477 * t473, t431 * t503 + t432 * t473 + (t442 * t502 - t443 * t470) * qJD(5), 0, -t427 * t473 - t431 * t470 + (t435 * t470 - t442 * t473) * qJD(5); 0, -t473 * t491 - t440 * t470 + (t473 * t484 + (-t470 * t489 - t473 * t496) * t471) * t466, t446 * t503 + t447 * t473 + (t457 * t502 - t458 * t470) * qJD(5), 0, -t437 * t473 - t446 * t470 + (t441 * t470 - t457 * t473) * qJD(5); 0, t451 * t467 + t478 * t464, -t433 * t464, 0, 0; 0, -t449 * t467 + t479 * t464, -t431 * t464, 0, 0; 0, (-t464 * t485 + (-t464 * t497 - t500) * qJD(2)) * t466, -t446 * t464, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end