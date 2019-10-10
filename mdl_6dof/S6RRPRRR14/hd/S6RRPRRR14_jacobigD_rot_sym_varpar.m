% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR14_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobigD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:58
	% EndTime: 2019-10-10 11:10:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t106 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t106, 0, 0, 0, 0; 0, sin(qJ(1)) * t106, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:10:59
	% EndTime: 2019-10-10 11:10:59
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (28->18), mult. (102->43), div. (0->0), fcn. (102->12), ass. (0->29)
	t241 = cos(pkin(14));
	t243 = cos(pkin(7));
	t262 = t243 * t241;
	t245 = sin(qJ(2));
	t246 = sin(qJ(1));
	t261 = t245 * t246;
	t248 = cos(qJ(1));
	t260 = t245 * t248;
	t247 = cos(qJ(2));
	t259 = t246 * t247;
	t258 = t247 * t248;
	t240 = sin(pkin(6));
	t257 = qJD(1) * t240;
	t237 = sin(pkin(14));
	t256 = qJD(2) * t237;
	t239 = sin(pkin(7));
	t255 = t240 * t239 * t241;
	t254 = t246 * t257;
	t253 = t248 * t257;
	t244 = cos(pkin(6));
	t252 = t244 * t258 - t261;
	t251 = -t244 * t259 - t260;
	t250 = -t244 * t260 - t259;
	t249 = t244 * t261 - t258;
	t242 = cos(pkin(8));
	t238 = sin(pkin(8));
	t236 = t251 * qJD(1) + t250 * qJD(2);
	t235 = -t252 * qJD(1) + t249 * qJD(2);
	t1 = [0, t253, 0, -(t235 * t262 - t251 * t256 + (-t250 * t237 + t248 * t255) * qJD(1)) * t238 + (-t235 * t239 + t243 * t253) * t242, 0, 0; 0, t254, 0, -(t236 * t262 - t252 * t256 + (t249 * t237 + t246 * t255) * qJD(1)) * t238 + (-t236 * t239 + t243 * t254) * t242, 0, 0; 0, 0, 0, (-(-t237 * t247 - t245 * t262) * t238 + t245 * t239 * t242) * t240 * qJD(2), 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:11:01
	% EndTime: 2019-10-10 11:11:01
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (110->49), mult. (371->110), div. (0->0), fcn. (396->14), ass. (0->46)
	t354 = sin(pkin(7));
	t359 = cos(pkin(6));
	t384 = t354 * t359;
	t361 = sin(qJ(2));
	t383 = t354 * t361;
	t355 = sin(pkin(6));
	t362 = sin(qJ(1));
	t382 = t355 * t362;
	t365 = cos(qJ(1));
	t381 = t355 * t365;
	t358 = cos(pkin(7));
	t380 = t358 * t361;
	t364 = cos(qJ(2));
	t379 = t358 * t364;
	t378 = t362 * t361;
	t377 = t362 * t364;
	t376 = t365 * t361;
	t375 = t365 * t364;
	t374 = qJD(1) * t355;
	t373 = qJD(2) * t355;
	t372 = t362 * t374;
	t371 = t365 * t374;
	t348 = t359 * t375 - t378;
	t370 = t348 * t358 - t354 * t381;
	t350 = -t359 * t377 - t376;
	t369 = t350 * t358 + t354 * t382;
	t349 = t359 * t376 + t377;
	t368 = t359 * t378 - t375;
	t343 = -t348 * qJD(1) + t368 * qJD(2);
	t367 = t343 * t358 + t354 * t371;
	t345 = t350 * qJD(1) - t349 * qJD(2);
	t366 = t345 * t358 + t354 * t372;
	t363 = cos(qJ(4));
	t360 = sin(qJ(4));
	t357 = cos(pkin(8));
	t356 = cos(pkin(14));
	t353 = sin(pkin(8));
	t352 = sin(pkin(14));
	t347 = (-t352 * t364 - t356 * t380) * t373;
	t346 = -t368 * qJD(1) + t348 * qJD(2);
	t344 = -t349 * qJD(1) + t350 * qJD(2);
	t342 = -t345 * t354 + t358 * t372;
	t341 = -t343 * t354 + t358 * t371;
	t340 = -t346 * t352 + t366 * t356;
	t339 = -t344 * t352 + t367 * t356;
	t1 = [0, t371, 0, -t339 * t353 + t341 * t357, (t344 * t356 + t367 * t352) * t360 + (-t339 * t357 - t341 * t353) * t363 + ((t369 * t352 - t356 * t368) * t363 + ((t352 * t368 + t369 * t356) * t357 + (-t350 * t354 + t358 * t382) * t353) * t360) * qJD(4), 0; 0, t372, 0, -t340 * t353 + t342 * t357, (t346 * t356 + t366 * t352) * t360 + (-t340 * t357 - t342 * t353) * t363 + ((t349 * t356 + t370 * t352) * t363 + ((-t349 * t352 + t370 * t356) * t357 + (-t348 * t354 - t358 * t381) * t353) * t360) * qJD(4), 0; 0, 0, 0, t357 * t373 * t383 - t347 * t353, -t347 * t357 * t363 + ((t355 * t361 * t356 + (t355 * t379 + t384) * t352) * t363 + ((t356 * t384 + (-t352 * t361 + t356 * t379) * t355) * t357 + (-t355 * t364 * t354 + t359 * t358) * t353) * t360) * qJD(4) + ((-t352 * t380 + t356 * t364) * t360 - t353 * t363 * t383) * t373, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:11:03
	% EndTime: 2019-10-10 11:11:03
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (278->73), mult. (900->150), div. (0->0), fcn. (999->16), ass. (0->73)
	t453 = sin(pkin(8));
	t461 = sin(qJ(4));
	t496 = t453 * t461;
	t454 = sin(pkin(7));
	t459 = cos(pkin(6));
	t495 = t454 * t459;
	t455 = sin(pkin(6));
	t463 = sin(qJ(1));
	t494 = t455 * t463;
	t467 = cos(qJ(1));
	t493 = t455 * t467;
	t457 = cos(pkin(8));
	t492 = t457 * t461;
	t458 = cos(pkin(7));
	t462 = sin(qJ(2));
	t491 = t458 * t462;
	t466 = cos(qJ(2));
	t490 = t458 * t466;
	t489 = t463 * t462;
	t488 = t463 * t466;
	t487 = t467 * t462;
	t486 = t467 * t466;
	t485 = qJD(1) * t455;
	t484 = qJD(2) * t455;
	t460 = sin(qJ(5));
	t483 = qJD(4) * t460;
	t482 = t463 * t485;
	t481 = t467 * t485;
	t480 = t454 * t462 * t484;
	t479 = t453 * t480;
	t449 = t459 * t487 + t488;
	t452 = sin(pkin(14));
	t456 = cos(pkin(14));
	t448 = t459 * t486 - t489;
	t475 = t448 * t458 - t454 * t493;
	t430 = -t449 * t452 + t475 * t456;
	t443 = -t448 * t454 - t458 * t493;
	t478 = t430 * t457 + t443 * t453;
	t473 = t459 * t489 - t486;
	t450 = -t459 * t488 - t487;
	t474 = t450 * t458 + t454 * t494;
	t432 = t452 * t473 + t474 * t456;
	t444 = -t450 * t454 + t458 * t494;
	t477 = t432 * t457 + t444 * t453;
	t441 = t456 * t495 + (-t452 * t462 + t456 * t490) * t455;
	t447 = -t455 * t466 * t454 + t459 * t458;
	t476 = t441 * t457 + t447 * t453;
	t437 = -t448 * qJD(1) + t473 * qJD(2);
	t472 = t437 * t458 + t454 * t481;
	t439 = t450 * qJD(1) - t449 * qJD(2);
	t471 = t439 * t458 + t454 * t482;
	t431 = t449 * t456 + t475 * t452;
	t465 = cos(qJ(4));
	t470 = t431 * t465 + t478 * t461;
	t433 = t474 * t452 - t456 * t473;
	t469 = t433 * t465 + t477 * t461;
	t442 = t455 * t462 * t456 + (t455 * t490 + t495) * t452;
	t468 = t442 * t465 + t476 * t461;
	t464 = cos(qJ(5));
	t446 = (-t452 * t491 + t456 * t466) * t484;
	t445 = (-t452 * t466 - t456 * t491) * t484;
	t440 = -t473 * qJD(1) + t448 * qJD(2);
	t438 = -t449 * qJD(1) + t450 * qJD(2);
	t436 = -t445 * t453 + t457 * t480;
	t435 = -t439 * t454 + t458 * t482;
	t434 = -t437 * t454 + t458 * t481;
	t429 = t440 * t456 + t471 * t452;
	t428 = -t440 * t452 + t471 * t456;
	t427 = t438 * t456 + t472 * t452;
	t426 = -t438 * t452 + t472 * t456;
	t425 = -t428 * t453 + t435 * t457;
	t424 = -t426 * t453 + t434 * t457;
	t1 = [0, t481, 0, t424, t427 * t461 + (-t426 * t457 - t434 * t453) * t465 + t469 * qJD(4), (t426 * t492 + t427 * t465 + t434 * t496) * t460 - t424 * t464 + (t469 * t464 + (-t432 * t453 + t444 * t457) * t460) * qJD(5) + (-t433 * t461 + t477 * t465) * t483; 0, t482, 0, t425, t429 * t461 + (-t428 * t457 - t435 * t453) * t465 + t470 * qJD(4), (t428 * t492 + t429 * t465 + t435 * t496) * t460 - t425 * t464 + (t470 * t464 + (-t430 * t453 + t443 * t457) * t460) * qJD(5) + (-t431 * t461 + t478 * t465) * t483; 0, 0, 0, t436, t446 * t461 + (-t445 * t457 - t479) * t465 + t468 * qJD(4), (t445 * t492 + t446 * t465 + t461 * t479) * t460 - t436 * t464 + (t468 * t464 + (-t441 * t453 + t447 * t457) * t460) * qJD(5) + (-t442 * t461 + t476 * t465) * t483;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end