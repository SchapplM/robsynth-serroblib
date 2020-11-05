% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRRR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:49
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:29
	% EndTime: 2020-11-04 22:49:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:29
	% EndTime: 2020-11-04 22:49:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t233 = cos(qJ(1));
	t232 = sin(qJ(1));
	t1 = [t233, -t232, 0, 0; t232, t233, 0, 0; 0, 0, 1, pkin(9) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:29
	% EndTime: 2020-11-04 22:49:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t234 = sin(pkin(6));
	t237 = sin(qJ(1));
	t245 = t237 * t234;
	t236 = sin(qJ(2));
	t244 = t237 * t236;
	t238 = cos(qJ(2));
	t243 = t237 * t238;
	t239 = cos(qJ(1));
	t242 = t239 * t234;
	t241 = t239 * t236;
	t240 = t239 * t238;
	t235 = cos(pkin(6));
	t1 = [-t235 * t244 + t240, -t235 * t243 - t241, t245, t239 * pkin(1) + pkin(10) * t245 + 0; t235 * t241 + t243, t235 * t240 - t244, -t242, t237 * pkin(1) - pkin(10) * t242 + 0; t234 * t236, t234 * t238, t235, t235 * pkin(10) + pkin(9) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:29
	% EndTime: 2020-11-04 22:49:30
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (104->58), div. (0->0), fcn. (138->10), ass. (0->29)
	t252 = sin(pkin(7));
	t255 = cos(pkin(6));
	t274 = t252 * t255;
	t260 = cos(qJ(2));
	t273 = t252 * t260;
	t253 = sin(pkin(6));
	t254 = cos(pkin(7));
	t272 = t253 * t254;
	t271 = t255 * t260;
	t256 = sin(qJ(3));
	t257 = sin(qJ(2));
	t270 = t256 * t257;
	t269 = t256 * t260;
	t259 = cos(qJ(3));
	t268 = t257 * t259;
	t258 = sin(qJ(1));
	t267 = t258 * t257;
	t266 = t259 * t260;
	t261 = cos(qJ(1));
	t265 = t261 * t257;
	t264 = t261 * t260;
	t263 = -pkin(2) * t257 + pkin(11) * t273;
	t247 = -t252 * t253 + t254 * t271;
	t262 = t247 * t256 + t255 * t268;
	t251 = t254 * pkin(11) + pkin(10);
	t249 = t252 * t257 * pkin(11) + pkin(2) * t260 + pkin(1);
	t248 = -t254 * t270 + t266;
	t246 = t253 * t251 + t263 * t255;
	t1 = [t261 * t248 - t262 * t258, (-t247 * t258 - t254 * t265) * t259 - (-t255 * t267 + t264) * t256, (t258 * t271 + t265) * t252 + t258 * t272, t246 * t258 + t249 * t261 + 0; t258 * t248 + t262 * t261, (t247 * t259 - t255 * t270) * t261 - t258 * (t254 * t268 + t269), -(t255 * t264 - t267) * t252 - t261 * t272, -t246 * t261 + t249 * t258 + 0; t256 * t274 + (t254 * t269 + t268) * t253, t259 * t274 + (t254 * t266 - t270) * t253, -t253 * t273 + t255 * t254, t251 * t255 - t263 * t253 + pkin(9) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:30
	% EndTime: 2020-11-04 22:49:30
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (133->51), mult. (336->93), div. (0->0), fcn. (425->14), ass. (0->46)
	t296 = sin(pkin(8));
	t299 = cos(pkin(8));
	t300 = cos(pkin(7));
	t297 = sin(pkin(7));
	t307 = cos(qJ(3));
	t321 = t297 * t307;
	t328 = t296 * t300 + t299 * t321;
	t295 = t299 * pkin(12) + pkin(11);
	t303 = sin(qJ(3));
	t326 = pkin(12) * t296;
	t312 = -pkin(3) * t303 + t307 * t326;
	t285 = t297 * t295 + t312 * t300;
	t289 = pkin(3) * t307 + t303 * t326 + pkin(2);
	t304 = sin(qJ(2));
	t308 = cos(qJ(2));
	t327 = -t285 * t308 + t289 * t304;
	t322 = t297 * t303;
	t320 = t299 * t300;
	t319 = t300 * t307;
	t301 = cos(pkin(6));
	t318 = t301 * t308;
	t317 = t303 * t304;
	t316 = t303 * t308;
	t314 = t301 * t317;
	t287 = -t297 * t296 + t299 * t319;
	t298 = sin(pkin(6));
	t278 = t287 * t318 - t328 * t298 - t299 * t314;
	t311 = t300 * t316 + t304 * t307;
	t281 = -t298 * t322 + t311 * t301;
	t302 = sin(qJ(4));
	t306 = cos(qJ(4));
	t313 = t278 * t302 + t306 * t281;
	t310 = t295 * t300 - t312 * t297 + pkin(10);
	t309 = cos(qJ(1));
	t305 = sin(qJ(1));
	t288 = t300 * t317 - t307 * t308;
	t286 = t296 * t319 + t297 * t299;
	t284 = t304 * t287 + t299 * t316;
	t283 = t304 * t286 + t296 * t316;
	t282 = t311 * t298 + t301 * t322;
	t280 = t285 * t304 + t289 * t308 + pkin(1);
	t279 = t284 * t302 + t306 * t288;
	t277 = t286 * t318 + t298 * t320 + (-t298 * t321 - t314) * t296;
	t276 = t328 * t301 + (t287 * t308 - t299 * t317) * t298;
	t275 = t298 * t310 - t327 * t301;
	t1 = [-t279 * t309 - t313 * t305, (-t278 * t305 - t309 * t284) * t306 + (t281 * t305 + t309 * t288) * t302, t277 * t305 + t309 * t283, t275 * t305 + t280 * t309 + 0; -t279 * t305 + t313 * t309, (t278 * t306 - t302 * t281) * t309 - (t284 * t306 - t302 * t288) * t305, -t277 * t309 + t305 * t283, -t275 * t309 + t280 * t305 + 0; t276 * t302 + t282 * t306, t276 * t306 - t282 * t302, -t298 * t286 * t308 + t301 * t320 + (t298 * t317 - t301 * t321) * t296, t327 * t298 + t310 * t301 + pkin(9) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:30
	% EndTime: 2020-11-04 22:49:30
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (256->78), mult. (642->143), div. (0->0), fcn. (786->16), ass. (0->65)
	t359 = cos(pkin(7));
	t362 = sin(qJ(4));
	t363 = sin(qJ(3));
	t367 = cos(qJ(4));
	t384 = t367 * t363;
	t358 = cos(pkin(8));
	t368 = cos(qJ(3));
	t391 = t358 * t368;
	t355 = sin(pkin(8));
	t356 = sin(pkin(7));
	t395 = t356 * t355;
	t341 = -t362 * t395 + (t362 * t391 + t384) * t359;
	t357 = sin(pkin(6));
	t360 = cos(pkin(6));
	t369 = cos(qJ(2));
	t364 = sin(qJ(2));
	t387 = t363 * t364;
	t381 = t358 * t387;
	t382 = t356 * t384;
	t392 = t358 * t362;
	t383 = t356 * t392;
	t385 = t364 * t367;
	t396 = t355 * t359;
	t332 = (-t341 * t369 + t362 * t381) * t360 + (t357 * t383 - t360 * t385) * t368 + t357 * (t362 * t396 + t382);
	t389 = t359 * t368;
	t349 = t355 * t389 + t356 * t358;
	t394 = t356 * t368;
	t336 = t360 * t349 * t369 + t358 * t357 * t359 + (-t357 * t394 - t360 * t387) * t355;
	t361 = sin(qJ(5));
	t366 = cos(qJ(5));
	t408 = t332 * t361 - t336 * t366;
	t404 = pkin(13) * t367;
	t380 = pkin(4) * t362 - t404;
	t407 = t358 * pkin(12) + t380 * t355 + pkin(11);
	t403 = t355 * pkin(12);
	t347 = -t380 * t358 + t403;
	t352 = t367 * pkin(4) + t362 * pkin(13) + pkin(3);
	t398 = t352 * t363;
	t378 = t347 * t368 - t398;
	t334 = t356 * t407 + t378 * t359;
	t339 = t347 * t363 + t352 * t368 + pkin(2);
	t406 = t334 * t369 - t339 * t364;
	t393 = t357 * t369;
	t390 = t359 * t360;
	t388 = t363 * t362;
	t386 = t363 * t369;
	t350 = t358 * t389 - t395;
	t376 = -t350 * t369 + t381;
	t375 = -t356 * t391 - t396;
	t374 = t359 * t386 + t364 * t368;
	t372 = t407 * t359 + pkin(10);
	t371 = t341 * t393 - (-t355 * t390 + t357 * t381) * t362 + (t357 * t385 + t360 * t383) * t368 + t360 * t382;
	t370 = cos(qJ(1));
	t365 = sin(qJ(1));
	t351 = t359 * t387 - t368 * t369;
	t343 = t364 * t350 + t358 * t386;
	t342 = t364 * t349 + t355 * t386;
	t340 = -t357 * t356 * t363 + t374 * t360;
	t338 = (t358 * t388 - t368 * t367) * t369 + t364 * t341;
	t337 = t375 * t357 - t376 * t360;
	t335 = -t349 * t393 + t358 * t390 + (t357 * t387 - t360 * t394) * t355;
	t333 = t338 * t361 + t342 * t366;
	t330 = t334 * t364 + t339 * t369 + pkin(1);
	t329 = t357 * (-t378 * t356 + t372) + t406 * t360;
	t1 = [(t332 * t365 - t338 * t370) * t366 + (t336 * t365 + t370 * t342) * t361, t370 * t333 - t408 * t365, (t337 * t365 + t370 * t343) * t367 - (t340 * t365 + t370 * t351) * t362, t329 * t365 + t330 * t370 + 0; (-t332 * t366 - t361 * t336) * t370 - t365 * (t338 * t366 - t342 * t361), t365 * t333 + t408 * t370, (-t337 * t367 + t362 * t340) * t370 + (t343 * t367 - t362 * t351) * t365, -t329 * t370 + t330 * t365 + 0; t335 * t361 + t371 * t366, t335 * t366 - t371 * t361, (t356 * t388 + t375 * t367) * t360 + (t374 * t362 + t376 * t367) * t357, pkin(9) + 0 - t406 * t357 + (((pkin(4) * t392 - t358 * t404 - t403) * t368 + t398) * t356 + t372) * t360; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:49:30
	% EndTime: 2020-11-04 22:49:31
	% DurationCPUTime: 0.85s
	% Computational Cost: add. (455->116), mult. (1081->206), div. (0->0), fcn. (1314->18), ass. (0->85)
	t449 = sin(pkin(8));
	t452 = cos(pkin(8));
	t456 = sin(qJ(5));
	t462 = cos(qJ(5));
	t520 = t456 * pkin(5) - pkin(14) * t462;
	t477 = -pkin(12) - t520;
	t457 = sin(qJ(4));
	t463 = cos(qJ(4));
	t513 = t463 * pkin(13);
	t521 = pkin(5) * t462 + pkin(14) * t456 + pkin(4);
	t522 = t457 * t521 - t513;
	t526 = t477 * t449 + t522 * t452;
	t524 = t463 * t521 + pkin(3);
	t435 = t457 * pkin(13) + t524;
	t497 = t457 * t462;
	t438 = -t449 * t456 + t452 * t497;
	t450 = sin(pkin(7));
	t453 = cos(pkin(7));
	t464 = cos(qJ(3));
	t458 = sin(qJ(3));
	t492 = t463 * t458;
	t421 = -t450 * (t449 * t497 + t452 * t456) + (t438 * t464 + t462 * t492) * t453;
	t454 = cos(pkin(6));
	t459 = sin(qJ(2));
	t495 = t458 * t459;
	t480 = t454 * t495;
	t451 = sin(pkin(6));
	t506 = t451 * t453;
	t432 = t449 * t506 + t452 * t480;
	t494 = t459 * t463;
	t479 = t462 * t494;
	t483 = t452 * t506;
	t508 = t450 * t451;
	t485 = t458 * t508;
	t465 = cos(qJ(2));
	t499 = t454 * t465;
	t412 = (t432 * t457 + t463 * t485) * t462 - t421 * t499 - (-t438 * t508 + t454 * t479) * t464 - t456 * (t449 * t480 - t483);
	t490 = t464 * t463;
	t496 = t458 * t457;
	t509 = t449 * t463;
	t428 = -t450 * t509 + (t452 * t490 - t496) * t453;
	t507 = t450 * t452;
	t488 = t451 * t507;
	t498 = t457 * t459;
	t416 = t428 * t499 + (-t454 * t498 - t463 * t488) * t464 - t432 * t463 + t457 * t485;
	t455 = sin(qJ(6));
	t461 = cos(qJ(6));
	t525 = t412 * t455 - t461 * t416;
	t491 = t464 * t457;
	t510 = t449 * t457;
	t427 = -t450 * t510 + (t452 * t491 + t492) * t453;
	t501 = t453 * t464;
	t437 = t449 * t501 + t507;
	t482 = t452 * t495;
	t523 = ((-t427 * t465 + t457 * t482) * t454 + (-t454 * t494 + t457 * t488) * t464 + t451 * (t450 * t492 + t453 * t510)) * t456 - (t437 * t499 + t483 + (-t464 * t508 - t480) * t449) * t462;
	t503 = t453 * t454;
	t442 = t452 * t503;
	t505 = t451 * t459;
	t486 = t449 * t505;
	t430 = t458 * t486 + t442;
	t489 = t449 * t503;
	t433 = -t451 * t482 + t489;
	t500 = t454 * t450;
	t481 = t463 * t500;
	t470 = t433 * t457 + t458 * t481;
	t504 = t451 * t465;
	t518 = t470 * t462 + t421 * t504 + (t438 * t500 + t451 * t479) * t464 + t456 * t430;
	t502 = t453 * t458;
	t487 = t457 * t500;
	t484 = t451 * t494;
	t418 = t421 * t459 - (-t458 * t438 + t462 * t490) * t465;
	t469 = pkin(13) * t509 + t477 * t452 - pkin(11);
	t466 = cos(qJ(1));
	t460 = sin(qJ(1));
	t445 = t452 * pkin(12) + pkin(11);
	t423 = t526 * t501;
	t422 = (t452 * t492 + t491) * t465 + t459 * t428;
	t419 = t435 * t464 - t458 * t526 + pkin(2);
	t417 = ((t452 * t496 - t490) * t465 + t459 * t427) * t456 + (t449 * t458 * t465 + t459 * t437) * t462;
	t415 = t428 * t504 + (-t451 * t498 + t452 * t481) * t464 + t433 * t463 - t458 * t487;
	t414 = -t423 - t435 * t502 - t450 * (-t510 * t521 + t469);
	t413 = t418 * t455 + t422 * t461;
	t410 = t414 * t459 + t419 * t465 + pkin(1);
	t409 = -(-pkin(10) + (-t435 * t458 - t526 * t464) * t450) * t451 + (t414 * t465 - t419 * t459) * t454 - (-t449 * t522 - t452 * t520 - t445) * t506;
	t1 = [(t412 * t460 - t466 * t418) * t461 + (t416 * t460 + t422 * t466) * t455, t466 * t413 - t525 * t460, -t466 * t417 + t523 * t460, t409 * t460 + t410 * t466 + 0; (-t412 * t461 - t416 * t455) * t466 + t460 * (-t418 * t461 + t422 * t455), t413 * t460 + t525 * t466, -t460 * t417 - t523 * t466, -t409 * t466 + t410 * t460 + 0; -t415 * t455 + t518 * t461, -t415 * t461 - t518 * t455, (t427 * t504 + (t452 * t487 + t484) * t464 + t470) * t456 - (-t449 * t464 * t500 - t437 * t504 + t430) * t462, -(-t423 + (t449 * t450 * t521 - pkin(13) * t502) * t457 - t524 * t502 - t450 * t469) * t504 + (t435 * t505 + t526 * t500) * t464 + ((-t452 * t505 * t521 + pkin(13) * t500) * t458 + t521 * t489) * t457 + (t452 * pkin(13) * t484 - t477 * t486 + t500 * t524) * t458 - t489 * t513 + pkin(2) * t505 + (t445 * t453 + pkin(10)) * t454 + 0 + pkin(9) + t520 * t442; 0, 0, 0, 1;];
	Tc_mdh = t1;
end