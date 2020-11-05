% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:06
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:02
	% EndTime: 2020-11-04 21:06:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:02
	% EndTime: 2020-11-04 21:06:02
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t196 = cos(pkin(13));
	t195 = sin(pkin(13));
	t1 = [t196, -t195, 0, 0; t195, t196, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:02
	% EndTime: 2020-11-04 21:06:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t197 = sin(pkin(13));
	t198 = sin(pkin(6));
	t206 = t197 * t198;
	t199 = cos(pkin(13));
	t205 = t199 * t198;
	t200 = cos(pkin(6));
	t201 = sin(qJ(2));
	t204 = t200 * t201;
	t202 = cos(qJ(2));
	t203 = t200 * t202;
	t1 = [-t197 * t204 + t199 * t202, -t197 * t203 - t199 * t201, t206, t199 * pkin(1) + pkin(9) * t206 + 0; t197 * t202 + t199 * t204, -t197 * t201 + t199 * t203, -t205, t197 * pkin(1) - pkin(9) * t205 + 0; t198 * t201, t198 * t202, t200, t200 * pkin(9) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:02
	% EndTime: 2020-11-04 21:06:02
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (44->42), mult. (116->82), div. (0->0), fcn. (150->10), ass. (0->32)
	t210 = sin(pkin(7));
	t237 = qJ(3) * t210;
	t209 = sin(pkin(13));
	t211 = sin(pkin(6));
	t236 = t209 * t211;
	t214 = cos(pkin(7));
	t235 = t209 * t214;
	t215 = cos(pkin(6));
	t234 = t209 * t215;
	t217 = cos(qJ(2));
	t233 = t210 * t217;
	t213 = cos(pkin(13));
	t232 = t211 * t213;
	t212 = cos(pkin(14));
	t231 = t212 * t209;
	t230 = t212 * t215;
	t229 = t213 * t214;
	t228 = t213 * t215;
	t227 = t214 * t217;
	t226 = t215 * t214;
	t225 = t215 * t217;
	t208 = sin(pkin(14));
	t224 = t208 * t234;
	t223 = t210 * t236;
	t222 = t209 * t237;
	t221 = t209 * t230;
	t220 = t213 * t237;
	t219 = t210 * t232;
	t218 = t213 * t226;
	t216 = sin(qJ(2));
	t207 = t214 * qJ(3) + pkin(9);
	t1 = [(t212 * t213 - t214 * t224) * t217 + (-t208 * t229 - t221) * t216 + t208 * t223, (-t208 * t213 - t214 * t221) * t217 + (-t212 * t229 + t224) * t216 + t212 * t223, (t209 * t225 + t213 * t216) * t210 + t211 * t235, (t213 * pkin(2) + t215 * t222) * t217 + (-pkin(2) * t234 + t220) * t216 + t207 * t236 + t213 * pkin(1) + 0; (t208 * t218 + t231) * t217 + (-t208 * t235 + t212 * t228) * t216 - t208 * t219, (-t209 * t208 + t212 * t218) * t217 + (-t208 * t228 - t214 * t231) * t216 - t212 * t219, -(-t209 * t216 + t213 * t225) * t210 - t211 * t229, (t209 * pkin(2) - t215 * t220) * t217 + (pkin(2) * t228 + t222) * t216 - t207 * t232 + t209 * pkin(1) + 0; t211 * t216 * t212 + (t210 * t215 + t211 * t227) * t208, t210 * t230 + (-t208 * t216 + t212 * t227) * t211, -t211 * t233 + t226, t207 * t215 + qJ(1) + 0 + (pkin(2) * t216 - qJ(3) * t233) * t211; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:02
	% EndTime: 2020-11-04 21:06:02
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (133->65), mult. (331->121), div. (0->0), fcn. (420->14), ass. (0->50)
	t254 = sin(pkin(14));
	t255 = sin(pkin(13));
	t287 = t254 * t255;
	t256 = sin(pkin(8));
	t286 = t254 * t256;
	t265 = sin(qJ(2));
	t285 = t254 * t265;
	t258 = sin(pkin(6));
	t284 = t255 * t258;
	t263 = cos(pkin(6));
	t283 = t255 * t263;
	t262 = cos(pkin(7));
	t282 = t256 * t262;
	t257 = sin(pkin(7));
	t281 = t257 * t256;
	t261 = cos(pkin(8));
	t280 = t257 * t261;
	t260 = cos(pkin(13));
	t279 = t258 * t260;
	t259 = cos(pkin(14));
	t278 = t259 * t255;
	t277 = t260 * t262;
	t276 = t260 * t263;
	t275 = t261 * t260;
	t274 = t262 * t261;
	t273 = t254 * t277;
	t272 = t254 * t283;
	t271 = t258 * t257 * t254;
	t270 = t263 * t278;
	t269 = t261 * t283;
	t268 = t254 * t275;
	t267 = cos(qJ(2));
	t266 = cos(qJ(4));
	t264 = sin(qJ(4));
	t253 = t261 * pkin(10) + qJ(3);
	t251 = t259 * t256 * pkin(10) - t254 * pkin(3);
	t250 = t259 * pkin(3) + pkin(10) * t286 + pkin(2);
	t249 = t259 * t274 - t281;
	t248 = t259 * t280 + t282;
	t247 = t259 * t281 - t274;
	t246 = t259 * t282 + t280;
	t245 = t251 * t262 + t257 * t253;
	t244 = -t251 * t257 + t253 * t262 + pkin(9);
	t243 = t258 * t265 * t259 + (t258 * t262 * t267 + t257 * t263) * t254;
	t242 = t263 * t248 + (t249 * t267 - t261 * t285) * t258;
	t241 = (-t259 * t260 + t262 * t272) * t267 + (t270 + t273) * t265 - t255 * t271;
	t240 = (t263 * t273 + t278) * t267 + (t259 * t276 - t262 * t287) * t265 - t260 * t271;
	t239 = (t249 * t276 - t261 * t287) * t267 + (-t255 * t249 - t263 * t268) * t265 - t248 * t279;
	t238 = (-t249 * t283 - t268) * t267 + (-t260 * t249 + t254 * t269) * t265 + t248 * t284;
	t1 = [t238 * t264 - t266 * t241, t238 * t266 + t264 * t241, ((t254 * t260 + t262 * t270) * t256 + t257 * t269) * t267 + ((t259 * t277 - t272) * t256 + t257 * t275) * t265 - t247 * t284, (t245 * t283 + t260 * t250) * t267 + (t245 * t260 - t250 * t283) * t265 + t244 * t284 + t260 * pkin(1) + 0; t239 * t264 + t240 * t266, t239 * t266 - t240 * t264, (-t246 * t276 + t255 * t286) * t267 + (t255 * t246 + t276 * t286) * t265 + t247 * t279, (-t245 * t276 + t255 * t250) * t267 + (t245 * t255 + t250 * t276) * t265 - t244 * t279 + t255 * pkin(1) + 0; t242 * t264 + t243 * t266, t242 * t266 - t243 * t264, -t263 * t247 + (-t246 * t267 + t256 * t285) * t258, t244 * t263 + qJ(1) + 0 + (-t245 * t267 + t250 * t265) * t258; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:02
	% EndTime: 2020-11-04 21:06:03
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (256->100), mult. (671->177), div. (0->0), fcn. (815->16), ass. (0->76)
	t339 = sin(qJ(2));
	t342 = cos(qJ(2));
	t332 = cos(pkin(14));
	t327 = sin(pkin(14));
	t329 = sin(pkin(8));
	t379 = t327 * t329;
	t315 = t332 * pkin(3) + pkin(10) * t379 + pkin(2);
	t334 = cos(pkin(8));
	t377 = t327 * t334;
	t317 = pkin(4) * t332 + pkin(11) * t377;
	t318 = pkin(4) * t377 - t332 * pkin(11);
	t338 = sin(qJ(4));
	t341 = cos(qJ(4));
	t344 = t317 * t341 - t318 * t338 + t315;
	t316 = t332 * t329 * pkin(10) - t327 * pkin(3);
	t326 = t334 * pkin(10) + qJ(3);
	t330 = sin(pkin(7));
	t335 = cos(pkin(7));
	t299 = t316 * t335 + t330 * t326;
	t366 = t332 * t334;
	t319 = -t327 * pkin(4) + pkin(11) * t366;
	t370 = t330 * t329;
	t302 = -pkin(11) * t370 + t319 * t335;
	t320 = pkin(4) * t366 + t327 * pkin(11);
	t347 = pkin(4) * t370 - t320 * t335;
	t345 = t302 * t341 + t347 * t338 + t299;
	t388 = t344 * t339 - t345 * t342;
	t371 = t329 * t335;
	t346 = t316 * t330 - t326 * t335 - pkin(9) + (pkin(11) * t371 + t330 * t319) * t341 - (pkin(4) * t371 + t330 * t320) * t338;
	t364 = t335 * t334;
	t313 = t332 * t364 - t370;
	t328 = sin(pkin(13));
	t336 = cos(pkin(6));
	t333 = cos(pkin(13));
	t353 = t333 * t377;
	t292 = t328 * t313 + t336 * t353;
	t354 = t328 * t377;
	t365 = t333 * t336;
	t294 = t313 * t365 - t354;
	t312 = t330 * t366 + t371;
	t374 = t327 * t341;
	t296 = t312 * t338 + t330 * t374;
	t376 = t327 * t335;
	t358 = t333 * t376;
	t367 = t332 * t328;
	t305 = t336 * t358 + t367;
	t357 = t328 * t376;
	t308 = t332 * t365 - t357;
	t331 = sin(pkin(6));
	t369 = t331 * t333;
	t384 = (t292 * t338 - t341 * t308) * t339 - (t294 * t338 + t341 * t305) * t342 + t296 * t369;
	t291 = -t333 * t313 + t336 * t354;
	t372 = t328 * t336;
	t293 = t313 * t372 + t353;
	t306 = -t332 * t333 + t336 * t357;
	t307 = t336 * t367 + t358;
	t373 = t328 * t331;
	t383 = (t291 * t338 - t341 * t307) * t339 - (t293 * t338 + t341 * t306) * t342 + t296 * t373;
	t378 = t327 * t330;
	t375 = t327 * t339;
	t368 = t331 * t342;
	t363 = t336 * t312;
	t309 = t331 * t339 * t332 + t336 * t378;
	t360 = t331 * t378;
	t359 = t334 * t375;
	t356 = t328 * t379;
	t355 = t333 * t379;
	t343 = -(t331 * t359 - t363) * t338 + (t313 * t338 + t335 * t374) * t368 + t341 * t309;
	t340 = cos(qJ(5));
	t337 = sin(qJ(5));
	t311 = -t332 * t370 + t364;
	t310 = t330 * t334 + t332 * t371;
	t290 = t336 * t311 + (-t310 * t342 + t329 * t375) * t331;
	t289 = (t310 * t365 - t356) * t342 + (-t328 * t310 - t336 * t355) * t339 + t311 * t369;
	t288 = (t310 * t372 + t355) * t342 + (t333 * t310 - t336 * t356) * t339 + t311 * t373;
	t1 = [t288 * t337 + t383 * t340, t288 * t340 - t383 * t337, (-t291 * t339 + t293 * t342 - t312 * t373) * t341 - t338 * (t306 * t342 + t307 * t339 - t328 * t360), 0 + (t345 * t339 + t344 * t342 + pkin(1)) * t333 + (-t346 * t331 - t388 * t336) * t328; -t289 * t337 - t384 * t340, -t340 * t289 + t384 * t337, (t292 * t339 - t294 * t342 + t312 * t369) * t341 + (t305 * t342 + t308 * t339 - t333 * t360) * t338, ((-t302 * t365 + t328 * t317) * t341 + (-t328 * t318 - t347 * t365) * t338 - t299 * t365 + t328 * t315) * t342 + ((t328 * t302 + t317 * t365) * t341 + (-t318 * t365 + t328 * t347) * t338 + t315 * t365 + t299 * t328) * t339 + t328 * pkin(1) + 0 + t346 * t369; t337 * t290 + t343 * t340, t290 * t340 - t343 * t337, (-t363 + (-t313 * t342 + t359) * t331) * t341 + (t368 * t376 + t309) * t338, t388 * t331 - t346 * t336 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:06:03
	% EndTime: 2020-11-04 21:06:04
	% DurationCPUTime: 0.63s
	% Computational Cost: add. (455->115), mult. (1241->195), div. (0->0), fcn. (1578->18), ass. (0->87)
	t450 = sin(qJ(2));
	t454 = cos(qJ(2));
	t442 = cos(pkin(14));
	t437 = sin(pkin(14));
	t439 = sin(pkin(8));
	t491 = t437 * t439;
	t425 = t442 * pkin(3) + pkin(10) * t491 + pkin(2);
	t444 = cos(pkin(8));
	t489 = t437 * t444;
	t427 = pkin(4) * t442 + pkin(11) * t489;
	t428 = pkin(4) * t489 - t442 * pkin(11);
	t449 = sin(qJ(4));
	t453 = cos(qJ(4));
	t456 = t427 * t453 - t428 * t449 + t425;
	t426 = t442 * t439 * pkin(10) - t437 * pkin(3);
	t436 = t444 * pkin(10) + qJ(3);
	t440 = sin(pkin(7));
	t445 = cos(pkin(7));
	t409 = t426 * t445 + t440 * t436;
	t478 = t442 * t444;
	t429 = -t437 * pkin(4) + pkin(11) * t478;
	t482 = t440 * t439;
	t412 = -pkin(11) * t482 + t429 * t445;
	t430 = pkin(4) * t478 + t437 * pkin(11);
	t459 = pkin(4) * t482 - t430 * t445;
	t457 = t412 * t453 + t459 * t449 + t409;
	t500 = t456 * t450 - t457 * t454;
	t483 = t439 * t445;
	t458 = t426 * t440 - t436 * t445 - pkin(9) + (pkin(11) * t483 + t440 * t429) * t453 - (pkin(4) * t483 + t440 * t430) * t449;
	t476 = t445 * t444;
	t423 = t442 * t476 - t482;
	t438 = sin(pkin(13));
	t446 = cos(pkin(6));
	t443 = cos(pkin(13));
	t465 = t443 * t489;
	t402 = t438 * t423 + t446 * t465;
	t466 = t438 * t489;
	t477 = t443 * t446;
	t404 = t423 * t477 - t466;
	t422 = t440 * t478 + t483;
	t486 = t437 * t453;
	t406 = t422 * t449 + t440 * t486;
	t488 = t437 * t445;
	t470 = t443 * t488;
	t479 = t442 * t438;
	t415 = t446 * t470 + t479;
	t469 = t438 * t488;
	t418 = t442 * t477 - t469;
	t441 = sin(pkin(6));
	t481 = t441 * t443;
	t496 = (t402 * t449 - t453 * t418) * t450 - (t404 * t449 + t453 * t415) * t454 + t406 * t481;
	t401 = -t443 * t423 + t446 * t466;
	t484 = t438 * t446;
	t403 = t423 * t484 + t465;
	t416 = -t442 * t443 + t446 * t469;
	t417 = t446 * t479 + t470;
	t485 = t438 * t441;
	t495 = (t401 * t449 - t453 * t417) * t450 - (t403 * t449 + t453 * t416) * t454 + t406 * t485;
	t490 = t437 * t440;
	t487 = t437 * t450;
	t480 = t441 * t454;
	t475 = t446 * t422;
	t419 = t441 * t450 * t442 + t446 * t490;
	t472 = t441 * t490;
	t471 = t444 * t487;
	t468 = t438 * t491;
	t467 = t443 * t491;
	t455 = -(t441 * t471 - t475) * t449 + (t423 * t449 + t445 * t486) * t480 + t453 * t419;
	t452 = cos(qJ(5));
	t451 = cos(qJ(6));
	t448 = sin(qJ(5));
	t447 = sin(qJ(6));
	t421 = -t442 * t482 + t476;
	t420 = t440 * t444 + t442 * t483;
	t400 = t446 * t421 + (-t420 * t454 + t439 * t487) * t441;
	t399 = (-t475 + (-t423 * t454 + t471) * t441) * t453 + (t480 * t488 + t419) * t449;
	t398 = (t420 * t477 - t468) * t454 + (-t438 * t420 - t446 * t467) * t450 + t421 * t481;
	t397 = (t420 * t484 + t467) * t454 + (t443 * t420 - t446 * t468) * t450 + t421 * t485;
	t396 = t400 * t452 - t455 * t448;
	t395 = t448 * t400 + t455 * t452;
	t394 = (t402 * t450 - t404 * t454 + t422 * t481) * t453 + (t415 * t454 + t418 * t450 - t443 * t472) * t449;
	t393 = (-t401 * t450 + t403 * t454 - t422 * t485) * t453 - t449 * (t416 * t454 + t417 * t450 - t438 * t472);
	t392 = -t452 * t398 + t496 * t448;
	t391 = -t398 * t448 - t496 * t452;
	t390 = t397 * t448 + t495 * t452;
	t389 = t397 * t452 - t495 * t448;
	t1 = [t390 * t451 + t393 * t447, -t390 * t447 + t393 * t451, -t389, t390 * pkin(5) - t389 * pkin(12) + 0 + (t457 * t450 + t456 * t454 + pkin(1)) * t443 + (-t458 * t441 - t500 * t446) * t438; t391 * t451 + t394 * t447, -t391 * t447 + t394 * t451, -t392, t391 * pkin(5) - t392 * pkin(12) + ((-t412 * t477 + t438 * t427) * t453 + (-t438 * t428 - t459 * t477) * t449 - t409 * t477 + t438 * t425) * t454 + ((t438 * t412 + t427 * t477) * t453 + (-t428 * t477 + t438 * t459) * t449 + t425 * t477 + t409 * t438) * t450 + t438 * pkin(1) + 0 + t458 * t481; t395 * t451 + t399 * t447, -t395 * t447 + t399 * t451, -t396, t395 * pkin(5) - t396 * pkin(12) + t500 * t441 - t458 * t446 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end