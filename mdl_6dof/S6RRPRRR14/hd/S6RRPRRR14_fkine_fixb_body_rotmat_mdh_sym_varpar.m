% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR14 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:21
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRR14_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:14
	% EndTime: 2020-11-04 22:21:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:14
	% EndTime: 2020-11-04 22:21:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t209 = cos(qJ(1));
	t208 = sin(qJ(1));
	t1 = [t209, -t208, 0, 0; t208, t209, 0, 0; 0, 0, 1, pkin(9) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:14
	% EndTime: 2020-11-04 22:21:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t210 = sin(pkin(6));
	t213 = sin(qJ(1));
	t221 = t213 * t210;
	t212 = sin(qJ(2));
	t220 = t213 * t212;
	t214 = cos(qJ(2));
	t219 = t213 * t214;
	t215 = cos(qJ(1));
	t218 = t215 * t210;
	t217 = t215 * t212;
	t216 = t215 * t214;
	t211 = cos(pkin(6));
	t1 = [-t211 * t220 + t216, -t211 * t219 - t217, t221, t215 * pkin(1) + pkin(10) * t221 + 0; t211 * t217 + t219, t211 * t216 - t220, -t218, t213 * pkin(1) - pkin(10) * t218 + 0; t210 * t212, t210 * t214, t211, t211 * pkin(10) + pkin(9) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:14
	% EndTime: 2020-11-04 22:21:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (44->30), mult. (108->57), div. (0->0), fcn. (142->10), ass. (0->27)
	t230 = sin(pkin(7));
	t234 = cos(pkin(6));
	t247 = t230 * t234;
	t237 = cos(qJ(2));
	t246 = t230 * t237;
	t231 = sin(pkin(6));
	t233 = cos(pkin(7));
	t245 = t231 * t233;
	t235 = sin(qJ(2));
	t244 = t233 * t235;
	t243 = t233 * t237;
	t242 = t234 * t235;
	t241 = t234 * t237;
	t240 = -pkin(2) * t235 + qJ(3) * t246;
	t239 = -t230 * t231 + t233 * t241;
	t238 = cos(qJ(1));
	t236 = sin(qJ(1));
	t232 = cos(pkin(14));
	t229 = sin(pkin(14));
	t228 = t233 * qJ(3) + pkin(10);
	t227 = t230 * t235 * qJ(3) + pkin(2) * t237 + pkin(1);
	t226 = t229 * t237 + t232 * t244;
	t225 = t229 * t244 - t232 * t237;
	t224 = t231 * t228 + t240 * t234;
	t223 = -t229 * t242 + t239 * t232;
	t222 = t239 * t229 + t232 * t242;
	t1 = [-t222 * t236 - t238 * t225, -t223 * t236 - t238 * t226, (t238 * t235 + t236 * t241) * t230 + t236 * t245, t224 * t236 + t227 * t238 + 0; t222 * t238 - t236 * t225, t223 * t238 - t236 * t226, -(-t236 * t235 + t238 * t241) * t230 - t238 * t245, -t224 * t238 + t227 * t236 + 0; t231 * t235 * t232 + (t231 * t243 + t247) * t229, t232 * t247 + (-t229 * t235 + t232 * t243) * t231, -t231 * t246 + t234 * t233, t228 * t234 - t240 * t231 + pkin(9) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:14
	% EndTime: 2020-11-04 22:21:14
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (133->50), mult. (316->90), div. (0->0), fcn. (405->14), ass. (0->48)
	t270 = sin(pkin(8));
	t296 = pkin(11) * t270;
	t269 = sin(pkin(14));
	t278 = sin(qJ(2));
	t295 = t269 * t278;
	t281 = cos(qJ(2));
	t294 = t269 * t281;
	t275 = cos(pkin(7));
	t293 = t270 * t275;
	t271 = sin(pkin(7));
	t292 = t271 * t270;
	t274 = cos(pkin(8));
	t291 = t271 * t274;
	t273 = cos(pkin(14));
	t290 = t273 * t278;
	t289 = t275 * t274;
	t276 = cos(pkin(6));
	t288 = t276 * t281;
	t287 = t274 * t295;
	t263 = t273 * t289 - t292;
	t272 = sin(pkin(6));
	t283 = t273 * t291 + t293;
	t251 = t263 * t288 - t272 * t283 - t276 * t287;
	t254 = t276 * t290 + (-t271 * t272 + t275 * t288) * t269;
	t277 = sin(qJ(4));
	t280 = cos(qJ(4));
	t286 = t251 * t277 + t280 * t254;
	t266 = -t269 * pkin(3) + t273 * t296;
	t268 = t274 * pkin(11) + qJ(3);
	t259 = t266 * t275 + t271 * t268;
	t265 = t273 * pkin(3) + t269 * t296 + pkin(2);
	t285 = t259 * t281 - t265 * t278;
	t261 = t273 * t293 + t291;
	t284 = -t261 * t281 + t270 * t295;
	t282 = cos(qJ(1));
	t279 = sin(qJ(1));
	t264 = -t273 * t281 + t275 * t295;
	t262 = t273 * t292 - t289;
	t258 = -t266 * t271 + t268 * t275 + pkin(10);
	t257 = t278 * t263 + t274 * t294;
	t256 = t278 * t261 + t270 * t294;
	t255 = t272 * t290 + (t272 * t275 * t281 + t271 * t276) * t269;
	t253 = t259 * t278 + t265 * t281 + pkin(1);
	t252 = t257 * t277 + t280 * t264;
	t250 = -t272 * t262 - t284 * t276;
	t249 = t276 * t283 + (t263 * t281 - t287) * t272;
	t248 = t258 * t272 + t285 * t276;
	t1 = [-t282 * t252 - t286 * t279, (-t251 * t279 - t257 * t282) * t280 + t277 * (t254 * t279 + t282 * t264), t250 * t279 + t256 * t282, t248 * t279 + t253 * t282 + 0; -t279 * t252 + t286 * t282, (t251 * t280 - t277 * t254) * t282 - (t257 * t280 - t277 * t264) * t279, -t250 * t282 + t256 * t279, -t248 * t282 + t253 * t279 + 0; t249 * t277 + t255 * t280, t249 * t280 - t255 * t277, -t276 * t262 + t284 * t272, t258 * t276 - t285 * t272 + pkin(9) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:14
	% EndTime: 2020-11-04 22:21:14
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (256->77), mult. (623->142), div. (0->0), fcn. (767->16), ass. (0->62)
	t333 = sin(pkin(14));
	t334 = sin(pkin(8));
	t337 = cos(pkin(14));
	t364 = t334 * t337;
	t322 = -pkin(3) * t333 + pkin(11) * t364;
	t338 = cos(pkin(8));
	t361 = t337 * t338;
	t325 = -pkin(4) * t333 + pkin(12) * t361;
	t326 = pkin(4) * t361 + pkin(12) * t333;
	t331 = pkin(11) * t338 + qJ(3);
	t335 = sin(pkin(7));
	t339 = cos(pkin(7));
	t342 = sin(qJ(4));
	t346 = cos(qJ(4));
	t363 = t334 * t339;
	t371 = -(pkin(12) * t363 + t325 * t335) * t346 + (pkin(4) * t363 + t326 * t335) * t342 + pkin(10) + t331 * t339 - t322 * t335;
	t369 = t333 * t338;
	t343 = sin(qJ(2));
	t368 = t333 * t343;
	t367 = t333 * t346;
	t347 = cos(qJ(2));
	t366 = t333 * t347;
	t365 = t334 * t335;
	t336 = sin(pkin(6));
	t362 = t336 * t347;
	t360 = t338 * t339;
	t319 = t335 * t361 + t363;
	t340 = cos(pkin(6));
	t359 = t340 * t319;
	t358 = t340 * t343;
	t357 = t340 * t347;
	t316 = t340 * t335 * t333 + t336 * t343 * t337;
	t356 = t338 * t368;
	t320 = t337 * t360 - t365;
	t312 = t320 * t342 + t339 * t367;
	t351 = -t337 * t346 + t342 * t369;
	t299 = t312 * t357 - (t319 * t342 + t335 * t367) * t336 - t351 * t358;
	t318 = -t335 * t364 + t360;
	t317 = t335 * t338 + t337 * t363;
	t353 = -t317 * t347 + t334 * t368;
	t303 = t336 * t318 - t353 * t340;
	t341 = sin(qJ(5));
	t345 = cos(qJ(5));
	t355 = t299 * t341 + t303 * t345;
	t301 = (-pkin(12) * t365 + t325 * t339) * t346 + (pkin(4) * t365 - t326 * t339) * t342 + t322 * t339 + t335 * t331;
	t306 = (pkin(4) * t337 + pkin(12) * t369) * t346 - (pkin(4) * t369 - pkin(12) * t337) * t342 + pkin(2) + t334 * t333 * pkin(11) + t337 * pkin(3);
	t354 = t301 * t347 - t306 * t343;
	t352 = -t320 * t347 + t356;
	t349 = -(t336 * t356 - t359) * t342 + t312 * t362 + t316 * t346;
	t348 = cos(qJ(1));
	t344 = sin(qJ(1));
	t321 = -t337 * t347 + t339 * t368;
	t311 = t320 * t343 + t338 * t366;
	t310 = t317 * t343 + t334 * t366;
	t309 = t337 * t358 + (-t335 * t336 + t339 * t357) * t333;
	t305 = t312 * t343 + t351 * t347;
	t304 = -t336 * t319 - t352 * t340;
	t302 = t340 * t318 + t353 * t336;
	t300 = t305 * t341 + t310 * t345;
	t298 = t301 * t343 + t306 * t347 + pkin(1);
	t297 = t336 * t371 + t354 * t340;
	t1 = [(-t299 * t344 - t305 * t348) * t345 + (t303 * t344 + t310 * t348) * t341, t300 * t348 + t355 * t344, (t304 * t344 + t311 * t348) * t346 - t342 * (t309 * t344 + t321 * t348), t297 * t344 + t298 * t348 + 0; (t299 * t345 - t303 * t341) * t348 - t344 * (t305 * t345 - t310 * t341), t300 * t344 - t355 * t348, (-t304 * t346 + t309 * t342) * t348 + (t311 * t346 - t321 * t342) * t344, -t297 * t348 + t298 * t344 + 0; t341 * t302 + t349 * t345, t345 * t302 - t349 * t341, (t352 * t336 - t359) * t346 + (t333 * t339 * t362 + t316) * t342, -t354 * t336 + t371 * t340 + pkin(9) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:21:14
	% EndTime: 2020-11-04 22:21:15
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (455->115), mult. (1181->194), div. (0->0), fcn. (1449->18), ass. (0->88)
	t417 = cos(pkin(7));
	t416 = cos(pkin(8));
	t408 = t416 * pkin(11) + qJ(3);
	t412 = sin(pkin(8));
	t421 = sin(qJ(4));
	t426 = cos(qJ(4));
	t477 = t408 + (pkin(4) * t421 - pkin(12) * t426) * t412;
	t480 = t477 * t417 + pkin(10);
	t411 = sin(pkin(14));
	t415 = cos(pkin(14));
	t455 = t415 * t416;
	t403 = -t411 * pkin(4) + pkin(12) * t455;
	t420 = sin(qJ(5));
	t474 = pkin(13) * t420;
	t440 = t411 * t474 - t403;
	t478 = t415 * t412 * pkin(11) - (pkin(4) * t455 + t411 * pkin(12)) * t421;
	t442 = t411 * pkin(3) - t478;
	t425 = cos(qJ(5));
	t465 = t411 * t425;
	t479 = (t426 * (pkin(5) * t465 + t440) + t442) * t417;
	t438 = pkin(5) * t425 + t474;
	t473 = pkin(13) * t425;
	t472 = t415 * pkin(12);
	t413 = sin(pkin(7));
	t463 = t412 * t417;
	t397 = t413 * t416 + t415 * t463;
	t427 = cos(qJ(2));
	t437 = -pkin(5) * t420 + t473;
	t454 = t417 * t416;
	t461 = t413 * t412;
	t400 = t415 * t454 - t461;
	t450 = t421 * t400;
	t471 = (t437 * t397 - t413 * t477 + t438 * t450 + t479) * t427;
	t399 = t413 * t455 + t463;
	t470 = t399 * t421;
	t468 = t411 * t412;
	t467 = t411 * t416;
	t466 = t411 * t421;
	t464 = t411 * t426;
	t462 = t412 * t420;
	t414 = sin(pkin(6));
	t460 = t413 * t414;
	t459 = t413 * t426;
	t398 = -t415 * t461 + t454;
	t458 = t414 * t398;
	t422 = sin(qJ(2));
	t457 = t414 * t422;
	t456 = t414 * t427;
	t418 = cos(pkin(6));
	t453 = t418 * t422;
	t452 = t418 * t427;
	t451 = t420 * t397;
	t449 = t421 * t425;
	t448 = t425 * t426;
	t445 = t411 * t457;
	t387 = t418 * t398 + t412 * t445;
	t446 = pkin(12) * t467;
	t444 = t411 * t453;
	t443 = pkin(4) + t474;
	t441 = t411 * pkin(5) * t462 + t415 * pkin(3) + pkin(11) * t468 + pkin(2);
	t390 = t417 * t464 + t450;
	t384 = t390 * t425 - t451;
	t396 = -t411 * t460 + t415 * t453;
	t432 = t414 * t399 + t416 * t444;
	t376 = t384 * t452 + t396 * t448 - t432 * t449 + t420 * (t412 * t444 - t458);
	t391 = t400 * t426 - t417 * t466;
	t379 = t391 * t452 - t421 * t396 - t426 * t432;
	t419 = sin(qJ(6));
	t424 = cos(qJ(6));
	t436 = t376 * t419 + t379 * t424;
	t401 = -t415 * t426 + t416 * t466;
	t435 = (-(t411 * t459 + t470) * t414 + (t390 * t427 - t401 * t422) * t418) * t420 + (t458 + (t397 * t427 - t422 * t468) * t418) * t425;
	t381 = t384 * t422 - (t415 * t448 - t411 * (t416 * t449 - t462)) * t427;
	t388 = t418 * t399 - t416 * t445;
	t395 = t418 * t413 * t411 + t415 * t457;
	t434 = t388 * t421 + t426 * t395;
	t433 = pkin(4) + t438;
	t429 = pkin(13) * t412 * t465 + (t433 * t467 - t472) * t421 - (t415 * t433 + t446) * t426 - t441;
	t428 = cos(qJ(1));
	t423 = sin(qJ(1));
	t385 = (t415 * t421 + t416 * t464) * t427 + t391 * t422;
	t382 = (t390 * t422 + t401 * t427) * t420 + t425 * (t422 * t397 + t427 * t468);
	t378 = t388 * t426 + t391 * t456 - t421 * t395;
	t377 = t381 * t419 + t385 * t424;
	t375 = t384 * t456 + t420 * t387 + t425 * t434;
	t373 = -t429 * t427 + pkin(1) + ((pkin(4) * t461 - t438 * t400) * t421 - t412 * pkin(12) * t459 - t397 * t473 + pkin(5) * t451 + t413 * t408 - t479) * t422;
	t372 = (t398 * t437 - t438 * t470 - t480) * t414 + (-t422 * t429 + t471) * t418 + (t403 * t426 + (-t426 * t438 - pkin(3)) * t411 + t478) * t460;
	t1 = [(-t376 * t423 - t381 * t428) * t424 + t419 * (t379 * t423 + t428 * t385), t428 * t377 + t423 * t436, -t382 * t428 - t423 * t435, -t372 * t423 + t373 * t428 + 0; (t376 * t424 - t419 * t379) * t428 + t423 * (-t381 * t424 + t419 * t385), t423 * t377 - t428 * t436, -t382 * t423 + t428 * t435, t372 * t428 + t373 * t423 + 0; t375 * t424 - t378 * t419, -t375 * t419 - t378 * t424, (t390 * t456 + t434) * t420 - t425 * (-t397 * t456 + t387), pkin(9) + 0 + (pkin(5) * t434 - t387 * pkin(13)) * t425 + ((pkin(5) * t398 + pkin(13) * t470) * t420 + (t426 * t440 + t442) * t413 + t480) * t418 + (t471 + (-(t443 * t467 - t472) * t421 + (t415 * t443 + t446) * t426 + t441) * t422) * t414; 0, 0, 0, 1;];
	Tc_mdh = t1;
end