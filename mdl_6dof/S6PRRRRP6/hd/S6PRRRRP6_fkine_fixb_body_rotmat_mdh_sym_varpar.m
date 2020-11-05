% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:20
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:54
	% EndTime: 2020-11-04 21:19:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:54
	% EndTime: 2020-11-04 21:19:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t235 = cos(pkin(12));
	t234 = sin(pkin(12));
	t1 = [t235, -t234, 0, 0; t234, t235, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:54
	% EndTime: 2020-11-04 21:19:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t236 = sin(pkin(12));
	t237 = sin(pkin(6));
	t245 = t236 * t237;
	t238 = cos(pkin(12));
	t244 = t238 * t237;
	t239 = cos(pkin(6));
	t240 = sin(qJ(2));
	t243 = t239 * t240;
	t241 = cos(qJ(2));
	t242 = t239 * t241;
	t1 = [-t236 * t243 + t238 * t241, -t236 * t242 - t238 * t240, t245, t238 * pkin(1) + pkin(8) * t245 + 0; t236 * t241 + t238 * t243, -t236 * t240 + t238 * t242, -t244, t236 * pkin(1) - pkin(8) * t244 + 0; t237 * t240, t237 * t241, t239, t239 * pkin(8) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:54
	% EndTime: 2020-11-04 21:19:54
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->35), mult. (112->66), div. (0->0), fcn. (146->10), ass. (0->30)
	t252 = sin(pkin(7));
	t274 = pkin(9) * t252;
	t251 = sin(pkin(12));
	t273 = t251 * pkin(2);
	t254 = cos(pkin(12));
	t272 = t254 * pkin(2);
	t256 = cos(pkin(6));
	t271 = t252 * t256;
	t260 = cos(qJ(2));
	t270 = t252 * t260;
	t255 = cos(pkin(7));
	t250 = t255 * pkin(9) + pkin(8);
	t253 = sin(pkin(6));
	t269 = t253 * t250;
	t268 = t253 * t255;
	t258 = sin(qJ(2));
	t267 = t255 * t258;
	t266 = t255 * t260;
	t265 = t256 * t258;
	t264 = t256 * t260;
	t263 = t251 * t274;
	t262 = t254 * t274;
	t261 = -t252 * t253 + t255 * t264;
	t259 = cos(qJ(3));
	t257 = sin(qJ(3));
	t249 = t251 * t260 + t254 * t265;
	t248 = t251 * t265 - t254 * t260;
	t247 = -t251 * t267 + t261 * t254;
	t246 = -t261 * t251 - t254 * t267;
	t1 = [t246 * t257 - t259 * t248, t246 * t259 + t257 * t248, (t251 * t264 + t254 * t258) * t252 + t251 * t268, (t256 * t263 + t272) * t260 + (-t256 * t273 + t262) * t258 + t251 * t269 + t254 * pkin(1) + 0; t247 * t257 + t249 * t259, t247 * t259 - t249 * t257, -(-t251 * t258 + t254 * t264) * t252 - t254 * t268, (-t256 * t262 + t273) * t260 + (t256 * t272 + t263) * t258 - t254 * t269 + t251 * pkin(1) + 0; t257 * t271 + (t257 * t266 + t258 * t259) * t253, t259 * t271 + (-t257 * t258 + t259 * t266) * t253, -t253 * t270 + t256 * t255, t250 * t256 + qJ(1) + 0 + (pkin(2) * t258 - pkin(9) * t270) * t253; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:54
	% EndTime: 2020-11-04 21:19:54
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (91->59), mult. (247->109), div. (0->0), fcn. (302->12), ass. (0->42)
	t286 = sin(pkin(6));
	t285 = sin(pkin(7));
	t288 = cos(pkin(7));
	t291 = sin(qJ(3));
	t294 = cos(qJ(3));
	t303 = pkin(3) * t291 - pkin(10) * t294;
	t300 = t288 * pkin(9) + t303 * t285 + pkin(8);
	t321 = t300 * t286;
	t292 = sin(qJ(2));
	t295 = cos(qJ(2));
	t319 = t285 * pkin(9);
	t299 = -t303 * t288 + t319;
	t302 = pkin(3) * t294 + pkin(10) * t291 + pkin(2);
	t320 = t302 * t292 - t299 * t295;
	t284 = sin(pkin(12));
	t318 = t284 * t288;
	t317 = t285 * t286;
	t316 = t286 * t288;
	t287 = cos(pkin(12));
	t289 = cos(pkin(6));
	t315 = t287 * t289;
	t314 = t288 * t291;
	t313 = t288 * t292;
	t312 = t288 * t295;
	t311 = t289 * t288;
	t310 = t289 * t291;
	t309 = t289 * t292;
	t308 = t289 * t294;
	t307 = t289 * t295;
	t306 = t287 * t311;
	t305 = t288 * t310;
	t304 = t291 * t317;
	t301 = t288 * t307 - t317;
	t298 = -(t284 * t305 - t287 * t294) * t295 - (t284 * t308 + t287 * t314) * t292 + t284 * t304;
	t297 = -(t284 * t294 + t287 * t305) * t295 - (-t284 * t314 + t287 * t308) * t292 + t287 * t304;
	t293 = cos(qJ(4));
	t290 = sin(qJ(4));
	t282 = t295 * t317 - t311;
	t277 = -t285 * t310 + (-t291 * t312 - t294 * t292) * t286;
	t276 = t287 * t316 + (-t284 * t292 + t287 * t307) * t285;
	t275 = t285 * t292 * t287 + (t285 * t307 + t316) * t284;
	t1 = [t290 * t275 + t298 * t293, t293 * t275 - t298 * t290, (t301 * t284 + t287 * t313) * t294 - t291 * (t284 * t309 - t287 * t295), 0 + (t299 * t292 + t302 * t295 + pkin(1)) * t287 + (-t320 * t289 + t321) * t284; -t290 * t276 - t297 * t293, -t293 * t276 + t297 * t290, (t284 * t313 - t301 * t287) * t294 + (t284 * t295 + t287 * t309) * t291, ((t284 * pkin(3) - pkin(10) * t306) * t294 + (pkin(3) * t306 + t284 * pkin(10)) * t291 - t315 * t319 + t284 * pkin(2)) * t295 + ((pkin(3) * t315 + pkin(10) * t318) * t294 + (-pkin(3) * t318 + pkin(10) * t315) * t291 + pkin(2) * t315 + t284 * t319) * t292 + t284 * pkin(1) + 0 - t287 * t321; -t277 * t293 - t290 * t282, t277 * t290 - t293 * t282, -t285 * t308 + (t291 * t292 - t294 * t312) * t286, t320 * t286 + t300 * t289 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:54
	% EndTime: 2020-11-04 21:19:55
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (167->81), mult. (469->151), div. (0->0), fcn. (558->14), ass. (0->62)
	t344 = cos(pkin(7));
	t351 = cos(qJ(4));
	t381 = t351 * t344;
	t347 = sin(qJ(4));
	t390 = t344 * t347;
	t411 = -pkin(4) * t390 - t344 * pkin(9) + pkin(11) * t381 - pkin(8);
	t410 = t351 * pkin(4) + t347 * pkin(11) + pkin(3);
	t341 = sin(pkin(7));
	t348 = sin(qJ(3));
	t352 = cos(qJ(3));
	t402 = t352 * pkin(10);
	t409 = (t348 * t410 - t402) * t344 - (t347 * pkin(4) - t351 * pkin(11) + pkin(9)) * t341;
	t340 = sin(pkin(12));
	t349 = sin(qJ(2));
	t343 = cos(pkin(12));
	t392 = t343 * t341;
	t342 = sin(pkin(6));
	t395 = t342 * t344;
	t328 = t340 * t395 + t349 * t392;
	t388 = t344 * t349;
	t399 = t341 * t342;
	t331 = -t340 * t399 + t343 * t388;
	t332 = -t341 * t347 + t348 * t381;
	t353 = cos(qJ(2));
	t345 = cos(pkin(6));
	t382 = t349 * t352;
	t374 = t345 * t382;
	t380 = t351 * t352;
	t401 = t340 * t345;
	t406 = (t331 * t348 + t340 * t374) * t351 - (-t332 * t401 + t343 * t380) * t353 - t347 * t328;
	t329 = t340 * t388 + t342 * t392;
	t335 = t343 * t395;
	t391 = t343 * t345;
	t400 = t340 * t349;
	t405 = (-t329 * t348 + t343 * t374) * t351 + (t332 * t391 + t340 * t380) * t353 - t347 * (-t341 * t400 + t335);
	t357 = pkin(10) * t348 + t352 * t410 + pkin(2);
	t404 = t411 * t342 + t399 * t402;
	t398 = t341 * t345;
	t397 = t341 * t348;
	t394 = t342 * t349;
	t393 = t342 * t353;
	t389 = t344 * t348;
	t387 = t344 * t353;
	t386 = t345 * t348;
	t385 = t345 * t352;
	t384 = t345 * t353;
	t383 = t348 * t349;
	t377 = t342 * t397;
	t376 = t344 * t386;
	t375 = t344 * t385;
	t373 = t345 * t383;
	t366 = t340 * t377;
	t365 = t343 * t377;
	t356 = t357 * t343;
	t354 = t409 * t340;
	t350 = cos(qJ(5));
	t346 = sin(qJ(5));
	t327 = t341 * t385 + (t352 * t387 - t383) * t342;
	t324 = -t345 * (t351 * t397 + t390) + (-t332 * t353 - t349 * t380) * t342;
	t323 = (-t340 * t348 + t343 * t375) * t353 - t329 * t352 - t343 * t373;
	t322 = (t340 * t375 + t343 * t348) * t353 + t331 * t352 - t340 * t373;
	t1 = [t346 * t322 - t406 * t350, t350 * t322 + t406 * t346, ((-t340 * t376 + t343 * t352) * t353 + (-t340 * t385 - t343 * t389) * t349 + t366) * t347 - t351 * (t341 * t340 * t384 + t328), (-t345 * t354 + t356) * t353 + (-t343 * t409 - t357 * t401) * t349 + t410 * t366 + t343 * pkin(1) + 0 - t404 * t340; -t346 * t323 + t405 * t350, -t323 * t350 - t405 * t346, ((t340 * t352 + t343 * t376) * t353 + (-t340 * t389 + t343 * t385) * t349 - t365) * t347 + t351 * (t335 + (t343 * t384 - t400) * t341), (t357 * t340 + t409 * t391) * t353 + (t345 * t356 - t354) * t349 - t410 * t365 + t340 * pkin(1) + 0 + t404 * t343; -t324 * t350 - t346 * t327, t324 * t346 - t350 * t327, (t341 * t386 + (t348 * t387 + t382) * t342) * t347 + t351 * (t341 * t393 - t345 * t344), t409 * t393 + (pkin(10) * t394 + t398 * t410) * t348 + (-pkin(10) * t398 + t394 * t410) * t352 + pkin(2) * t394 + qJ(1) + 0 - t411 * t345; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:55
	% EndTime: 2020-11-04 21:19:55
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (235->87), mult. (569->152), div. (0->0), fcn. (658->14), ass. (0->64)
	t436 = sin(qJ(5));
	t440 = cos(qJ(5));
	t422 = pkin(5) * t440 + qJ(6) * t436 + pkin(4);
	t434 = cos(pkin(7));
	t437 = sin(qJ(4));
	t478 = t434 * t437;
	t441 = cos(qJ(4));
	t488 = pkin(11) * t441;
	t495 = -t422 * t478 + t434 * t488;
	t431 = sin(pkin(7));
	t438 = sin(qJ(3));
	t442 = cos(qJ(3));
	t450 = pkin(5) * t436 - qJ(6) * t440 + pkin(10);
	t486 = t450 * t442;
	t493 = t437 * pkin(11) + t422 * t441 + pkin(3);
	t492 = (-t422 * t437 - pkin(9) + t488) * t431 + (t438 * t493 - t486) * t434;
	t435 = cos(pkin(6));
	t439 = sin(qJ(2));
	t473 = t435 * t439;
	t443 = cos(qJ(2));
	t471 = t435 * t443;
	t432 = sin(pkin(6));
	t485 = t431 * t432;
	t490 = t495 * t432 + t485 * t486;
	t447 = t450 * t438 + t442 * t493 + pkin(2);
	t484 = t431 * t435;
	t483 = t431 * t439;
	t433 = cos(pkin(12));
	t482 = t432 * t433;
	t481 = t432 * t434;
	t480 = t432 * t439;
	t479 = t432 * t443;
	t477 = t434 * t438;
	t476 = t434 * t439;
	t475 = t434 * t443;
	t474 = t435 * t438;
	t472 = t435 * t442;
	t470 = t438 * t439;
	t469 = t438 * t441;
	t468 = t439 * t442;
	t430 = sin(pkin(12));
	t467 = t430 * t481 + t433 * t483;
	t465 = t430 * t471;
	t464 = t438 * t485;
	t463 = t433 * t471;
	t461 = t434 * t474;
	t460 = t434 * t472;
	t459 = t435 * t470;
	t455 = t430 * t464;
	t453 = t433 * t464;
	t446 = t447 * t433;
	t444 = t492 * t430;
	t427 = pkin(9) * t434 + pkin(8);
	t425 = t433 * t481;
	t421 = -t431 * t437 + t434 * t469;
	t419 = -t430 * t485 + t433 * t476;
	t418 = t430 * t476 + t431 * t482;
	t417 = t431 * t472 + (t442 * t475 - t470) * t432;
	t416 = t435 * (t431 * t469 + t478) + (t421 * t443 + t441 * t468) * t432;
	t415 = (-t430 * t438 + t433 * t460) * t443 - t418 * t442 - t433 * t459;
	t414 = (t430 * t460 + t433 * t438) * t443 + t419 * t442 - t430 * t459;
	t413 = -t421 * t465 + t437 * t467 + (-t438 * t419 + (-t430 * t473 + t433 * t443) * t442) * t441;
	t412 = t421 * t463 - t437 * (-t430 * t483 + t425) + (-t438 * t418 + (t430 * t443 + t433 * t473) * t442) * t441;
	t1 = [t413 * t440 + t414 * t436, ((-t430 * t461 + t433 * t442) * t443 + (-t430 * t472 - t433 * t477) * t439 + t455) * t437 - t441 * (t431 * t465 + t467), t413 * t436 - t414 * t440, (-t435 * t444 + t446) * t443 + t493 * t455 + 0 + (-t439 * t492 + pkin(1)) * t433 + (t432 * t427 - t447 * t473 - t490) * t430; t412 * t440 - t415 * t436, ((t430 * t442 + t433 * t461) * t443 + (-t430 * t477 + t433 * t472) * t439 - t453) * t437 + t441 * (t425 + (-t430 * t439 + t463) * t431), t412 * t436 + t415 * t440, (t435 * t446 - t444) * t439 - t493 * t453 - t427 * t482 + 0 + (t447 * t443 + pkin(1)) * t430 + (t492 * t471 + t490) * t433; t416 * t440 - t417 * t436, (t431 * t474 + (t438 * t475 + t468) * t432) * t437 + t441 * (t431 * t479 - t435 * t434), t416 * t436 + t417 * t440, t492 * t479 + (t450 * t480 + t484 * t493) * t438 + (-t450 * t484 + t480 * t493) * t442 + pkin(2) * t480 + qJ(1) + 0 + (t427 - t495) * t435; 0, 0, 0, 1;];
	Tc_mdh = t1;
end