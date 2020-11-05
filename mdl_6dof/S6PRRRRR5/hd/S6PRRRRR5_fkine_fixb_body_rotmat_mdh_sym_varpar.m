% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:22
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:59
	% EndTime: 2020-11-04 21:21:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:59
	% EndTime: 2020-11-04 21:21:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t211 = cos(pkin(13));
	t210 = sin(pkin(13));
	t1 = [t211, -t210, 0, 0; t210, t211, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:59
	% EndTime: 2020-11-04 21:21:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t212 = sin(pkin(13));
	t213 = sin(pkin(6));
	t221 = t212 * t213;
	t214 = cos(pkin(13));
	t220 = t214 * t213;
	t215 = cos(pkin(6));
	t216 = sin(qJ(2));
	t219 = t215 * t216;
	t217 = cos(qJ(2));
	t218 = t215 * t217;
	t1 = [-t212 * t219 + t214 * t217, -t212 * t218 - t214 * t216, t221, t214 * pkin(1) + pkin(8) * t221 + 0; t212 * t217 + t214 * t219, -t212 * t216 + t214 * t218, -t220, t212 * pkin(1) - pkin(8) * t220 + 0; t213 * t216, t213 * t217, t215, t215 * pkin(8) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:59
	% EndTime: 2020-11-04 21:21:59
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->35), mult. (112->66), div. (0->0), fcn. (146->10), ass. (0->30)
	t228 = sin(pkin(7));
	t250 = pkin(9) * t228;
	t227 = sin(pkin(13));
	t249 = t227 * pkin(2);
	t230 = cos(pkin(13));
	t248 = t230 * pkin(2);
	t232 = cos(pkin(6));
	t247 = t228 * t232;
	t236 = cos(qJ(2));
	t246 = t228 * t236;
	t231 = cos(pkin(7));
	t226 = t231 * pkin(9) + pkin(8);
	t229 = sin(pkin(6));
	t245 = t229 * t226;
	t244 = t229 * t231;
	t234 = sin(qJ(2));
	t243 = t231 * t234;
	t242 = t231 * t236;
	t241 = t232 * t234;
	t240 = t232 * t236;
	t239 = t227 * t250;
	t238 = t230 * t250;
	t237 = -t228 * t229 + t231 * t240;
	t235 = cos(qJ(3));
	t233 = sin(qJ(3));
	t225 = t227 * t236 + t230 * t241;
	t224 = t227 * t241 - t230 * t236;
	t223 = -t227 * t243 + t237 * t230;
	t222 = -t237 * t227 - t230 * t243;
	t1 = [t222 * t233 - t235 * t224, t222 * t235 + t233 * t224, (t227 * t240 + t230 * t234) * t228 + t227 * t244, (t232 * t239 + t248) * t236 + (-t232 * t249 + t238) * t234 + t227 * t245 + t230 * pkin(1) + 0; t223 * t233 + t225 * t235, t223 * t235 - t225 * t233, -(-t227 * t234 + t230 * t240) * t228 - t230 * t244, (-t232 * t238 + t249) * t236 + (t232 * t248 + t239) * t234 - t230 * t245 + t227 * pkin(1) + 0; t233 * t247 + (t233 * t242 + t234 * t235) * t229, t235 * t247 + (-t233 * t234 + t235 * t242) * t229, -t229 * t246 + t232 * t231, t226 * t232 + qJ(1) + 0 + (pkin(2) * t234 - pkin(9) * t246) * t229; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:59
	% EndTime: 2020-11-04 21:21:59
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (91->59), mult. (247->109), div. (0->0), fcn. (302->12), ass. (0->42)
	t262 = sin(pkin(6));
	t261 = sin(pkin(7));
	t264 = cos(pkin(7));
	t267 = sin(qJ(3));
	t270 = cos(qJ(3));
	t279 = pkin(3) * t267 - pkin(10) * t270;
	t276 = t264 * pkin(9) + t279 * t261 + pkin(8);
	t297 = t276 * t262;
	t268 = sin(qJ(2));
	t271 = cos(qJ(2));
	t295 = t261 * pkin(9);
	t275 = -t279 * t264 + t295;
	t278 = pkin(3) * t270 + pkin(10) * t267 + pkin(2);
	t296 = t278 * t268 - t275 * t271;
	t260 = sin(pkin(13));
	t294 = t260 * t264;
	t293 = t261 * t262;
	t292 = t262 * t264;
	t263 = cos(pkin(13));
	t265 = cos(pkin(6));
	t291 = t263 * t265;
	t290 = t264 * t267;
	t289 = t264 * t268;
	t288 = t264 * t271;
	t287 = t265 * t264;
	t286 = t265 * t267;
	t285 = t265 * t268;
	t284 = t265 * t270;
	t283 = t265 * t271;
	t282 = t263 * t287;
	t281 = t264 * t286;
	t280 = t267 * t293;
	t277 = t264 * t283 - t293;
	t274 = -(t260 * t281 - t263 * t270) * t271 - (t260 * t284 + t263 * t290) * t268 + t260 * t280;
	t273 = -(t260 * t270 + t263 * t281) * t271 - (-t260 * t290 + t263 * t284) * t268 + t263 * t280;
	t269 = cos(qJ(4));
	t266 = sin(qJ(4));
	t258 = t271 * t293 - t287;
	t253 = -t261 * t286 + (-t267 * t288 - t270 * t268) * t262;
	t252 = t263 * t292 + (-t260 * t268 + t263 * t283) * t261;
	t251 = t261 * t268 * t263 + (t261 * t283 + t292) * t260;
	t1 = [t266 * t251 + t274 * t269, t269 * t251 - t274 * t266, (t277 * t260 + t263 * t289) * t270 - t267 * (t260 * t285 - t263 * t271), 0 + (t275 * t268 + t278 * t271 + pkin(1)) * t263 + (-t296 * t265 + t297) * t260; -t266 * t252 - t273 * t269, -t269 * t252 + t273 * t266, (t260 * t289 - t277 * t263) * t270 + (t260 * t271 + t263 * t285) * t267, ((t260 * pkin(3) - pkin(10) * t282) * t270 + (pkin(3) * t282 + t260 * pkin(10)) * t267 - t291 * t295 + t260 * pkin(2)) * t271 + ((pkin(3) * t291 + pkin(10) * t294) * t270 + (-pkin(3) * t294 + pkin(10) * t291) * t267 + pkin(2) * t291 + t260 * t295) * t268 + t260 * pkin(1) + 0 - t263 * t297; -t253 * t269 - t266 * t258, t253 * t266 - t269 * t258, -t261 * t284 + (t267 * t268 - t270 * t288) * t262, t296 * t262 + t276 * t265 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:59
	% EndTime: 2020-11-04 21:21:59
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (167->81), mult. (469->151), div. (0->0), fcn. (558->14), ass. (0->62)
	t320 = cos(pkin(7));
	t327 = cos(qJ(4));
	t357 = t327 * t320;
	t323 = sin(qJ(4));
	t366 = t320 * t323;
	t387 = -pkin(4) * t366 - t320 * pkin(9) + pkin(11) * t357 - pkin(8);
	t386 = t327 * pkin(4) + t323 * pkin(11) + pkin(3);
	t317 = sin(pkin(7));
	t324 = sin(qJ(3));
	t328 = cos(qJ(3));
	t378 = t328 * pkin(10);
	t385 = (t324 * t386 - t378) * t320 - (t323 * pkin(4) - t327 * pkin(11) + pkin(9)) * t317;
	t316 = sin(pkin(13));
	t325 = sin(qJ(2));
	t319 = cos(pkin(13));
	t368 = t319 * t317;
	t318 = sin(pkin(6));
	t371 = t318 * t320;
	t304 = t316 * t371 + t325 * t368;
	t364 = t320 * t325;
	t375 = t317 * t318;
	t307 = -t316 * t375 + t319 * t364;
	t308 = -t317 * t323 + t324 * t357;
	t329 = cos(qJ(2));
	t321 = cos(pkin(6));
	t358 = t325 * t328;
	t350 = t321 * t358;
	t356 = t327 * t328;
	t377 = t316 * t321;
	t382 = (t307 * t324 + t316 * t350) * t327 - (-t308 * t377 + t319 * t356) * t329 - t323 * t304;
	t305 = t316 * t364 + t318 * t368;
	t311 = t319 * t371;
	t367 = t319 * t321;
	t376 = t316 * t325;
	t381 = (-t305 * t324 + t319 * t350) * t327 + (t308 * t367 + t316 * t356) * t329 - t323 * (-t317 * t376 + t311);
	t333 = pkin(10) * t324 + t328 * t386 + pkin(2);
	t380 = t387 * t318 + t375 * t378;
	t374 = t317 * t321;
	t373 = t317 * t324;
	t370 = t318 * t325;
	t369 = t318 * t329;
	t365 = t320 * t324;
	t363 = t320 * t329;
	t362 = t321 * t324;
	t361 = t321 * t328;
	t360 = t321 * t329;
	t359 = t324 * t325;
	t353 = t318 * t373;
	t352 = t320 * t362;
	t351 = t320 * t361;
	t349 = t321 * t359;
	t342 = t316 * t353;
	t341 = t319 * t353;
	t332 = t333 * t319;
	t330 = t385 * t316;
	t326 = cos(qJ(5));
	t322 = sin(qJ(5));
	t303 = t317 * t361 + (t328 * t363 - t359) * t318;
	t300 = -t321 * (t327 * t373 + t366) + (-t308 * t329 - t325 * t356) * t318;
	t299 = (-t316 * t324 + t319 * t351) * t329 - t305 * t328 - t319 * t349;
	t298 = (t316 * t351 + t319 * t324) * t329 + t307 * t328 - t316 * t349;
	t1 = [t322 * t298 - t382 * t326, t326 * t298 + t382 * t322, ((-t316 * t352 + t319 * t328) * t329 + (-t316 * t361 - t319 * t365) * t325 + t342) * t323 - t327 * (t317 * t316 * t360 + t304), (-t321 * t330 + t332) * t329 + (-t319 * t385 - t333 * t377) * t325 + t386 * t342 + t319 * pkin(1) + 0 - t380 * t316; -t322 * t299 + t381 * t326, -t326 * t299 - t381 * t322, ((t316 * t328 + t319 * t352) * t329 + (-t316 * t365 + t319 * t361) * t325 - t341) * t323 + t327 * (t311 + (t319 * t360 - t376) * t317), (t333 * t316 + t385 * t367) * t329 + (t321 * t332 - t330) * t325 - t386 * t341 + t316 * pkin(1) + 0 + t380 * t319; -t300 * t326 - t322 * t303, t300 * t322 - t326 * t303, (t317 * t362 + (t324 * t363 + t358) * t318) * t323 + t327 * (t317 * t369 - t321 * t320), t385 * t369 + (pkin(10) * t370 + t374 * t386) * t324 + (-pkin(10) * t374 + t370 * t386) * t328 + pkin(2) * t370 + qJ(1) + 0 - t387 * t321; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:21:59
	% EndTime: 2020-11-04 21:22:00
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (198->80), mult. (504->132), div. (0->0), fcn. (639->16), ass. (0->57)
	t412 = sin(pkin(6));
	t411 = sin(pkin(7));
	t414 = cos(pkin(7));
	t418 = sin(qJ(3));
	t421 = cos(qJ(3));
	t431 = pkin(3) * t418 - pkin(10) * t421;
	t428 = t414 * pkin(9) + t431 * t411 + pkin(8);
	t450 = t428 * t412;
	t419 = sin(qJ(2));
	t422 = cos(qJ(2));
	t447 = t411 * pkin(9);
	t427 = -t431 * t414 + t447;
	t430 = pkin(3) * t421 + pkin(10) * t418 + pkin(2);
	t449 = t430 * t419 - t427 * t422;
	t448 = pkin(5) * sin(qJ(5));
	t410 = sin(pkin(13));
	t446 = t410 * t414;
	t445 = t411 * t412;
	t444 = t412 * t414;
	t413 = cos(pkin(13));
	t415 = cos(pkin(6));
	t443 = t413 * t415;
	t442 = t414 * t418;
	t441 = t414 * t419;
	t440 = t414 * t422;
	t439 = t415 * t414;
	t438 = t415 * t418;
	t437 = t415 * t419;
	t436 = t415 * t421;
	t435 = t415 * t422;
	t434 = t413 * t439;
	t433 = t414 * t438;
	t432 = t418 * t445;
	t429 = t414 * t435 - t445;
	t426 = -(t410 * t433 - t413 * t421) * t422 - (t410 * t436 + t413 * t442) * t419 + t410 * t432;
	t425 = -(t410 * t421 + t413 * t433) * t422 - (-t410 * t442 + t413 * t436) * t419 + t413 * t432;
	t423 = -pkin(12) - pkin(11);
	t420 = cos(qJ(4));
	t417 = sin(qJ(4));
	t409 = qJ(5) + qJ(6);
	t408 = cos(t409);
	t407 = sin(t409);
	t406 = cos(qJ(5)) * pkin(5) + pkin(4);
	t404 = t422 * t445 - t439;
	t399 = -t411 * t436 + (t418 * t419 - t421 * t440) * t412;
	t398 = -t411 * t438 + (-t418 * t440 - t421 * t419) * t412;
	t397 = t413 * t444 + (-t410 * t419 + t413 * t435) * t411;
	t396 = t411 * t419 * t413 + (t411 * t435 + t444) * t410;
	t395 = -t398 * t420 - t417 * t404;
	t394 = t398 * t417 - t420 * t404;
	t393 = (t410 * t441 - t429 * t413) * t421 + (t410 * t422 + t413 * t437) * t418;
	t392 = (t429 * t410 + t413 * t441) * t421 - t418 * (t410 * t437 - t413 * t422);
	t391 = -t420 * t397 + t425 * t417;
	t390 = -t417 * t397 - t425 * t420;
	t389 = t417 * t396 + t426 * t420;
	t388 = t420 * t396 - t426 * t417;
	t1 = [t389 * t408 + t392 * t407, -t389 * t407 + t392 * t408, -t388, t392 * t448 + t388 * t423 + t389 * t406 + 0 + (t427 * t419 + t430 * t422 + pkin(1)) * t413 + (-t449 * t415 + t450) * t410; t390 * t408 + t393 * t407, -t390 * t407 + t393 * t408, -t391, t390 * t406 + t391 * t423 + t393 * t448 + ((t410 * pkin(3) - pkin(10) * t434) * t421 + (pkin(3) * t434 + t410 * pkin(10)) * t418 - t443 * t447 + t410 * pkin(2)) * t422 + ((pkin(3) * t443 + pkin(10) * t446) * t421 + (-pkin(3) * t446 + pkin(10) * t443) * t418 + pkin(2) * t443 + t410 * t447) * t419 + t410 * pkin(1) + 0 - t413 * t450; t395 * t408 + t399 * t407, -t395 * t407 + t399 * t408, -t394, t394 * t423 + t395 * t406 + t399 * t448 + t449 * t412 + t428 * t415 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end