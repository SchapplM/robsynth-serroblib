% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:17
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:28
	% EndTime: 2020-11-04 21:17:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:28
	% EndTime: 2020-11-04 21:17:28
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t214 = cos(pkin(12));
	t213 = sin(pkin(12));
	t1 = [t214, -t213, 0, 0; t213, t214, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:28
	% EndTime: 2020-11-04 21:17:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t215 = sin(pkin(12));
	t216 = sin(pkin(6));
	t224 = t215 * t216;
	t217 = cos(pkin(12));
	t223 = t217 * t216;
	t218 = cos(pkin(6));
	t219 = sin(qJ(2));
	t222 = t218 * t219;
	t220 = cos(qJ(2));
	t221 = t218 * t220;
	t1 = [-t215 * t222 + t217 * t220, -t215 * t221 - t217 * t219, t224, t217 * pkin(1) + pkin(8) * t224 + 0; t215 * t220 + t217 * t222, -t215 * t219 + t217 * t221, -t223, t215 * pkin(1) - pkin(8) * t223 + 0; t216 * t219, t216 * t220, t218, t218 * pkin(8) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:28
	% EndTime: 2020-11-04 21:17:28
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->35), mult. (112->66), div. (0->0), fcn. (146->10), ass. (0->30)
	t231 = sin(pkin(7));
	t253 = pkin(9) * t231;
	t230 = sin(pkin(12));
	t252 = t230 * pkin(2);
	t233 = cos(pkin(12));
	t251 = t233 * pkin(2);
	t235 = cos(pkin(6));
	t250 = t231 * t235;
	t239 = cos(qJ(2));
	t249 = t231 * t239;
	t234 = cos(pkin(7));
	t229 = t234 * pkin(9) + pkin(8);
	t232 = sin(pkin(6));
	t248 = t232 * t229;
	t247 = t232 * t234;
	t237 = sin(qJ(2));
	t246 = t234 * t237;
	t245 = t234 * t239;
	t244 = t235 * t237;
	t243 = t235 * t239;
	t242 = t230 * t253;
	t241 = t233 * t253;
	t240 = -t231 * t232 + t234 * t243;
	t238 = cos(qJ(3));
	t236 = sin(qJ(3));
	t228 = t230 * t244 - t233 * t239;
	t227 = t230 * t239 + t233 * t244;
	t226 = -t230 * t246 + t240 * t233;
	t225 = -t240 * t230 - t233 * t246;
	t1 = [t225 * t236 - t238 * t228, t225 * t238 + t236 * t228, (t230 * t243 + t233 * t237) * t231 + t230 * t247, (t235 * t242 + t251) * t239 + (-t235 * t252 + t241) * t237 + t230 * t248 + t233 * pkin(1) + 0; t226 * t236 + t227 * t238, t226 * t238 - t227 * t236, -(-t230 * t237 + t233 * t243) * t231 - t233 * t247, (-t235 * t241 + t252) * t239 + (t235 * t251 + t242) * t237 - t233 * t248 + t230 * pkin(1) + 0; t236 * t250 + (t236 * t245 + t237 * t238) * t232, t238 * t250 + (-t236 * t237 + t238 * t245) * t232, -t232 * t249 + t235 * t234, t229 * t235 + qJ(1) + 0 + (pkin(2) * t237 - pkin(9) * t249) * t232; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:28
	% EndTime: 2020-11-04 21:17:28
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (91->59), mult. (247->109), div. (0->0), fcn. (302->12), ass. (0->42)
	t265 = sin(pkin(6));
	t264 = sin(pkin(7));
	t267 = cos(pkin(7));
	t270 = sin(qJ(3));
	t273 = cos(qJ(3));
	t282 = pkin(3) * t270 - pkin(10) * t273;
	t279 = t267 * pkin(9) + t282 * t264 + pkin(8);
	t300 = t279 * t265;
	t271 = sin(qJ(2));
	t274 = cos(qJ(2));
	t298 = t264 * pkin(9);
	t278 = -t282 * t267 + t298;
	t281 = pkin(3) * t273 + pkin(10) * t270 + pkin(2);
	t299 = t281 * t271 - t278 * t274;
	t263 = sin(pkin(12));
	t297 = t263 * t267;
	t296 = t264 * t265;
	t295 = t265 * t267;
	t266 = cos(pkin(12));
	t268 = cos(pkin(6));
	t294 = t266 * t268;
	t293 = t267 * t270;
	t292 = t267 * t271;
	t291 = t267 * t274;
	t290 = t268 * t267;
	t289 = t268 * t270;
	t288 = t268 * t271;
	t287 = t268 * t273;
	t286 = t268 * t274;
	t285 = t266 * t290;
	t284 = t267 * t289;
	t283 = t270 * t296;
	t280 = t267 * t286 - t296;
	t277 = -(t263 * t284 - t266 * t273) * t274 - (t263 * t287 + t266 * t293) * t271 + t263 * t283;
	t276 = -(t263 * t273 + t266 * t284) * t274 - (-t263 * t293 + t266 * t287) * t271 + t266 * t283;
	t272 = cos(qJ(4));
	t269 = sin(qJ(4));
	t261 = t274 * t296 - t290;
	t256 = t264 * t289 + (t270 * t291 + t273 * t271) * t265;
	t255 = t266 * t295 + (-t263 * t271 + t266 * t286) * t264;
	t254 = t264 * t271 * t266 + (t264 * t286 + t295) * t263;
	t1 = [t269 * t254 + t277 * t272, t272 * t254 - t277 * t269, (t280 * t263 + t266 * t292) * t273 - t270 * (t263 * t288 - t266 * t274), 0 + (t278 * t271 + t281 * t274 + pkin(1)) * t266 + (-t299 * t268 + t300) * t263; -t269 * t255 - t276 * t272, -t272 * t255 + t276 * t269, (t263 * t292 - t280 * t266) * t273 + (t263 * t274 + t266 * t288) * t270, ((t263 * pkin(3) - pkin(10) * t285) * t273 + (pkin(3) * t285 + t263 * pkin(10)) * t270 - t294 * t298 + t263 * pkin(2)) * t274 + ((pkin(3) * t294 + pkin(10) * t297) * t273 + (-pkin(3) * t297 + pkin(10) * t294) * t270 + pkin(2) * t294 + t263 * t298) * t271 + t263 * pkin(1) + 0 - t266 * t300; t256 * t272 - t269 * t261, -t256 * t269 - t272 * t261, -t264 * t287 + (t270 * t271 - t273 * t291) * t265, t299 * t265 + t279 * t268 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:28
	% EndTime: 2020-11-04 21:17:29
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (133->62), mult. (362->109), div. (0->0), fcn. (417->12), ass. (0->45)
	t319 = sin(qJ(4));
	t322 = cos(qJ(4));
	t370 = t319 * pkin(4) - qJ(5) * t322 + pkin(9);
	t317 = cos(pkin(7));
	t369 = -t370 * t317 - pkin(8);
	t368 = pkin(4) * t322 + qJ(5) * t319 + pkin(3);
	t314 = sin(pkin(7));
	t320 = sin(qJ(3));
	t323 = cos(qJ(3));
	t361 = t323 * pkin(10);
	t367 = (t320 * t368 - t361) * t317 - t370 * t314;
	t318 = cos(pkin(6));
	t321 = sin(qJ(2));
	t349 = t318 * t321;
	t324 = cos(qJ(2));
	t347 = t318 * t324;
	t328 = pkin(10) * t320 + t323 * t368 + pkin(2);
	t315 = sin(pkin(6));
	t359 = t314 * t315;
	t364 = t369 * t315 + t359 * t361;
	t358 = t314 * t318;
	t356 = t315 * t317;
	t355 = t315 * t321;
	t354 = t315 * t324;
	t353 = t317 * t320;
	t352 = t317 * t321;
	t351 = t317 * t324;
	t350 = t318 * t320;
	t348 = t318 * t323;
	t344 = t320 * t359;
	t343 = t317 * t350;
	t313 = sin(pkin(12));
	t336 = t313 * t344;
	t316 = cos(pkin(12));
	t335 = t316 * t344;
	t332 = t317 * t347 - t359;
	t330 = -(t313 * t343 - t316 * t323) * t324 - (t313 * t348 + t316 * t353) * t321 + t336;
	t329 = -(t313 * t323 + t316 * t343) * t324 - (-t313 * t353 + t316 * t348) * t321 + t335;
	t327 = t328 * t316;
	t325 = t367 * t313;
	t308 = t314 * t354 - t318 * t317;
	t303 = t314 * t350 + (t320 * t351 + t323 * t321) * t315;
	t302 = t316 * t356 + (-t313 * t321 + t316 * t347) * t314;
	t301 = t314 * t321 * t316 + (t314 * t347 + t356) * t313;
	t1 = [(t332 * t313 + t316 * t352) * t323 - t320 * (t313 * t349 - t316 * t324), -t319 * t301 - t330 * t322, -t322 * t301 + t330 * t319, (-t318 * t325 + t327) * t324 + t368 * t336 + 0 + (-t321 * t367 + pkin(1)) * t316 + (-t328 * t349 - t364) * t313; (t313 * t352 - t332 * t316) * t323 + (t313 * t324 + t316 * t349) * t320, t319 * t302 + t329 * t322, t322 * t302 - t329 * t319, (t318 * t327 - t325) * t321 - t368 * t335 + 0 + (t328 * t324 + pkin(1)) * t313 + (t367 * t347 + t364) * t316; -t314 * t348 + (t320 * t321 - t323 * t351) * t315, -t303 * t322 + t319 * t308, t303 * t319 + t322 * t308, t367 * t354 + (pkin(10) * t355 + t358 * t368) * t320 + (-pkin(10) * t358 + t355 * t368) * t323 + pkin(2) * t355 + qJ(1) + 0 - t369 * t318; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:29
	% EndTime: 2020-11-04 21:17:29
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (201->80), mult. (465->140), div. (0->0), fcn. (554->14), ass. (0->60)
	t388 = cos(pkin(7));
	t391 = sin(qJ(4));
	t399 = pkin(4) + pkin(11);
	t424 = t399 * t391;
	t395 = cos(qJ(4));
	t427 = t395 * t388;
	t449 = qJ(5) * t427 - t388 * t424;
	t448 = qJ(5) * t391 + t399 * t395 + pkin(3);
	t385 = sin(pkin(7));
	t392 = sin(qJ(3));
	t396 = cos(qJ(3));
	t398 = pkin(5) + pkin(10);
	t425 = t398 * t396;
	t447 = (t392 * t448 - t425) * t388 + t385 * (qJ(5) * t395 - pkin(9) - t424);
	t389 = cos(pkin(6));
	t397 = cos(qJ(2));
	t431 = t389 * t397;
	t386 = sin(pkin(6));
	t442 = t385 * t386;
	t445 = t449 * t386 + t425 * t442;
	t430 = t391 * t392;
	t378 = t385 * t430 - t427;
	t377 = t385 * t395 + t388 * t430;
	t393 = sin(qJ(2));
	t428 = t393 * t396;
	t405 = t377 * t397 + t391 * t428;
	t444 = -t386 * t378 + t405 * t389;
	t403 = t398 * t392 + t396 * t448 + pkin(2);
	t441 = t385 * t389;
	t387 = cos(pkin(12));
	t439 = t386 * t387;
	t438 = t386 * t388;
	t437 = t386 * t393;
	t436 = t386 * t397;
	t435 = t388 * t392;
	t434 = t388 * t393;
	t433 = t389 * t392;
	t432 = t389 * t396;
	t429 = t392 * t393;
	t426 = t396 * t397;
	t421 = t392 * t442;
	t420 = t388 * t433;
	t419 = t388 * t432;
	t417 = t389 * t429;
	t384 = sin(pkin(12));
	t412 = t384 * t421;
	t410 = t387 * t421;
	t406 = -t377 * t393 + t391 * t426;
	t402 = t403 * t387;
	t400 = t447 * t384;
	t394 = cos(qJ(6));
	t390 = sin(qJ(6));
	t381 = t388 * pkin(9) + pkin(8);
	t376 = t385 * t432 + (t388 * t426 - t429) * t386;
	t375 = t389 * t378 + t405 * t386;
	t374 = (t384 * t419 + t387 * t392) * t397 + (-t384 * t442 + t387 * t434) * t396 - t384 * t417;
	t373 = (-t384 * t392 + t387 * t419) * t397 + (-t384 * t434 - t385 * t439) * t396 - t387 * t417;
	t372 = -t444 * t384 + t406 * t387;
	t371 = t406 * t384 + t444 * t387;
	t1 = [t372 * t390 + t394 * t374, t372 * t394 - t390 * t374, ((-t384 * t420 + t387 * t396) * t397 + (-t384 * t432 - t387 * t435) * t393 + t412) * t395 + t391 * (t385 * t393 * t387 + (t385 * t431 + t438) * t384), (-t389 * t400 + t402) * t397 + t448 * t412 + 0 + (-t393 * t447 + pkin(1)) * t387 + (-t403 * t389 * t393 + t386 * t381 - t445) * t384; t371 * t390 - t373 * t394, t371 * t394 + t390 * t373, ((t384 * t396 + t387 * t420) * t397 + (-t384 * t435 + t387 * t432) * t393 - t410) * t395 - t391 * (t387 * t438 + (-t384 * t393 + t387 * t431) * t385), (t389 * t402 - t400) * t393 - t448 * t410 - t381 * t439 + 0 + (t403 * t397 + pkin(1)) * t384 + (t447 * t431 + t445) * t387; t375 * t390 - t394 * t376, t375 * t394 + t390 * t376, (t385 * t433 + (t397 * t435 + t428) * t386) * t395 - t391 * (t385 * t436 - t389 * t388), t447 * t436 + (t398 * t437 + t441 * t448) * t392 + (-t398 * t441 + t437 * t448) * t396 + pkin(2) * t437 + qJ(1) + 0 + (t381 - t449) * t389; 0, 0, 0, 1;];
	Tc_mdh = t1;
end