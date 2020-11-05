% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRP5 (for one body)
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
% Datum: 2020-11-04 21:19
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:25
	% EndTime: 2020-11-04 21:19:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:25
	% EndTime: 2020-11-04 21:19:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t207 = cos(pkin(12));
	t206 = sin(pkin(12));
	t1 = [t207, -t206, 0, 0; t206, t207, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:25
	% EndTime: 2020-11-04 21:19:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t208 = sin(pkin(12));
	t209 = sin(pkin(6));
	t217 = t208 * t209;
	t210 = cos(pkin(12));
	t216 = t210 * t209;
	t211 = cos(pkin(6));
	t212 = sin(qJ(2));
	t215 = t211 * t212;
	t213 = cos(qJ(2));
	t214 = t211 * t213;
	t1 = [-t208 * t215 + t210 * t213, -t208 * t214 - t210 * t212, t217, t210 * pkin(1) + pkin(8) * t217 + 0; t208 * t213 + t210 * t215, -t208 * t212 + t210 * t214, -t216, t208 * pkin(1) - pkin(8) * t216 + 0; t209 * t212, t209 * t213, t211, t211 * pkin(8) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:25
	% EndTime: 2020-11-04 21:19:25
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->35), mult. (112->66), div. (0->0), fcn. (146->10), ass. (0->30)
	t224 = sin(pkin(7));
	t246 = pkin(9) * t224;
	t223 = sin(pkin(12));
	t245 = t223 * pkin(2);
	t226 = cos(pkin(12));
	t244 = t226 * pkin(2);
	t228 = cos(pkin(6));
	t243 = t224 * t228;
	t232 = cos(qJ(2));
	t242 = t224 * t232;
	t227 = cos(pkin(7));
	t222 = t227 * pkin(9) + pkin(8);
	t225 = sin(pkin(6));
	t241 = t225 * t222;
	t240 = t225 * t227;
	t230 = sin(qJ(2));
	t239 = t227 * t230;
	t238 = t227 * t232;
	t237 = t228 * t230;
	t236 = t228 * t232;
	t235 = t223 * t246;
	t234 = t226 * t246;
	t233 = -t224 * t225 + t227 * t236;
	t231 = cos(qJ(3));
	t229 = sin(qJ(3));
	t221 = t223 * t232 + t226 * t237;
	t220 = t223 * t237 - t226 * t232;
	t219 = -t223 * t239 + t233 * t226;
	t218 = -t233 * t223 - t226 * t239;
	t1 = [t218 * t229 - t231 * t220, t218 * t231 + t229 * t220, (t223 * t236 + t226 * t230) * t224 + t223 * t240, (t228 * t235 + t244) * t232 + (-t228 * t245 + t234) * t230 + t223 * t241 + t226 * pkin(1) + 0; t219 * t229 + t221 * t231, t219 * t231 - t221 * t229, -(-t223 * t230 + t226 * t236) * t224 - t226 * t240, (-t228 * t234 + t245) * t232 + (t228 * t244 + t235) * t230 - t226 * t241 + t223 * pkin(1) + 0; t229 * t243 + (t229 * t238 + t230 * t231) * t225, t231 * t243 + (-t229 * t230 + t231 * t238) * t225, -t225 * t242 + t228 * t227, t222 * t228 + qJ(1) + 0 + (pkin(2) * t230 - pkin(9) * t242) * t225; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:25
	% EndTime: 2020-11-04 21:19:26
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (91->59), mult. (247->109), div. (0->0), fcn. (302->12), ass. (0->42)
	t258 = sin(pkin(6));
	t257 = sin(pkin(7));
	t260 = cos(pkin(7));
	t263 = sin(qJ(3));
	t266 = cos(qJ(3));
	t275 = pkin(3) * t263 - pkin(10) * t266;
	t272 = t260 * pkin(9) + t275 * t257 + pkin(8);
	t293 = t272 * t258;
	t264 = sin(qJ(2));
	t267 = cos(qJ(2));
	t291 = t257 * pkin(9);
	t271 = -t275 * t260 + t291;
	t274 = pkin(3) * t266 + pkin(10) * t263 + pkin(2);
	t292 = t274 * t264 - t271 * t267;
	t256 = sin(pkin(12));
	t290 = t256 * t260;
	t289 = t257 * t258;
	t288 = t258 * t260;
	t259 = cos(pkin(12));
	t261 = cos(pkin(6));
	t287 = t259 * t261;
	t286 = t260 * t263;
	t285 = t260 * t264;
	t284 = t260 * t267;
	t283 = t261 * t260;
	t282 = t261 * t263;
	t281 = t261 * t264;
	t280 = t261 * t266;
	t279 = t261 * t267;
	t278 = t259 * t283;
	t277 = t260 * t282;
	t276 = t263 * t289;
	t273 = t260 * t279 - t289;
	t270 = -(t256 * t277 - t259 * t266) * t267 - (t256 * t280 + t259 * t286) * t264 + t256 * t276;
	t269 = -(t256 * t266 + t259 * t277) * t267 - (-t256 * t286 + t259 * t280) * t264 + t259 * t276;
	t265 = cos(qJ(4));
	t262 = sin(qJ(4));
	t254 = t267 * t289 - t283;
	t249 = -t257 * t282 + (-t263 * t284 - t266 * t264) * t258;
	t248 = t259 * t288 + (-t256 * t264 + t259 * t279) * t257;
	t247 = t257 * t264 * t259 + (t257 * t279 + t288) * t256;
	t1 = [t262 * t247 + t270 * t265, t265 * t247 - t270 * t262, (t273 * t256 + t259 * t285) * t266 - t263 * (t256 * t281 - t259 * t267), 0 + (t271 * t264 + t274 * t267 + pkin(1)) * t259 + (-t292 * t261 + t293) * t256; -t262 * t248 - t269 * t265, -t265 * t248 + t269 * t262, (t256 * t285 - t273 * t259) * t266 + (t256 * t267 + t259 * t281) * t263, ((t256 * pkin(3) - pkin(10) * t278) * t266 + (pkin(3) * t278 + t256 * pkin(10)) * t263 - t287 * t291 + t256 * pkin(2)) * t267 + ((pkin(3) * t287 + pkin(10) * t290) * t266 + (-pkin(3) * t290 + pkin(10) * t287) * t263 + pkin(2) * t287 + t256 * t291) * t264 + t256 * pkin(1) + 0 - t259 * t293; -t249 * t265 - t262 * t254, t249 * t262 - t265 * t254, -t257 * t280 + (t263 * t264 - t266 * t284) * t258, t292 * t258 + t272 * t261 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:26
	% EndTime: 2020-11-04 21:19:26
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (167->81), mult. (469->151), div. (0->0), fcn. (558->14), ass. (0->62)
	t316 = cos(pkin(7));
	t323 = cos(qJ(4));
	t353 = t323 * t316;
	t319 = sin(qJ(4));
	t356 = t319 * t316;
	t383 = -pkin(4) * t356 - t316 * pkin(9) + pkin(11) * t353 - pkin(8);
	t382 = t323 * pkin(4) + t319 * pkin(11) + pkin(3);
	t313 = sin(pkin(7));
	t320 = sin(qJ(3));
	t324 = cos(qJ(3));
	t374 = t324 * pkin(10);
	t381 = (t320 * t382 - t374) * t316 - (t319 * pkin(4) - t323 * pkin(11) + pkin(9)) * t313;
	t312 = sin(pkin(12));
	t321 = sin(qJ(2));
	t315 = cos(pkin(12));
	t364 = t315 * t313;
	t314 = sin(pkin(6));
	t367 = t314 * t316;
	t300 = t312 * t367 + t321 * t364;
	t361 = t316 * t321;
	t371 = t313 * t314;
	t303 = -t312 * t371 + t315 * t361;
	t304 = -t313 * t319 + t320 * t353;
	t325 = cos(qJ(2));
	t317 = cos(pkin(6));
	t354 = t321 * t324;
	t346 = t317 * t354;
	t352 = t323 * t324;
	t373 = t312 * t317;
	t378 = t323 * (t303 * t320 + t312 * t346) - (-t304 * t373 + t315 * t352) * t325 - t319 * t300;
	t301 = t312 * t361 + t314 * t364;
	t307 = t315 * t367;
	t363 = t315 * t317;
	t372 = t312 * t321;
	t377 = t323 * (-t301 * t320 + t315 * t346) + (t304 * t363 + t312 * t352) * t325 - t319 * (-t313 * t372 + t307);
	t329 = pkin(10) * t320 + t324 * t382 + pkin(2);
	t376 = t383 * t314 + t371 * t374;
	t370 = t313 * t317;
	t369 = t313 * t320;
	t366 = t314 * t321;
	t365 = t314 * t325;
	t362 = t316 * t320;
	t360 = t316 * t325;
	t359 = t317 * t320;
	t358 = t317 * t324;
	t357 = t317 * t325;
	t355 = t320 * t321;
	t349 = t314 * t369;
	t348 = t316 * t359;
	t347 = t316 * t358;
	t345 = t317 * t355;
	t338 = t312 * t349;
	t337 = t315 * t349;
	t328 = t329 * t315;
	t326 = t381 * t312;
	t322 = cos(qJ(5));
	t318 = sin(qJ(5));
	t299 = t313 * t358 + (t324 * t360 - t355) * t314;
	t296 = -t317 * (t323 * t369 + t356) + (-t304 * t325 - t321 * t352) * t314;
	t295 = (-t312 * t320 + t315 * t347) * t325 - t301 * t324 - t315 * t345;
	t294 = (t312 * t347 + t315 * t320) * t325 + t303 * t324 - t312 * t345;
	t1 = [t318 * t294 - t378 * t322, t322 * t294 + t378 * t318, ((-t312 * t348 + t315 * t324) * t325 + (-t312 * t358 - t315 * t362) * t321 + t338) * t319 - t323 * (t312 * t313 * t357 + t300), (-t317 * t326 + t328) * t325 + (-t315 * t381 - t329 * t373) * t321 + t382 * t338 + t315 * pkin(1) + 0 - t376 * t312; -t318 * t295 + t377 * t322, -t322 * t295 - t377 * t318, ((t312 * t324 + t315 * t348) * t325 + (-t312 * t362 + t315 * t358) * t321 - t337) * t319 + t323 * (t307 + (t315 * t357 - t372) * t313), (t329 * t312 + t381 * t363) * t325 + (t317 * t328 - t326) * t321 - t382 * t337 + t312 * pkin(1) + 0 + t376 * t315; -t296 * t322 - t318 * t299, t296 * t318 - t322 * t299, (t313 * t359 + (t320 * t360 + t354) * t314) * t319 + t323 * (t313 * t365 - t317 * t316), t381 * t365 + (pkin(10) * t366 + t370 * t382) * t320 + (-pkin(10) * t370 + t366 * t382) * t324 + pkin(2) * t366 + qJ(1) + 0 - t383 * t317; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:19:26
	% EndTime: 2020-11-04 21:19:26
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (186->79), mult. (504->131), div. (0->0), fcn. (639->14), ass. (0->58)
	t405 = sin(pkin(6));
	t404 = sin(pkin(7));
	t407 = cos(pkin(7));
	t412 = sin(qJ(3));
	t416 = cos(qJ(3));
	t425 = pkin(3) * t412 - pkin(10) * t416;
	t422 = t407 * pkin(9) + t425 * t404 + pkin(8);
	t446 = t422 * t405;
	t413 = sin(qJ(2));
	t417 = cos(qJ(2));
	t444 = t404 * pkin(9);
	t421 = -t425 * t407 + t444;
	t424 = pkin(3) * t416 + pkin(10) * t412 + pkin(2);
	t445 = t424 * t413 - t421 * t417;
	t403 = sin(pkin(12));
	t406 = cos(pkin(12));
	t408 = cos(pkin(6));
	t429 = t408 * t417;
	t439 = t404 * t405;
	t423 = t407 * t429 - t439;
	t431 = t408 * t413;
	t435 = t407 * t413;
	t388 = (t423 * t403 + t406 * t435) * t416 - t412 * (t403 * t431 - t406 * t417);
	t410 = sin(qJ(5));
	t443 = t388 * t410;
	t389 = (t403 * t435 - t423 * t406) * t416 + (t403 * t417 + t406 * t431) * t412;
	t442 = t389 * t410;
	t430 = t408 * t416;
	t434 = t407 * t417;
	t395 = -t404 * t430 + (t412 * t413 - t416 * t434) * t405;
	t441 = t395 * t410;
	t440 = t403 * t407;
	t438 = t405 * t407;
	t437 = t406 * t408;
	t436 = t407 * t412;
	t433 = t408 * t407;
	t432 = t408 * t412;
	t428 = t406 * t433;
	t427 = t407 * t432;
	t426 = t412 * t439;
	t420 = -(t403 * t427 - t406 * t416) * t417 - (t403 * t430 + t406 * t436) * t413 + t403 * t426;
	t419 = -(t403 * t416 + t406 * t427) * t417 - (-t403 * t436 + t406 * t430) * t413 + t406 * t426;
	t415 = cos(qJ(4));
	t414 = cos(qJ(5));
	t411 = sin(qJ(4));
	t409 = -qJ(6) - pkin(11);
	t402 = t414 * pkin(5) + pkin(4);
	t400 = t417 * t439 - t433;
	t394 = -t404 * t432 + (-t412 * t434 - t416 * t413) * t405;
	t393 = t406 * t438 + (-t403 * t413 + t406 * t429) * t404;
	t392 = t404 * t413 * t406 + (t404 * t429 + t438) * t403;
	t391 = -t394 * t415 - t411 * t400;
	t390 = t394 * t411 - t415 * t400;
	t387 = -t415 * t393 + t419 * t411;
	t386 = -t411 * t393 - t419 * t415;
	t385 = t411 * t392 + t420 * t415;
	t384 = t415 * t392 - t420 * t411;
	t1 = [t385 * t414 + t443, -t385 * t410 + t388 * t414, -t384, pkin(5) * t443 + t384 * t409 + t385 * t402 + 0 + (t421 * t413 + t424 * t417 + pkin(1)) * t406 + (-t445 * t408 + t446) * t403; t386 * t414 + t442, -t386 * t410 + t389 * t414, -t387, t386 * t402 + t387 * t409 + pkin(5) * t442 + ((t403 * pkin(3) - pkin(10) * t428) * t416 + (pkin(3) * t428 + t403 * pkin(10)) * t412 - t437 * t444 + t403 * pkin(2)) * t417 + ((pkin(3) * t437 + pkin(10) * t440) * t416 + (-pkin(3) * t440 + pkin(10) * t437) * t412 + pkin(2) * t437 + t403 * t444) * t413 + t403 * pkin(1) + 0 - t406 * t446; t391 * t414 + t441, -t391 * t410 + t395 * t414, -t390, pkin(5) * t441 + t390 * t409 + t391 * t402 + t445 * t405 + t422 * t408 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end