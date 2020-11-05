% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPRR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:13
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:01
	% EndTime: 2020-11-04 21:13:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:01
	% EndTime: 2020-11-04 21:13:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t201 = cos(pkin(12));
	t200 = sin(pkin(12));
	t1 = [t201, -t200, 0, 0; t200, t201, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:01
	% EndTime: 2020-11-04 21:13:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t202 = sin(pkin(12));
	t203 = sin(pkin(6));
	t211 = t202 * t203;
	t204 = cos(pkin(12));
	t210 = t204 * t203;
	t205 = cos(pkin(6));
	t206 = sin(qJ(2));
	t209 = t205 * t206;
	t207 = cos(qJ(2));
	t208 = t205 * t207;
	t1 = [-t202 * t209 + t204 * t207, -t202 * t208 - t204 * t206, t211, t204 * pkin(1) + pkin(8) * t211 + 0; t202 * t207 + t204 * t209, -t202 * t206 + t204 * t208, -t210, t202 * pkin(1) - pkin(8) * t210 + 0; t203 * t206, t203 * t207, t205, t205 * pkin(8) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:02
	% EndTime: 2020-11-04 21:13:02
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (44->35), mult. (112->66), div. (0->0), fcn. (146->10), ass. (0->30)
	t218 = sin(pkin(7));
	t240 = pkin(9) * t218;
	t217 = sin(pkin(12));
	t239 = t217 * pkin(2);
	t220 = cos(pkin(12));
	t238 = t220 * pkin(2);
	t222 = cos(pkin(6));
	t237 = t218 * t222;
	t226 = cos(qJ(2));
	t236 = t218 * t226;
	t221 = cos(pkin(7));
	t216 = t221 * pkin(9) + pkin(8);
	t219 = sin(pkin(6));
	t235 = t219 * t216;
	t234 = t219 * t221;
	t224 = sin(qJ(2));
	t233 = t221 * t224;
	t232 = t221 * t226;
	t231 = t222 * t224;
	t230 = t222 * t226;
	t229 = t217 * t240;
	t228 = t220 * t240;
	t227 = -t218 * t219 + t221 * t230;
	t225 = cos(qJ(3));
	t223 = sin(qJ(3));
	t215 = t217 * t226 + t220 * t231;
	t214 = t217 * t231 - t220 * t226;
	t213 = -t217 * t233 + t227 * t220;
	t212 = -t227 * t217 - t220 * t233;
	t1 = [t212 * t223 - t225 * t214, t212 * t225 + t223 * t214, (t217 * t230 + t220 * t224) * t218 + t217 * t234, (t222 * t229 + t238) * t226 + (-t222 * t239 + t228) * t224 + t217 * t235 + t220 * pkin(1) + 0; t213 * t223 + t215 * t225, t213 * t225 - t215 * t223, -(-t217 * t224 + t220 * t230) * t218 - t220 * t234, (-t222 * t228 + t239) * t226 + (t222 * t238 + t229) * t224 - t220 * t235 + t217 * pkin(1) + 0; t223 * t237 + (t223 * t232 + t224 * t225) * t219, t225 * t237 + (-t223 * t224 + t225 * t232) * t219, -t219 * t236 + t222 * t221, t216 * t222 + qJ(1) + 0 + (pkin(2) * t224 - pkin(9) * t236) * t219; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:02
	% EndTime: 2020-11-04 21:13:02
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (70->47), mult. (181->83), div. (0->0), fcn. (215->10), ass. (0->34)
	t248 = sin(pkin(6));
	t247 = sin(pkin(7));
	t250 = cos(pkin(7));
	t252 = sin(qJ(3));
	t254 = cos(qJ(3));
	t261 = pkin(3) * t252 - qJ(4) * t254;
	t258 = t250 * pkin(9) + t261 * t247 + pkin(8);
	t275 = t258 * t248;
	t253 = sin(qJ(2));
	t255 = cos(qJ(2));
	t273 = t247 * pkin(9);
	t257 = -t261 * t250 + t273;
	t260 = pkin(3) * t254 + qJ(4) * t252 + pkin(2);
	t274 = t260 * t253 - t257 * t255;
	t246 = sin(pkin(12));
	t272 = t246 * t250;
	t271 = t247 * t248;
	t251 = cos(pkin(6));
	t270 = t247 * t251;
	t269 = t248 * t250;
	t249 = cos(pkin(12));
	t268 = t249 * t251;
	t267 = t250 * t253;
	t266 = t250 * t255;
	t265 = t251 * t250;
	t264 = t251 * t253;
	t263 = t251 * t255;
	t262 = t249 * t265;
	t259 = t250 * t263 - t271;
	t244 = t246 * t255 + t249 * t264;
	t243 = t246 * t264 - t249 * t255;
	t242 = t246 * t267 - t259 * t249;
	t241 = t259 * t246 + t249 * t267;
	t1 = [(t246 * t263 + t249 * t253) * t247 + t246 * t269, t241 * t252 + t254 * t243, t241 * t254 - t252 * t243, 0 + (t257 * t253 + t260 * t255 + pkin(1)) * t249 + (-t274 * t251 + t275) * t246; -(-t246 * t253 + t249 * t263) * t247 - t249 * t269, t242 * t252 - t244 * t254, t242 * t254 + t244 * t252, ((t246 * pkin(3) - qJ(4) * t262) * t254 + (pkin(3) * t262 + qJ(4) * t246) * t252 - t268 * t273 + t246 * pkin(2)) * t255 + ((pkin(3) * t268 + qJ(4) * t272) * t254 + (-pkin(3) * t272 + qJ(4) * t268) * t252 + pkin(2) * t268 + t246 * t273) * t253 + t246 * pkin(1) + 0 - t249 * t275; -t255 * t271 + t265, -t252 * t270 + (-t252 * t266 - t254 * t253) * t248, -t254 * t270 + (t252 * t253 - t254 * t266) * t248, t274 * t248 + t258 * t251 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:02
	% EndTime: 2020-11-04 21:13:02
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (112->61), mult. (248->104), div. (0->0), fcn. (303->12), ass. (0->47)
	t285 = sin(pkin(6));
	t284 = sin(pkin(7));
	t287 = cos(pkin(7));
	t295 = pkin(4) + pkin(9);
	t290 = sin(qJ(3));
	t293 = cos(qJ(3));
	t296 = pkin(3) + pkin(10);
	t302 = qJ(4) * t293 - t290 * t296;
	t323 = t302 * t284 - t295 * t287 - pkin(8);
	t324 = t285 * t323;
	t291 = sin(qJ(2));
	t294 = cos(qJ(2));
	t318 = t284 * t295;
	t298 = t302 * t287 + t318;
	t301 = qJ(4) * t290 + t296 * t293 + pkin(2);
	t322 = t301 * t291 - t298 * t294;
	t283 = sin(pkin(12));
	t321 = qJ(4) * t283;
	t320 = t283 * t296;
	t319 = t284 * t285;
	t317 = t285 * t287;
	t286 = cos(pkin(12));
	t316 = t286 * t284;
	t288 = cos(pkin(6));
	t315 = t286 * t288;
	t314 = t287 * t291;
	t313 = t287 * t294;
	t312 = t288 * t291;
	t311 = t288 * t293;
	t310 = t288 * t294;
	t309 = t290 * t291;
	t308 = qJ(4) * t315;
	t307 = t283 * t314;
	t306 = t286 * t314;
	t305 = t287 * t311;
	t304 = t296 * t315;
	t303 = t288 * t309;
	t300 = t287 * t310 - t319;
	t292 = cos(qJ(5));
	t289 = sin(qJ(5));
	t281 = -t288 * t287 + t294 * t319;
	t280 = -t284 * t311 + (-t293 * t313 + t309) * t285;
	t279 = t286 * t317 + (-t283 * t291 + t286 * t310) * t284;
	t278 = t291 * t316 + (t284 * t310 + t317) * t283;
	t277 = (t283 * t290 - t286 * t305) * t294 + (t285 * t316 + t307) * t293 + t286 * t303;
	t276 = (t283 * t305 + t286 * t290) * t294 + (-t283 * t319 + t306) * t293 - t283 * t303;
	t1 = [t276 * t289 + t292 * t278, t276 * t292 - t289 * t278, (-t300 * t283 - t306) * t290 - t293 * (t283 * t312 - t286 * t294), 0 + (t298 * t291 + t301 * t294 + pkin(1)) * t286 + (-t322 * t288 - t324) * t283; t277 * t289 - t292 * t279, t277 * t292 + t289 * t279, (t300 * t286 - t307) * t290 + (t283 * t294 + t286 * t312) * t293, ((t287 * t304 + t321) * t290 + (-t287 * t308 + t320) * t293 - t315 * t318 + t283 * pkin(2)) * t294 + ((-t287 * t320 + t308) * t290 + (t287 * t321 + t304) * t293 + pkin(2) * t315 + t283 * t318) * t291 + t283 * pkin(1) + 0 + t286 * t324; t280 * t289 - t292 * t281, t280 * t292 + t289 * t281, t288 * t284 * t290 + (t290 * t313 + t293 * t291) * t285, t322 * t285 - t288 * t323 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:13:02
	% EndTime: 2020-11-04 21:13:03
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (188->84), mult. (475->147), div. (0->0), fcn. (564->14), ass. (0->62)
	t346 = sin(qJ(3));
	t350 = cos(qJ(3));
	t353 = pkin(3) + pkin(10);
	t345 = sin(qJ(5));
	t349 = cos(qJ(5));
	t363 = t345 * pkin(5) - pkin(11) * t349 + qJ(4);
	t401 = -t353 * t346 + t363 * t350;
	t339 = sin(pkin(7));
	t342 = cos(pkin(7));
	t352 = pkin(4) + pkin(9);
	t397 = pkin(11) * t345;
	t355 = t401 * t342 + t339 * (pkin(5) * t349 + t352 + t397);
	t343 = cos(pkin(6));
	t351 = cos(qJ(2));
	t383 = t343 * t351;
	t340 = sin(pkin(6));
	t377 = t342 * t397;
	t395 = t339 * t340;
	t399 = -t340 * t377 + t401 * t395;
	t359 = t363 * t346 + t353 * t350 + pkin(2);
	t394 = t339 * t343;
	t393 = t339 * t349;
	t341 = cos(pkin(12));
	t392 = t340 * t341;
	t391 = t340 * t342;
	t347 = sin(qJ(2));
	t390 = t340 * t347;
	t389 = t340 * t351;
	t388 = t342 * t347;
	t387 = t342 * t349;
	t386 = t342 * t351;
	t385 = t343 * t346;
	t384 = t343 * t350;
	t382 = t345 * t346;
	t381 = t345 * t350;
	t380 = t346 * t347;
	t379 = t347 * t350;
	t376 = t342 * t385;
	t375 = t342 * t384;
	t374 = t343 * t382;
	t373 = t343 * t380;
	t372 = t343 * t379;
	t371 = t340 * t387;
	t338 = sin(pkin(12));
	t366 = t338 * t371;
	t365 = t341 * t371;
	t361 = t342 * t381 + t393;
	t360 = t361 * t343;
	t358 = t359 * t341;
	t334 = -t338 * t395 + t341 * t388;
	t357 = (t338 * t360 + t341 * t382) * t351 - (t338 * t374 - t341 * t393) * t347 + t334 * t381 + t366;
	t333 = t338 * t388 + t339 * t392;
	t356 = (-t338 * t382 + t341 * t360) * t351 - (t338 * t393 + t341 * t374) * t347 - t333 * t381 + t365;
	t354 = t355 * t338;
	t348 = cos(qJ(6));
	t344 = sin(qJ(6));
	t336 = t342 * t352 + pkin(8);
	t330 = t339 * t385 + (t346 * t386 + t379) * t340;
	t327 = (t339 * t381 - t387) * t343 + (-t345 * t380 + t361 * t351) * t340;
	t326 = (t338 * t350 + t341 * t376) * t351 + t341 * t372 - t346 * t333;
	t325 = (t338 * t376 - t341 * t350) * t351 + t338 * t372 + t346 * t334;
	t1 = [-t344 * t325 + t357 * t348, -t348 * t325 - t357 * t344, ((-t338 * t375 - t341 * t346) * t351 - t334 * t350 + t338 * t373) * t349 + t345 * (t341 * t339 * t347 + (t339 * t383 + t391) * t338), (t343 * t354 + t358) * t351 + pkin(5) * t366 + 0 + (t355 * t347 + pkin(1)) * t341 + (-t359 * t343 * t347 + t340 * t336 - t399) * t338; t326 * t344 - t356 * t348, t326 * t348 + t356 * t344, ((-t338 * t346 + t341 * t375) * t351 - t333 * t350 - t341 * t373) * t349 - t345 * (t341 * t391 + (-t338 * t347 + t341 * t383) * t339), (t343 * t358 + t354) * t347 - pkin(5) * t365 - t336 * t392 + 0 + (t359 * t351 + pkin(1)) * t338 + (-t355 * t383 + t399) * t341; -t327 * t348 + t330 * t344, t327 * t344 + t330 * t348, (t339 * t384 + (t350 * t386 - t380) * t340) * t349 - t345 * (t339 * t389 - t342 * t343), -t355 * t389 + (t353 * t390 - t363 * t394) * t350 + (t353 * t394 + t363 * t390) * t346 + pkin(2) * t390 + qJ(1) + 0 + (pkin(5) * t387 + t336 + t377) * t343; 0, 0, 0, 1;];
	Tc_mdh = t1;
end