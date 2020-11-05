% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:54
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRP12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:38
	% EndTime: 2020-11-04 21:54:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:38
	% EndTime: 2020-11-04 21:54:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t201 = cos(qJ(1));
	t200 = sin(qJ(1));
	t1 = [t201, -t200, 0, 0; t200, t201, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:38
	% EndTime: 2020-11-04 21:54:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t202 = sin(pkin(12));
	t206 = sin(qJ(1));
	t213 = t206 * t202;
	t203 = sin(pkin(6));
	t212 = t206 * t203;
	t204 = cos(pkin(12));
	t211 = t206 * t204;
	t207 = cos(qJ(1));
	t210 = t207 * t202;
	t209 = t207 * t203;
	t208 = t207 * t204;
	t205 = cos(pkin(6));
	t1 = [-t205 * t213 + t208, -t205 * t211 - t210, t212, t207 * pkin(1) + qJ(2) * t212 + 0; t205 * t210 + t211, t205 * t208 - t213, -t209, t206 * pkin(1) - qJ(2) * t209 + 0; t203 * t202, t203 * t204, t205, t205 * qJ(2) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:38
	% EndTime: 2020-11-04 21:54:38
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (101->58), div. (0->0), fcn. (135->10), ass. (0->31)
	t222 = sin(pkin(7));
	t245 = pkin(9) * t222;
	t221 = sin(pkin(12));
	t227 = sin(qJ(3));
	t244 = t221 * t227;
	t229 = cos(qJ(3));
	t243 = t221 * t229;
	t223 = sin(pkin(6));
	t242 = t223 * t222;
	t225 = cos(pkin(7));
	t241 = t223 * t225;
	t226 = cos(pkin(6));
	t240 = t226 * t225;
	t239 = t226 * t227;
	t238 = t226 * t229;
	t224 = cos(pkin(12));
	t237 = t227 * t224;
	t228 = sin(qJ(1));
	t236 = t228 * t221;
	t235 = t229 * t224;
	t230 = cos(qJ(1));
	t234 = t230 * t221;
	t233 = t230 * t224;
	t217 = -t221 * pkin(2) + t224 * t245;
	t219 = t225 * pkin(9) + qJ(2);
	t232 = t217 * t226 + t223 * t219;
	t214 = t224 * t240 - t242;
	t231 = t214 * t227 + t221 * t238;
	t216 = t224 * pkin(2) + t221 * t245 + pkin(1);
	t215 = t225 * t244 - t235;
	t1 = [-t230 * t215 - t231 * t228, (-t214 * t228 - t225 * t234) * t229 + t227 * (t226 * t236 - t233), (t228 * t226 * t224 + t234) * t222 + t228 * t241, t216 * t230 + t232 * t228 + 0; -t228 * t215 + t231 * t230, (t214 * t229 - t221 * t239) * t230 - t228 * (t225 * t243 + t237), -(t226 * t233 - t236) * t222 - t230 * t241, t216 * t228 - t232 * t230 + 0; t222 * t239 + (t225 * t237 + t243) * t223, t222 * t238 + (t225 * t235 - t244) * t223, -t224 * t242 + t240, -t217 * t223 + t219 * t226 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:38
	% EndTime: 2020-11-04 21:54:38
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (91->53), mult. (213->97), div. (0->0), fcn. (268->12), ass. (0->38)
	t265 = cos(pkin(12));
	t266 = cos(pkin(7));
	t267 = cos(pkin(6));
	t281 = t267 * t266;
	t263 = sin(pkin(7));
	t264 = sin(pkin(6));
	t284 = t264 * t263;
	t253 = t265 * t281 - t284;
	t262 = sin(pkin(12));
	t269 = sin(qJ(3));
	t272 = cos(qJ(3));
	t280 = t267 * t272;
	t247 = t253 * t269 + t262 * t280;
	t282 = t267 * t263;
	t251 = t264 * t266 + t265 * t282;
	t268 = sin(qJ(4));
	t271 = cos(qJ(4));
	t291 = t247 * t268 + t271 * t251;
	t255 = t265 * t263 * pkin(9) - t262 * pkin(2);
	t283 = t265 * t266;
	t256 = -t262 * pkin(3) + pkin(10) * t283;
	t257 = pkin(3) * t283 + t262 * pkin(10);
	t259 = t266 * pkin(9) + qJ(2);
	t290 = (pkin(3) * t284 - t257 * t267) * t269 - (pkin(10) * t284 - t256 * t267) * t272 + t255 * t267 + t264 * t259;
	t289 = t262 * t263;
	t288 = t262 * t266;
	t287 = t262 * t269;
	t286 = t262 * t272;
	t273 = cos(qJ(1));
	t285 = t262 * t273;
	t278 = t272 * t265;
	t274 = (t264 * t283 + t282) * t269 + t264 * t286;
	t270 = sin(qJ(1));
	t254 = -t266 * t287 + t278;
	t250 = t265 * t284 - t281;
	t249 = t254 * t268 - t271 * t289;
	t246 = (t265 * pkin(3) + pkin(10) * t288) * t272 + (-pkin(3) * t288 + t265 * pkin(10)) * t269 + pkin(9) * t289 + t265 * pkin(2) + pkin(1);
	t1 = [(-t247 * t270 + t273 * t254) * t271 + t268 * (t251 * t270 + t263 * t285), -t273 * t249 + t291 * t270, (t253 * t270 + t266 * t285) * t272 - t269 * (t270 * t267 * t262 - t273 * t265), t246 * t273 + t290 * t270 + 0; (t247 * t271 - t268 * t251) * t273 + t270 * (t254 * t271 + t268 * t289), -t270 * t249 - t291 * t273, (-t253 * t272 + t267 * t287) * t273 + t270 * (t269 * t265 + t266 * t286), t246 * t270 - t290 * t273 + 0; -t268 * t250 + t274 * t271, -t271 * t250 - t274 * t268, -t263 * t280 + (-t266 * t278 + t287) * t264, (-pkin(10) * t282 - t256 * t264) * t272 + (pkin(3) * t282 + t257 * t264) * t269 - t255 * t264 + t259 * t267 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:38
	% EndTime: 2020-11-04 21:54:38
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (167->66), mult. (417->110), div. (0->0), fcn. (518->14), ass. (0->53)
	t317 = cos(pkin(12));
	t318 = cos(pkin(7));
	t319 = cos(pkin(6));
	t342 = t319 * t318;
	t315 = sin(pkin(7));
	t316 = sin(pkin(6));
	t347 = t316 * t315;
	t305 = t317 * t342 - t347;
	t311 = t318 * pkin(9) + qJ(2);
	t322 = sin(qJ(3));
	t314 = sin(pkin(12));
	t326 = cos(qJ(3));
	t345 = t317 * t318;
	t325 = cos(qJ(4));
	t349 = t314 * t325;
	t321 = sin(qJ(4));
	t352 = pkin(11) * t321;
	t329 = (pkin(4) * t349 - pkin(10) * t345 + (pkin(3) + t352) * t314) * t326 - t317 * t315 * pkin(9) + t314 * pkin(2) + (pkin(3) * t345 + t314 * pkin(10)) * t322;
	t337 = pkin(4) * t325 + t352;
	t343 = t319 * t315;
	t303 = t316 * t318 + t317 * t343;
	t340 = t325 * t303;
	t341 = t321 * t303;
	t346 = t316 * t326;
	t358 = -t315 * pkin(10) * t346 + pkin(4) * t341 - pkin(11) * t340 + t316 * t311 + (pkin(3) * t347 - t337 * t305) * t322 - t319 * t329;
	t348 = t314 * t326;
	t332 = t305 * t322 + t319 * t348;
	t356 = t321 * t332 + t340;
	t353 = pkin(10) * t326;
	t304 = t316 * t345 + t343;
	t351 = t304 * t322;
	t350 = t314 * t322;
	t344 = t318 * t322;
	t339 = t326 * t317;
	t336 = -pkin(4) * t321 + pkin(11) * t325;
	t293 = t325 * t332 - t341;
	t297 = t305 * t326 - t319 * t350;
	t320 = sin(qJ(5));
	t324 = cos(qJ(5));
	t335 = t293 * t320 + t297 * t324;
	t334 = pkin(3) + t337;
	t333 = t314 * t346 + t351;
	t327 = cos(qJ(1));
	t323 = sin(qJ(1));
	t306 = t322 * t317 + t318 * t348;
	t302 = t317 * t347 - t342;
	t299 = t325 * t339 - t314 * (-t315 * t321 + t325 * t344);
	t298 = (-t314 * t344 + t339) * t321 - t315 * t349;
	t296 = t304 * t326 - t316 * t350;
	t295 = -t299 * t320 + t324 * t306;
	t294 = -t321 * t302 + t325 * t333;
	t292 = pkin(1) + (pkin(10) * t322 + t326 * t334 + pkin(2)) * t317 + ((pkin(9) - t336) * t315 + (-t322 * t334 + t353) * t318) * t314;
	t1 = [(-t293 * t323 + t299 * t327) * t324 + t320 * (t297 * t323 + t327 * t306), t327 * t295 + t323 * t335, t327 * t298 - t356 * t323, t292 * t327 + t358 * t323 + 0; (t293 * t324 - t297 * t320) * t327 + t323 * (t299 * t324 + t320 * t306), t323 * t295 - t327 * t335, t323 * t298 + t356 * t327, t292 * t323 - t358 * t327 + 0; t294 * t324 - t320 * t296, -t294 * t320 - t296 * t324, t325 * t302 + t321 * t333, pkin(8) + 0 + t337 * t351 + (t311 + (pkin(3) * t322 - t353) * t315) * t319 + t336 * t302 + t329 * t316; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:54:38
	% EndTime: 2020-11-04 21:54:39
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (223->67), mult. (515->111), div. (0->0), fcn. (622->14), ass. (0->53)
	t390 = sin(qJ(3));
	t394 = cos(qJ(3));
	t388 = sin(qJ(5));
	t392 = cos(qJ(5));
	t408 = pkin(5) * t388 - qJ(6) * t392;
	t403 = pkin(10) + t408;
	t377 = pkin(5) * t392 + qJ(6) * t388 + pkin(4);
	t389 = sin(qJ(4));
	t393 = cos(qJ(4));
	t406 = pkin(11) * t389 + t377 * t393;
	t404 = pkin(3) + t406;
	t429 = t403 * t390 + t404 * t394 + pkin(2);
	t428 = pkin(3) * t390 - pkin(10) * t394;
	t385 = cos(pkin(12));
	t386 = cos(pkin(7));
	t387 = cos(pkin(6));
	t412 = t387 * t386;
	t383 = sin(pkin(7));
	t384 = sin(pkin(6));
	t416 = t384 * t383;
	t372 = t385 * t412 - t416;
	t379 = t386 * pkin(9) + qJ(2);
	t382 = sin(pkin(12));
	t415 = t385 * t386;
	t396 = -t385 * t383 * pkin(9) + t429 * t382 + t428 * t415;
	t413 = t387 * t383;
	t370 = t384 * t386 + t385 * t413;
	t410 = t393 * t370;
	t411 = t389 * t370;
	t427 = (pkin(3) * t416 - t406 * t372) * t390 + (-pkin(10) * t416 + t408 * t372) * t394 - pkin(11) * t410 + t377 * t411 + t384 * t379 - t396 * t387;
	t417 = t382 * t394;
	t401 = t372 * t390 + t387 * t417;
	t425 = t401 * t389 + t410;
	t418 = t382 * t390;
	t414 = t386 * t390;
	t409 = t394 * t385;
	t407 = pkin(11) * t393 - t377 * t389;
	t360 = t401 * t393 - t411;
	t364 = t372 * t394 - t387 * t418;
	t405 = t360 * t388 + t364 * t392;
	t371 = t384 * t415 + t413;
	t402 = t371 * t390 + t384 * t417;
	t395 = cos(qJ(1));
	t391 = sin(qJ(1));
	t373 = t390 * t385 + t386 * t417;
	t369 = t385 * t416 - t412;
	t366 = t393 * t409 - t382 * (-t383 * t389 + t393 * t414);
	t365 = (-t382 * t414 + t409) * t389 - t383 * t382 * t393;
	t363 = t371 * t394 - t384 * t418;
	t362 = -t366 * t388 + t392 * t373;
	t361 = -t389 * t369 + t402 * t393;
	t359 = pkin(1) + t429 * t385 + ((pkin(9) - t407) * t383 + (-t404 * t390 + t403 * t394) * t386) * t382;
	t1 = [(-t360 * t391 + t366 * t395) * t392 + t388 * (t364 * t391 + t395 * t373), t395 * t365 - t425 * t391, -t395 * t362 - t405 * t391, t359 * t395 + t427 * t391 + 0; (t360 * t392 - t364 * t388) * t395 + t391 * (t366 * t392 + t388 * t373), t391 * t365 + t425 * t395, -t391 * t362 + t405 * t395, t359 * t391 - t427 * t395 + 0; t361 * t392 - t388 * t363, t393 * t369 + t402 * t389, t361 * t388 + t363 * t392, pkin(8) + 0 + (t428 * t383 + t379) * t387 + (t406 * t390 - t408 * t394) * t371 + t407 * t369 + t396 * t384; 0, 0, 0, 1;];
	Tc_mdh = t1;
end