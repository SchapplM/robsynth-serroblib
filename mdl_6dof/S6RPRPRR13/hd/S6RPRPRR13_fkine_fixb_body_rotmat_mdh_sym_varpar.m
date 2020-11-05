% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR13 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:43
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRR13_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:29
	% EndTime: 2020-11-04 21:43:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:29
	% EndTime: 2020-11-04 21:43:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t195 = cos(qJ(1));
	t194 = sin(qJ(1));
	t1 = [t195, -t194, 0, 0; t194, t195, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:29
	% EndTime: 2020-11-04 21:43:29
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t196 = sin(pkin(12));
	t200 = sin(qJ(1));
	t207 = t200 * t196;
	t197 = sin(pkin(6));
	t206 = t200 * t197;
	t198 = cos(pkin(12));
	t205 = t200 * t198;
	t201 = cos(qJ(1));
	t204 = t201 * t196;
	t203 = t201 * t197;
	t202 = t201 * t198;
	t199 = cos(pkin(6));
	t1 = [-t199 * t207 + t202, -t199 * t205 - t204, t206, t201 * pkin(1) + qJ(2) * t206 + 0; t199 * t204 + t205, t199 * t202 - t207, -t203, t200 * pkin(1) - qJ(2) * t203 + 0; t197 * t196, t197 * t198, t199, t199 * qJ(2) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:29
	% EndTime: 2020-11-04 21:43:29
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (101->58), div. (0->0), fcn. (135->10), ass. (0->31)
	t216 = sin(pkin(7));
	t239 = pkin(9) * t216;
	t215 = sin(pkin(12));
	t221 = sin(qJ(3));
	t238 = t215 * t221;
	t223 = cos(qJ(3));
	t237 = t215 * t223;
	t217 = sin(pkin(6));
	t236 = t217 * t216;
	t219 = cos(pkin(7));
	t235 = t217 * t219;
	t220 = cos(pkin(6));
	t234 = t220 * t219;
	t233 = t220 * t221;
	t232 = t220 * t223;
	t218 = cos(pkin(12));
	t231 = t221 * t218;
	t222 = sin(qJ(1));
	t230 = t222 * t215;
	t229 = t223 * t218;
	t224 = cos(qJ(1));
	t228 = t224 * t215;
	t227 = t224 * t218;
	t211 = -t215 * pkin(2) + t218 * t239;
	t213 = t219 * pkin(9) + qJ(2);
	t226 = t211 * t220 + t217 * t213;
	t208 = t218 * t234 - t236;
	t225 = t208 * t221 + t215 * t232;
	t210 = t218 * pkin(2) + t215 * t239 + pkin(1);
	t209 = t219 * t238 - t229;
	t1 = [-t224 * t209 - t225 * t222, (-t208 * t222 - t219 * t228) * t223 + t221 * (t220 * t230 - t227), (t222 * t220 * t218 + t228) * t216 + t222 * t235, t210 * t224 + t226 * t222 + 0; -t222 * t209 + t225 * t224, (t208 * t223 - t215 * t233) * t224 - t222 * (t219 * t237 + t231), -(t220 * t227 - t230) * t216 - t224 * t235, t210 * t222 - t226 * t224 + 0; t216 * t233 + (t219 * t231 + t237) * t217, t216 * t232 + (t219 * t229 - t238) * t217, -t218 * t236 + t234, -t211 * t217 + t213 * t220 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:29
	% EndTime: 2020-11-04 21:43:29
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (70->46), mult. (159->83), div. (0->0), fcn. (193->10), ass. (0->36)
	t250 = sin(pkin(12));
	t253 = cos(pkin(12));
	t254 = cos(pkin(7));
	t272 = t253 * t254;
	t243 = pkin(3) * t272 + qJ(4) * t250;
	t244 = -t250 * pkin(3) + qJ(4) * t272;
	t251 = sin(pkin(7));
	t279 = t251 * pkin(9);
	t245 = -t250 * pkin(2) + t253 * t279;
	t247 = t254 * pkin(9) + qJ(2);
	t252 = sin(pkin(6));
	t255 = cos(pkin(6));
	t256 = sin(qJ(3));
	t258 = cos(qJ(3));
	t274 = t252 * t251;
	t280 = (pkin(3) * t274 - t243 * t255) * t256 - (qJ(4) * t274 - t244 * t255) * t258 + t245 * t255 + t252 * t247;
	t278 = t250 * t254;
	t277 = t250 * t256;
	t276 = t250 * t258;
	t275 = t251 * t255;
	t273 = t252 * t254;
	t271 = t255 * t254;
	t270 = t255 * t256;
	t269 = t255 * t258;
	t268 = t256 * t253;
	t257 = sin(qJ(1));
	t267 = t257 * t250;
	t266 = t258 * t253;
	t259 = cos(qJ(1));
	t265 = t259 * t250;
	t264 = t259 * t253;
	t241 = t253 * t271 - t274;
	t260 = t241 * t256 + t250 * t269;
	t242 = t254 * t277 - t266;
	t240 = (pkin(3) * t253 + qJ(4) * t278) * t258 + (-pkin(3) * t278 + qJ(4) * t253) * t256 + t250 * t279 + t253 * pkin(2) + pkin(1);
	t1 = [(t257 * t255 * t253 + t265) * t251 + t257 * t273, t259 * t242 + t260 * t257, (t241 * t257 + t254 * t265) * t258 - t256 * (t255 * t267 - t264), t240 * t259 + t280 * t257 + 0; -(t255 * t264 - t267) * t251 - t259 * t273, t257 * t242 - t260 * t259, (-t241 * t258 + t250 * t270) * t259 + t257 * (t254 * t276 + t268), t240 * t257 - t280 * t259 + 0; -t253 * t274 + t271, -t251 * t270 + (-t254 * t268 - t276) * t252, -t251 * t269 + (-t254 * t266 + t277) * t252, (-qJ(4) * t275 - t244 * t252) * t258 + (pkin(3) * t275 + t243 * t252) * t256 - t245 * t252 + t247 * t255 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:29
	% EndTime: 2020-11-04 21:43:29
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (112->53), mult. (213->90), div. (0->0), fcn. (268->12), ass. (0->44)
	t301 = cos(pkin(12));
	t302 = cos(pkin(7));
	t303 = cos(pkin(6));
	t320 = t303 * t302;
	t299 = sin(pkin(7));
	t300 = sin(pkin(6));
	t323 = t300 * t299;
	t288 = t301 * t320 - t323;
	t298 = sin(pkin(12));
	t308 = cos(qJ(3));
	t305 = sin(qJ(3));
	t319 = t303 * t305;
	t284 = -t288 * t308 + t298 * t319;
	t321 = t303 * t299;
	t287 = t300 * t302 + t301 * t321;
	t304 = sin(qJ(5));
	t307 = cos(qJ(5));
	t331 = t284 * t304 - t307 * t287;
	t311 = pkin(3) + pkin(10);
	t316 = t311 * t301;
	t329 = qJ(4) * t298;
	t291 = t302 * t316 + t329;
	t322 = t301 * t302;
	t325 = t298 * t311;
	t292 = qJ(4) * t322 - t325;
	t310 = pkin(4) + pkin(9);
	t324 = t299 * t310;
	t293 = -t298 * pkin(2) + t301 * t324;
	t294 = t310 * t302 + qJ(2);
	t330 = (-t291 * t303 + t311 * t323) * t305 - (qJ(4) * t323 - t292 * t303) * t308 + t293 * t303 + t294 * t300;
	t328 = t298 * t299;
	t327 = t298 * t305;
	t326 = t298 * t308;
	t318 = t305 * t301;
	t313 = t288 * t305 + t303 * t326;
	t309 = cos(qJ(1));
	t306 = sin(qJ(1));
	t290 = t302 * t326 + t318;
	t289 = -t308 * t301 + t302 * t327;
	t286 = t301 * t323 - t320;
	t285 = t290 * t304 + t307 * t328;
	t282 = (-t300 * t322 - t321) * t308 + t300 * t327;
	t281 = (qJ(4) * t301 - t302 * t325) * t305 + (t302 * t329 + t316) * t308 + t301 * pkin(2) + t298 * t324 + pkin(1);
	t1 = [t309 * t285 - t331 * t306, (-t284 * t306 + t309 * t290) * t307 - (t287 * t306 + t309 * t328) * t304, -t309 * t289 - t313 * t306, t281 * t309 + t330 * t306 + 0; t306 * t285 + t331 * t309, (t284 * t307 + t304 * t287) * t309 + (t290 * t307 - t304 * t328) * t306, -t306 * t289 + t313 * t309, t281 * t306 - t330 * t309 + 0; t282 * t304 - t307 * t286, t282 * t307 + t304 * t286, t299 * t319 + (t302 * t318 + t326) * t300, (t291 * t300 + t311 * t321) * t305 + (-qJ(4) * t321 - t292 * t300) * t308 - t293 * t300 + t294 * t303 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:43:29
	% EndTime: 2020-11-04 21:43:30
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (188->71), mult. (421->114), div. (0->0), fcn. (522->14), ass. (0->54)
	t356 = cos(pkin(12));
	t357 = cos(pkin(7));
	t358 = cos(pkin(6));
	t386 = t358 * t357;
	t354 = sin(pkin(7));
	t355 = sin(pkin(6));
	t389 = t355 * t354;
	t342 = t356 * t386 - t389;
	t367 = pkin(4) + pkin(9);
	t349 = t367 * t357 + qJ(2);
	t365 = cos(qJ(3));
	t353 = sin(pkin(12));
	t361 = sin(qJ(3));
	t368 = pkin(3) + pkin(10);
	t388 = t356 * t357;
	t360 = sin(qJ(5));
	t392 = t353 * t360;
	t364 = cos(qJ(5));
	t394 = pkin(11) * t364;
	t370 = (pkin(5) * t392 + t368 * t388 + (qJ(4) - t394) * t353) * t361 - (qJ(4) * t388 - t353 * t368) * t365 - t354 * t367 * t356 + t353 * pkin(2);
	t377 = -pkin(5) * t360 + t394;
	t387 = t358 * t354;
	t340 = t355 * t357 + t356 * t387;
	t381 = t364 * t340;
	t382 = t361 * t368;
	t385 = t360 * t340;
	t398 = pkin(5) * t381 + pkin(11) * t385 + t349 * t355 + (-qJ(4) * t389 - t377 * t342) * t365 + t382 * t389 - t370 * t358;
	t341 = t355 * t388 + t387;
	t393 = t341 * t365;
	t391 = t353 * t361;
	t390 = t353 * t365;
	t384 = t360 * t365;
	t383 = t361 * t356;
	t380 = t358 * t391;
	t378 = -pkin(5) * t364 - pkin(11) * t360;
	t334 = t342 * t384 - t360 * t380 + t381;
	t338 = t342 * t361 + t358 * t390;
	t359 = sin(qJ(6));
	t363 = cos(qJ(6));
	t376 = t334 * t359 + t363 * t338;
	t375 = qJ(4) - t377;
	t374 = -t355 * t391 + t393;
	t373 = -t342 * t365 + t380;
	t366 = cos(qJ(1));
	t362 = sin(qJ(1));
	t344 = t357 * t390 + t383;
	t343 = -t365 * t356 + t357 * t391;
	t339 = t356 * t389 - t386;
	t337 = t361 * t341 + t355 * t390;
	t336 = t360 * t383 + (t354 * t364 + t357 * t384) * t353;
	t335 = t336 * t359 + t363 * t343;
	t333 = t364 * t339 + t374 * t360;
	t332 = pkin(1) + (t375 * t361 + t368 * t365 + pkin(2)) * t356 + ((t367 - t378) * t354 + (t375 * t365 - t382) * t357) * t353;
	t1 = [(t334 * t362 + t366 * t336) * t363 - t359 * (t338 * t362 + t366 * t343), -t335 * t366 - t376 * t362, (-t366 * t344 + t373 * t362) * t364 + (t354 * t353 * t366 + t340 * t362) * t360, t332 * t366 + t398 * t362 + 0; (-t334 * t363 + t359 * t338) * t366 + (t336 * t363 - t359 * t343) * t362, -t335 * t362 + t376 * t366, (-t373 * t364 - t385) * t366 - (t344 * t364 - t354 * t392) * t362, t332 * t362 - t398 * t366 + 0; -t333 * t363 + t337 * t359, t333 * t359 + t337 * t363, -t360 * t339 + t374 * t364, pkin(8) + 0 + t377 * t393 + (t349 + (-qJ(4) * t365 + t382) * t354) * t358 + t378 * t339 + t370 * t355; 0, 0, 0, 1;];
	Tc_mdh = t1;
end