% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR12 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:50
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:33
	% EndTime: 2020-11-04 21:50:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:33
	% EndTime: 2020-11-04 21:50:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t192 = cos(qJ(1));
	t191 = sin(qJ(1));
	t1 = [t192, -t191, 0, 0; t191, t192, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:33
	% EndTime: 2020-11-04 21:50:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t193 = sin(pkin(12));
	t197 = sin(qJ(1));
	t204 = t197 * t193;
	t194 = sin(pkin(6));
	t203 = t197 * t194;
	t195 = cos(pkin(12));
	t202 = t197 * t195;
	t198 = cos(qJ(1));
	t201 = t198 * t193;
	t200 = t198 * t194;
	t199 = t198 * t195;
	t196 = cos(pkin(6));
	t1 = [-t196 * t204 + t199, -t196 * t202 - t201, t203, t198 * pkin(1) + qJ(2) * t203 + 0; t196 * t201 + t202, t196 * t199 - t204, -t200, t197 * pkin(1) - qJ(2) * t200 + 0; t194 * t193, t194 * t195, t196, t196 * qJ(2) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:33
	% EndTime: 2020-11-04 21:50:33
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (101->58), div. (0->0), fcn. (135->10), ass. (0->31)
	t213 = sin(pkin(7));
	t236 = pkin(9) * t213;
	t212 = sin(pkin(12));
	t218 = sin(qJ(3));
	t235 = t212 * t218;
	t220 = cos(qJ(3));
	t234 = t212 * t220;
	t214 = sin(pkin(6));
	t233 = t214 * t213;
	t216 = cos(pkin(7));
	t232 = t214 * t216;
	t217 = cos(pkin(6));
	t231 = t217 * t216;
	t230 = t217 * t218;
	t229 = t217 * t220;
	t215 = cos(pkin(12));
	t228 = t218 * t215;
	t219 = sin(qJ(1));
	t227 = t219 * t212;
	t226 = t220 * t215;
	t221 = cos(qJ(1));
	t225 = t221 * t212;
	t224 = t221 * t215;
	t208 = -t212 * pkin(2) + t215 * t236;
	t210 = t216 * pkin(9) + qJ(2);
	t223 = t208 * t217 + t214 * t210;
	t205 = t215 * t231 - t233;
	t222 = t205 * t218 + t212 * t229;
	t207 = t215 * pkin(2) + t212 * t236 + pkin(1);
	t206 = t216 * t235 - t226;
	t1 = [-t221 * t206 - t222 * t219, (-t205 * t219 - t216 * t225) * t220 + t218 * (t217 * t227 - t224), (t219 * t217 * t215 + t225) * t213 + t219 * t232, t207 * t221 + t223 * t219 + 0; -t219 * t206 + t222 * t221, (t205 * t220 - t212 * t230) * t221 - t219 * (t216 * t234 + t228), -(t217 * t224 - t227) * t213 - t221 * t232, t207 * t219 - t223 * t221 + 0; t213 * t230 + (t216 * t228 + t234) * t214, t213 * t229 + (t216 * t226 - t235) * t214, -t215 * t233 + t231, -t208 * t214 + t210 * t217 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:33
	% EndTime: 2020-11-04 21:50:33
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (91->53), mult. (213->97), div. (0->0), fcn. (268->12), ass. (0->38)
	t256 = cos(pkin(12));
	t257 = cos(pkin(7));
	t258 = cos(pkin(6));
	t272 = t258 * t257;
	t254 = sin(pkin(7));
	t255 = sin(pkin(6));
	t275 = t255 * t254;
	t244 = t256 * t272 - t275;
	t253 = sin(pkin(12));
	t260 = sin(qJ(3));
	t263 = cos(qJ(3));
	t271 = t258 * t263;
	t238 = t244 * t260 + t253 * t271;
	t273 = t258 * t254;
	t242 = t255 * t257 + t256 * t273;
	t259 = sin(qJ(4));
	t262 = cos(qJ(4));
	t282 = t238 * t259 + t262 * t242;
	t246 = t256 * t254 * pkin(9) - t253 * pkin(2);
	t274 = t256 * t257;
	t247 = -t253 * pkin(3) + pkin(10) * t274;
	t248 = pkin(3) * t274 + t253 * pkin(10);
	t250 = t257 * pkin(9) + qJ(2);
	t281 = (pkin(3) * t275 - t248 * t258) * t260 - (pkin(10) * t275 - t247 * t258) * t263 + t246 * t258 + t255 * t250;
	t280 = t253 * t254;
	t279 = t253 * t257;
	t278 = t253 * t260;
	t277 = t253 * t263;
	t264 = cos(qJ(1));
	t276 = t253 * t264;
	t269 = t263 * t256;
	t265 = (t255 * t274 + t273) * t260 + t255 * t277;
	t261 = sin(qJ(1));
	t245 = -t257 * t278 + t269;
	t241 = t256 * t275 - t272;
	t240 = t245 * t259 - t262 * t280;
	t237 = (t256 * pkin(3) + pkin(10) * t279) * t263 + (-pkin(3) * t279 + t256 * pkin(10)) * t260 + pkin(9) * t280 + t256 * pkin(2) + pkin(1);
	t1 = [(-t238 * t261 + t264 * t245) * t262 + t259 * (t242 * t261 + t254 * t276), -t240 * t264 + t282 * t261, (t244 * t261 + t257 * t276) * t263 - t260 * (t261 * t258 * t253 - t264 * t256), t237 * t264 + t281 * t261 + 0; (t238 * t262 - t259 * t242) * t264 + (t245 * t262 + t259 * t280) * t261, -t240 * t261 - t282 * t264, (-t244 * t263 + t258 * t278) * t264 + t261 * (t260 * t256 + t257 * t277), t237 * t261 - t281 * t264 + 0; -t259 * t241 + t265 * t262, -t262 * t241 - t265 * t259, -t254 * t271 + (-t257 * t269 + t278) * t255, (-pkin(10) * t273 - t247 * t255) * t263 + (pkin(3) * t273 + t248 * t255) * t260 - t246 * t255 + t250 * t258 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:33
	% EndTime: 2020-11-04 21:50:33
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (133->63), mult. (323->106), div. (0->0), fcn. (390->12), ass. (0->44)
	t303 = cos(pkin(12));
	t304 = cos(pkin(7));
	t305 = cos(pkin(6));
	t326 = t305 * t304;
	t301 = sin(pkin(7));
	t302 = sin(pkin(6));
	t330 = t302 * t301;
	t290 = t303 * t326 - t330;
	t300 = sin(pkin(12));
	t307 = sin(qJ(3));
	t310 = cos(qJ(3));
	t324 = t305 * t310;
	t284 = t290 * t307 + t300 * t324;
	t306 = sin(qJ(4));
	t327 = t305 * t301;
	t288 = t302 * t304 + t303 * t327;
	t309 = cos(qJ(4));
	t322 = t309 * t288;
	t340 = t284 * t306 + t322;
	t319 = pkin(4) * t309 + qJ(5) * t306;
	t297 = t304 * pkin(9) + qJ(2);
	t323 = t306 * t288;
	t329 = t302 * t310;
	t339 = -t301 * pkin(10) * t329 + pkin(4) * t323 - qJ(5) * t322 + t302 * t297 + (pkin(3) * t330 - t319 * t290) * t307;
	t328 = t303 * t304;
	t332 = t300 * t309;
	t334 = t300 * t306;
	t313 = (t300 * pkin(3) + pkin(4) * t332 - pkin(10) * t328 + qJ(5) * t334) * t310 - t303 * t301 * pkin(9) + t300 * pkin(2) + (pkin(3) * t328 + t300 * pkin(10)) * t307;
	t337 = pkin(10) * t310;
	t335 = (t302 * t328 + t327) * t307;
	t333 = t300 * t307;
	t311 = cos(qJ(1));
	t331 = t300 * t311;
	t308 = sin(qJ(1));
	t325 = t305 * t308;
	t321 = t310 * t303;
	t318 = -pkin(4) * t306 + qJ(5) * t309;
	t317 = pkin(3) + t319;
	t316 = t300 * t329 + t335;
	t291 = -t304 * t333 + t321;
	t287 = t303 * t330 - t326;
	t286 = t291 * t306 - t301 * t332;
	t283 = pkin(1) + (pkin(10) * t307 + t317 * t310 + pkin(2)) * t303 + ((pkin(9) - t318) * t301 + (-t317 * t307 + t337) * t304) * t300;
	t1 = [(t290 * t308 + t304 * t331) * t310 - t307 * (t300 * t325 - t311 * t303), (t284 * t308 - t311 * t291) * t309 - t306 * (t288 * t308 + t301 * t331), t286 * t311 - t340 * t308, t283 * t311 + t339 * t308 - t313 * t325 + 0; (-t290 * t310 + t305 * t333) * t311 + t308 * (t300 * t310 * t304 + t307 * t303), (-t284 * t309 + t323) * t311 - (t291 * t309 + t301 * t334) * t308, t286 * t308 + t340 * t311, t283 * t308 + 0 + (t313 * t305 - t339) * t311; -t301 * t324 + (-t304 * t321 + t333) * t302, t306 * t287 - t316 * t309, t309 * t287 + t316 * t306, pkin(8) + 0 + t319 * t335 + (t297 + (pkin(3) * t307 - t337) * t301) * t305 + t318 * t287 + t313 * t302; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:50:33
	% EndTime: 2020-11-04 21:50:34
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (195->70), mult. (417->113), div. (0->0), fcn. (518->14), ass. (0->53)
	t364 = sin(pkin(6));
	t365 = cos(pkin(12));
	t366 = cos(pkin(7));
	t363 = sin(pkin(7));
	t367 = cos(pkin(6));
	t395 = t367 * t363;
	t350 = t364 * t366 + t365 * t395;
	t394 = t367 * t366;
	t398 = t364 * t363;
	t352 = t365 * t394 - t398;
	t359 = t366 * pkin(9) + qJ(2);
	t370 = sin(qJ(3));
	t362 = sin(pkin(12));
	t374 = cos(qJ(3));
	t376 = pkin(5) + pkin(10);
	t373 = cos(qJ(4));
	t377 = pkin(4) + pkin(11);
	t392 = t373 * t377;
	t397 = t365 * t366;
	t369 = sin(qJ(4));
	t401 = t362 * t369;
	t379 = (qJ(5) * t401 - t376 * t397 + (pkin(3) + t392) * t362) * t374 + (pkin(3) * t397 + t362 * t376) * t370 - t365 * t363 * pkin(9) + t362 * pkin(2);
	t386 = qJ(5) * t369 + t392;
	t389 = t377 * t369;
	t390 = t374 * t376;
	t393 = t373 * t350;
	t406 = -qJ(5) * t393 + t350 * t389 + t364 * t359 + (pkin(3) * t398 - t386 * t352) * t370 - t390 * t398 - t379 * t367;
	t351 = t364 * t397 + t395;
	t402 = t351 * t370;
	t400 = t362 * t370;
	t399 = t362 * t374;
	t396 = t366 * t370;
	t391 = t374 * t365;
	t387 = qJ(5) * t373 - t389;
	t382 = t352 * t370 + t367 * t399;
	t342 = t382 * t369 + t393;
	t346 = t352 * t374 - t367 * t400;
	t368 = sin(qJ(6));
	t372 = cos(qJ(6));
	t385 = t342 * t368 - t346 * t372;
	t384 = pkin(3) + t386;
	t383 = t364 * t399 + t402;
	t375 = cos(qJ(1));
	t371 = sin(qJ(1));
	t354 = t370 * t365 + t366 * t399;
	t353 = -t362 * t396 + t391;
	t349 = t365 * t398 - t394;
	t347 = t369 * t391 - t362 * (t363 * t373 + t369 * t396);
	t345 = t351 * t374 - t364 * t400;
	t344 = t347 * t368 + t372 * t354;
	t343 = t373 * t349 + t383 * t369;
	t341 = pkin(1) + (t376 * t370 + t384 * t374 + pkin(2)) * t365 + ((pkin(9) - t387) * t363 + (-t384 * t370 + t390) * t366) * t362;
	t1 = [t375 * t344 - t385 * t371, (-t342 * t371 + t375 * t347) * t372 - (t346 * t371 + t375 * t354) * t368, (t375 * t353 - t382 * t371) * t373 + t369 * (t363 * t362 * t375 + t350 * t371), t341 * t375 + t406 * t371 + 0; t371 * t344 + t385 * t375, (t342 * t372 + t346 * t368) * t375 - t371 * (-t347 * t372 + t368 * t354), (-t369 * t350 + t382 * t373) * t375 + (t353 * t373 + t363 * t401) * t371, t341 * t371 - t406 * t375 + 0; t343 * t368 - t345 * t372, t343 * t372 + t368 * t345, -t369 * t349 + t383 * t373, 0 + pkin(8) + t386 * t402 + (t359 + (pkin(3) * t370 - t390) * t363) * t367 + t387 * t349 + t379 * t364; 0, 0, 0, 1;];
	Tc_mdh = t1;
end