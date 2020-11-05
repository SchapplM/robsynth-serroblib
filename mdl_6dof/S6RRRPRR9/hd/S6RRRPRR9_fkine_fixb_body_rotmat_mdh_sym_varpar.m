% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:32
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRRPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:35
	% EndTime: 2020-11-04 22:32:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:35
	% EndTime: 2020-11-04 22:32:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t160 = cos(qJ(1));
	t159 = sin(qJ(1));
	t1 = [t160, -t159, 0, 0; t159, t160, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:35
	% EndTime: 2020-11-04 22:32:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t161 = sin(pkin(6));
	t164 = sin(qJ(1));
	t172 = t164 * t161;
	t163 = sin(qJ(2));
	t171 = t164 * t163;
	t165 = cos(qJ(2));
	t170 = t164 * t165;
	t166 = cos(qJ(1));
	t169 = t166 * t161;
	t168 = t166 * t163;
	t167 = t166 * t165;
	t162 = cos(pkin(6));
	t1 = [-t162 * t171 + t167, -t162 * t170 - t168, t172, t166 * pkin(1) + pkin(9) * t172 + 0; t162 * t168 + t170, t162 * t167 - t171, -t169, t164 * pkin(1) - pkin(9) * t169 + 0; t161 * t163, t161 * t165, t162, t162 * pkin(9) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:35
	% EndTime: 2020-11-04 22:32:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (104->58), div. (0->0), fcn. (138->10), ass. (0->29)
	t179 = sin(pkin(7));
	t182 = cos(pkin(6));
	t201 = t179 * t182;
	t187 = cos(qJ(2));
	t200 = t179 * t187;
	t180 = sin(pkin(6));
	t181 = cos(pkin(7));
	t199 = t180 * t181;
	t198 = t182 * t187;
	t183 = sin(qJ(3));
	t184 = sin(qJ(2));
	t197 = t183 * t184;
	t196 = t183 * t187;
	t186 = cos(qJ(3));
	t195 = t184 * t186;
	t185 = sin(qJ(1));
	t194 = t185 * t184;
	t193 = t186 * t187;
	t188 = cos(qJ(1));
	t192 = t188 * t184;
	t191 = t188 * t187;
	t190 = -pkin(2) * t184 + pkin(10) * t200;
	t174 = -t179 * t180 + t181 * t198;
	t189 = t174 * t183 + t182 * t195;
	t178 = t181 * pkin(10) + pkin(9);
	t176 = t179 * t184 * pkin(10) + pkin(2) * t187 + pkin(1);
	t175 = -t181 * t197 + t193;
	t173 = t180 * t178 + t190 * t182;
	t1 = [t188 * t175 - t189 * t185, (-t174 * t185 - t181 * t192) * t186 - (-t182 * t194 + t191) * t183, (t185 * t198 + t192) * t179 + t185 * t199, t173 * t185 + t176 * t188 + 0; t185 * t175 + t189 * t188, (t174 * t186 - t182 * t197) * t188 - t185 * (t181 * t195 + t196), -(t182 * t191 - t194) * t179 - t188 * t199, -t173 * t188 + t176 * t185 + 0; t183 * t201 + (t181 * t196 + t195) * t180, t186 * t201 + (t181 * t193 - t197) * t180, -t180 * t200 + t182 * t181, t178 * t182 - t190 * t180 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:35
	% EndTime: 2020-11-04 22:32:35
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (131->42), mult. (130->62), div. (0->0), fcn. (164->16), ass. (0->45)
	t223 = qJ(3) + pkin(13);
	t219 = pkin(7) + t223;
	t213 = sin(t219);
	t220 = pkin(7) - t223;
	t214 = sin(t220);
	t209 = -t213 + t214;
	t247 = -t209 / 0.2e1;
	t215 = cos(t219);
	t216 = cos(t220);
	t211 = t215 + t216;
	t246 = t211 / 0.2e1;
	t227 = cos(pkin(6));
	t245 = t227 / 0.2e1;
	t244 = pkin(3) * sin(qJ(3));
	t225 = sin(pkin(6));
	t231 = sin(qJ(1));
	t243 = t225 * t231;
	t233 = cos(qJ(1));
	t242 = t225 * t233;
	t230 = sin(qJ(2));
	t241 = t231 * t230;
	t232 = cos(qJ(2));
	t240 = t231 * t232;
	t239 = t233 * t230;
	t238 = t233 * t232;
	t237 = t243 / 0.2e1;
	t236 = -t242 / 0.2e1;
	t224 = sin(pkin(7));
	t226 = cos(pkin(7));
	t228 = pkin(10) + qJ(4);
	t235 = t224 * t244 + t226 * t228 + pkin(9);
	t208 = t224 * t228 - t226 * t244;
	t218 = cos(qJ(3)) * pkin(3) + pkin(2);
	t234 = t208 * t232 - t218 * t230;
	t222 = cos(t223);
	t221 = sin(t223);
	t212 = t216 - t215;
	t210 = t214 + t213;
	t207 = t227 * t239 + t240;
	t206 = t227 * t238 - t241;
	t205 = t227 * t240 + t239;
	t204 = t227 * t241 - t238;
	t203 = t208 * t230 + t218 * t232 + pkin(1);
	t202 = t225 * t235 + t234 * t227;
	t1 = [-t204 * t222 + t205 * t209 / 0.2e1 + t212 * t237, t204 * t221 - t205 * t211 / 0.2e1 + t210 * t237, t205 * t224 + t226 * t243, t202 * t231 + t203 * t233 + 0; t206 * t247 + t207 * t222 + t212 * t236, t206 * t246 - t207 * t221 + t210 * t236, -t206 * t224 - t226 * t242, -t202 * t233 + t203 * t231 + 0; t212 * t245 + (t230 * t222 + t232 * t247) * t225, t210 * t245 + (-t230 * t221 + t232 * t246) * t225, -t225 * t232 * t224 + t227 * t226, -t234 * t225 + t235 * t227 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:35
	% EndTime: 2020-11-04 22:32:35
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (183->53), mult. (265->84), div. (0->0), fcn. (320->20), ass. (0->51)
	t275 = sin(pkin(13));
	t278 = cos(pkin(13));
	t265 = t278 * pkin(4) + t275 * pkin(11) + pkin(3);
	t266 = -t275 * pkin(4) + t278 * pkin(11);
	t283 = sin(qJ(3));
	t287 = cos(qJ(3));
	t308 = t265 * t283 - t266 * t287;
	t255 = t265 * t287 + t266 * t283 + pkin(2);
	t284 = sin(qJ(2));
	t288 = cos(qJ(2));
	t276 = sin(pkin(7));
	t279 = cos(pkin(7));
	t281 = pkin(10) + qJ(4);
	t292 = -t276 * t281 + t308 * t279;
	t307 = t255 * t284 + t292 * t288;
	t274 = qJ(3) + pkin(13);
	t270 = pkin(7) + t274;
	t271 = pkin(7) - t274;
	t264 = cos(t270) + cos(t271);
	t306 = -t264 / 0.2e1;
	t277 = sin(pkin(6));
	t302 = t276 * t277;
	t263 = sin(t271) + sin(t270);
	t301 = t277 * t263;
	t300 = t279 * t277;
	t280 = cos(pkin(6));
	t299 = t280 * t288;
	t285 = sin(qJ(1));
	t298 = t285 * t284;
	t297 = t285 * t288;
	t289 = cos(qJ(1));
	t296 = t289 * t284;
	t295 = t289 * t288;
	t259 = t279 * t299 - t302;
	t261 = t280 * t298 - t295;
	t272 = sin(t274);
	t273 = cos(t274);
	t294 = (t259 * t285 + t279 * t296) * t272 + t261 * t273;
	t262 = t280 * t296 + t297;
	t293 = (t259 * t289 - t279 * t298) * t272 + t262 * t273;
	t291 = t273 * t277 * t284 + (t280 * t276 + t288 * t300) * t272;
	t290 = t308 * t276 + t279 * t281 + pkin(9);
	t286 = cos(qJ(5));
	t282 = sin(qJ(5));
	t257 = -t280 * t279 + t288 * t302;
	t256 = t276 * t299 + t300;
	t252 = t256 * t289 - t276 * t298;
	t251 = t256 * t285 + t276 * t296;
	t249 = t255 * t288 - t284 * t292 + pkin(1);
	t248 = t290 * t277 - t307 * t280;
	t1 = [t251 * t282 - t294 * t286, t251 * t286 + t294 * t282, -t261 * t272 + (t280 * t297 + t296) * t264 / 0.2e1 - t285 * t301 / 0.2e1, t248 * t285 + t249 * t289 + 0; -t282 * t252 + t293 * t286, -t252 * t286 - t293 * t282, t262 * t272 + (t280 * t295 - t298) * t306 + t289 * t301 / 0.2e1, -t248 * t289 + t249 * t285 + 0; -t282 * t257 + t291 * t286, -t286 * t257 - t291 * t282, -t280 * t263 / 0.2e1 + (t284 * t272 + t288 * t306) * t277, t307 * t277 + t290 * t280 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:32:35
	% EndTime: 2020-11-04 22:32:35
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (258->75), mult. (494->131), div. (0->0), fcn. (583->18), ass. (0->58)
	t340 = sin(qJ(5));
	t345 = cos(qJ(5));
	t370 = -t340 * pkin(5) + pkin(12) * t345 - pkin(10) - qJ(4);
	t333 = sin(pkin(7));
	t332 = sin(pkin(13));
	t335 = cos(pkin(13));
	t350 = pkin(5) * t345 + pkin(12) * t340 + pkin(4);
	t318 = t332 * pkin(11) + t350 * t335 + pkin(3);
	t319 = -t335 * pkin(11) + t350 * t332;
	t341 = sin(qJ(3));
	t346 = cos(qJ(3));
	t351 = t318 * t341 + t319 * t346;
	t369 = -t351 * t333 - pkin(9);
	t334 = sin(pkin(6));
	t368 = t333 * t334;
	t336 = cos(pkin(7));
	t367 = t336 * t334;
	t337 = cos(pkin(6));
	t347 = cos(qJ(2));
	t366 = t337 * t347;
	t339 = sin(qJ(6));
	t365 = t339 * t340;
	t342 = sin(qJ(2));
	t364 = t339 * t342;
	t363 = t339 * t345;
	t344 = cos(qJ(6));
	t362 = t340 * t344;
	t343 = sin(qJ(1));
	t361 = t342 * t343;
	t360 = t342 * t344;
	t348 = cos(qJ(1));
	t359 = t342 * t348;
	t358 = t344 * t345;
	t357 = t347 * t339;
	t356 = t347 * t344;
	t355 = t342 * t363;
	t354 = t342 * t358;
	t311 = t333 * t370 + t351 * t336;
	t312 = t318 * t346 - t319 * t341 + pkin(2);
	t352 = t311 * t347 + t312 * t342;
	t331 = qJ(3) + pkin(13);
	t330 = cos(t331);
	t329 = sin(t331);
	t326 = t337 * t361 - t348 * t347;
	t325 = t336 * t366 - t368;
	t324 = t337 * t333 + t347 * t367;
	t323 = -t337 * t336 + t347 * t368;
	t322 = t333 * t366 + t367;
	t321 = t336 * t360 - t345 * t357;
	t320 = t336 * t355 + t356;
	t317 = t325 * t343 + t336 * t359;
	t316 = t322 * t348 - t333 * t361;
	t315 = t322 * t343 + t333 * t359;
	t314 = t325 * t363 - t337 * t360;
	t313 = t344 * t325 + t337 * t355;
	t310 = -t311 * t342 + t312 * t347 + pkin(1);
	t309 = t369 * t334 + t352 * t337 + t367 * t370;
	t1 = [(t317 * t339 - t326 * t358) * t330 + (-t317 * t358 - t339 * t326) * t329 + t315 * t362, (t313 * t343 + t348 * t321) * t330 + (t314 * t343 + t348 * t320) * t329 - t315 * t365, -t315 * t345 + (-t317 * t329 - t326 * t330) * t340, -t309 * t343 + t310 * t348 + 0; ((-t339 * t325 + t337 * t354) * t348 + t343 * (t336 * t364 + t345 * t356)) * t330 + ((t325 * t358 + t337 * t364) * t348 - t343 * (t336 * t354 - t357)) * t329 - t316 * t362, (-t313 * t348 + t343 * t321) * t330 + (-t314 * t348 + t343 * t320) * t329 + t316 * t365, t316 * t345 + ((t325 * t348 - t336 * t361) * t329 + (t337 * t359 + t343 * t347) * t330) * t340, t309 * t348 + t310 * t343 + 0; (-t339 * t324 + t334 * t354) * t330 + (t324 * t358 + t334 * t364) * t329 - t323 * t362, (-t344 * t324 - t334 * t355) * t330 + (-t324 * t363 + t334 * t360) * t329 + t323 * t365, t345 * t323 + (t330 * t334 * t342 + t324 * t329) * t340, pkin(8) + 0 + (-t336 * t370 - t369) * t337 + t352 * t334; 0, 0, 0, 1;];
	Tc_mdh = t1;
end