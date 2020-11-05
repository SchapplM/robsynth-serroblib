% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRRPR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:54
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PPRRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:15
	% EndTime: 2020-11-04 20:54:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:15
	% EndTime: 2020-11-04 20:54:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t160 = cos(pkin(11));
	t159 = sin(pkin(11));
	t1 = [t160, -t159, 0, 0; t159, t160, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:15
	% EndTime: 2020-11-04 20:54:15
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t162 = sin(pkin(11));
	t163 = sin(pkin(6));
	t170 = t162 * t163;
	t166 = cos(pkin(6));
	t169 = t162 * t166;
	t165 = cos(pkin(11));
	t168 = t165 * t163;
	t167 = t165 * t166;
	t164 = cos(pkin(12));
	t161 = sin(pkin(12));
	t1 = [-t161 * t169 + t165 * t164, -t165 * t161 - t164 * t169, t170, t165 * pkin(1) + qJ(2) * t170 + 0; t161 * t167 + t162 * t164, -t162 * t161 + t164 * t167, -t168, t162 * pkin(1) - qJ(2) * t168 + 0; t163 * t161, t163 * t164, t166, t166 * qJ(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:15
	% EndTime: 2020-11-04 20:54:15
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (44->30), mult. (101->54), div. (0->0), fcn. (135->10), ass. (0->29)
	t182 = sin(pkin(7));
	t199 = pkin(8) * t182;
	t180 = sin(pkin(12));
	t185 = cos(pkin(11));
	t198 = t180 * t185;
	t181 = sin(pkin(11));
	t197 = t181 * t180;
	t187 = cos(pkin(6));
	t196 = t182 * t187;
	t183 = sin(pkin(6));
	t195 = t183 * t182;
	t184 = cos(pkin(12));
	t186 = cos(pkin(7));
	t194 = t184 * t186;
	t193 = t185 * t187;
	t192 = t186 * t183;
	t191 = t187 * t186;
	t177 = -t180 * pkin(2) + t184 * t199;
	t178 = t186 * pkin(8) + qJ(2);
	t190 = t177 * t187 + t183 * t178;
	t189 = cos(qJ(3));
	t188 = sin(qJ(3));
	t176 = t184 * pkin(2) + t180 * t199 + pkin(1);
	t175 = t185 * t184 - t187 * t197;
	t174 = t180 * t193 + t181 * t184;
	t173 = t184 * t191 - t195;
	t172 = -t173 * t181 - t186 * t198;
	t171 = t173 * t185 - t186 * t197;
	t1 = [t172 * t188 + t175 * t189, t172 * t189 - t175 * t188, (t184 * t196 + t192) * t181 + t182 * t198, t176 * t185 + t190 * t181 + 0; t171 * t188 + t174 * t189, t171 * t189 - t174 * t188, (-t184 * t193 + t197) * t182 - t185 * t192, t176 * t181 - t190 * t185 + 0; t188 * t196 + (t180 * t189 + t188 * t194) * t183, t189 * t196 + (-t180 * t188 + t189 * t194) * t183, -t184 * t195 + t191, -t177 * t183 + t178 * t187 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:15
	% EndTime: 2020-11-04 20:54:15
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (91->53), mult. (217->95), div. (0->0), fcn. (272->12), ass. (0->43)
	t220 = sin(pkin(12));
	t222 = sin(pkin(7));
	t245 = t220 * t222;
	t226 = cos(pkin(7));
	t244 = t220 * t226;
	t227 = cos(pkin(6));
	t243 = t220 * t227;
	t223 = sin(pkin(6));
	t242 = t223 * t222;
	t224 = cos(pkin(12));
	t241 = t224 * t226;
	t240 = t226 * t223;
	t239 = t227 * t222;
	t238 = t227 * t226;
	t207 = t224 * t238 - t242;
	t221 = sin(pkin(11));
	t225 = cos(pkin(11));
	t202 = t207 * t221 + t225 * t244;
	t209 = -t221 * t243 + t225 * t224;
	t229 = sin(qJ(3));
	t231 = cos(qJ(3));
	t237 = t202 * t229 - t209 * t231;
	t203 = -t207 * t225 + t221 * t244;
	t208 = t221 * t224 + t225 * t243;
	t236 = t203 * t229 - t208 * t231;
	t211 = t224 * t222 * pkin(8) - t220 * pkin(2);
	t217 = t226 * pkin(8) + qJ(2);
	t235 = t211 * t227 + t223 * t217;
	t215 = pkin(3) * t241 + t220 * pkin(9);
	t234 = pkin(3) * t242 - t215 * t227;
	t214 = -t220 * pkin(3) + pkin(9) * t241;
	t233 = pkin(9) * t242 - t214 * t227;
	t232 = t223 * t220 * t231 + (t224 * t240 + t239) * t229;
	t230 = cos(qJ(4));
	t228 = sin(qJ(4));
	t213 = pkin(3) * t244 - t224 * pkin(9);
	t212 = t224 * pkin(3) + pkin(9) * t244;
	t210 = t224 * pkin(2) + pkin(8) * t245 + pkin(1);
	t205 = t224 * t242 - t238;
	t204 = t224 * t239 + t240;
	t201 = t204 * t225 - t221 * t245;
	t200 = t204 * t221 + t225 * t245;
	t1 = [t200 * t228 - t237 * t230, t200 * t230 + t237 * t228, t202 * t231 + t209 * t229, (t225 * t212 - t233 * t221) * t231 + (-t225 * t213 + t234 * t221) * t229 + t235 * t221 + t210 * t225 + 0; -t201 * t228 - t236 * t230, -t201 * t230 + t236 * t228, t203 * t231 + t208 * t229, (t221 * t212 + t233 * t225) * t231 + (-t221 * t213 - t234 * t225) * t229 - t235 * t225 + t210 * t221 + 0; -t228 * t205 + t232 * t230, -t230 * t205 - t232 * t228, -t231 * t239 + (t220 * t229 - t231 * t241) * t223, (-pkin(9) * t239 - t214 * t223) * t231 + (pkin(3) * t239 + t215 * t223) * t229 - t211 * t223 + t217 * t227 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:15
	% EndTime: 2020-11-04 20:54:15
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (133->66), mult. (327->109), div. (0->0), fcn. (406->12), ass. (0->49)
	t273 = cos(pkin(12));
	t275 = cos(pkin(7));
	t276 = cos(pkin(6));
	t287 = t276 * t275;
	t271 = sin(pkin(7));
	t272 = sin(pkin(6));
	t291 = t272 * t271;
	t255 = t273 * t287 - t291;
	t270 = sin(pkin(11));
	t274 = cos(pkin(11));
	t269 = sin(pkin(12));
	t293 = t269 * t275;
	t250 = -t255 * t274 + t270 * t293;
	t292 = t269 * t276;
	t256 = t270 * t273 + t274 * t292;
	t278 = sin(qJ(3));
	t280 = cos(qJ(3));
	t303 = t250 * t278 - t256 * t280;
	t277 = sin(qJ(4));
	t279 = cos(qJ(4));
	t286 = pkin(4) * t279 + qJ(5) * t277;
	t288 = t276 * t271;
	t289 = t275 * t272;
	t252 = t273 * t288 + t289;
	t294 = t269 * t271;
	t246 = t252 * t270 + t274 * t294;
	t300 = t246 * t277;
	t299 = t246 * t279;
	t247 = t252 * t274 - t270 * t294;
	t298 = t247 * t277;
	t297 = t247 * t279;
	t296 = (t273 * t289 + t288) * t278;
	t290 = t273 * t275;
	t248 = t255 * t270 + t274 * t293;
	t251 = t270 * t292 - t274 * t273;
	t285 = t248 * t278 + t251 * t280;
	t258 = t273 * t271 * pkin(8) - t269 * pkin(2);
	t266 = t275 * pkin(8) + qJ(2);
	t284 = t258 * t276 + t272 * t266;
	t262 = pkin(3) * t290 + t269 * pkin(9);
	t283 = pkin(3) * t291 - t262 * t276;
	t261 = -t269 * pkin(3) + pkin(9) * t290;
	t282 = pkin(9) * t291 - t261 * t276;
	t281 = t272 * t269 * t280 + t296;
	t260 = pkin(3) * t293 - t273 * pkin(9);
	t259 = t273 * pkin(3) + pkin(9) * t293;
	t257 = t273 * pkin(2) + pkin(8) * t294 + pkin(1);
	t253 = t273 * t291 - t287;
	t1 = [t248 * t280 - t251 * t278, t285 * t279 - t300, -t285 * t277 - t299, (-t286 * t248 - t274 * t260 + t283 * t270) * t278 + (-t286 * t251 + t274 * t259 - t282 * t270) * t280 - qJ(5) * t299 + pkin(4) * t300 + t284 * t270 + t257 * t274 + 0; t250 * t280 + t256 * t278, t303 * t279 + t298, -t303 * t277 + t297, (-t250 * t286 - t270 * t260 - t283 * t274) * t278 + (t286 * t256 + t270 * t259 + t282 * t274) * t280 + qJ(5) * t297 - pkin(4) * t298 - t284 * t274 + t257 * t270 + 0; -t280 * t288 + (t269 * t278 - t280 * t290) * t272, t277 * t253 - t281 * t279, t279 * t253 + t281 * t277, qJ(1) + 0 + t286 * t296 + (t266 + (pkin(3) * t278 - pkin(9) * t280) * t271) * t276 + (-pkin(4) * t277 + qJ(5) * t279) * t253 + (t262 * t278 + (t286 * t269 - t261) * t280 - t258) * t272; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:54:15
	% EndTime: 2020-11-04 20:54:15
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (189->73), mult. (421->120), div. (0->0), fcn. (534->14), ass. (0->61)
	t340 = sin(qJ(4));
	t343 = cos(qJ(4));
	t346 = pkin(4) + pkin(10);
	t353 = qJ(5) * t340 + t346 * t343;
	t334 = sin(pkin(6));
	t335 = cos(pkin(12));
	t337 = cos(pkin(7));
	t333 = sin(pkin(7));
	t338 = cos(pkin(6));
	t359 = t333 * t338;
	t315 = t334 * t337 + t335 * t359;
	t332 = sin(pkin(11));
	t336 = cos(pkin(11));
	t331 = sin(pkin(12));
	t365 = t331 * t333;
	t309 = t315 * t332 + t336 * t365;
	t368 = t309 * t343;
	t310 = t315 * t336 - t332 * t365;
	t367 = t310 * t343;
	t358 = t335 * t337;
	t317 = t334 * t358 + t359;
	t341 = sin(qJ(3));
	t366 = t317 * t341;
	t364 = t331 * t334;
	t363 = t331 * t337;
	t362 = t331 * t338;
	t361 = t333 * t334;
	t360 = t333 * t335;
	t345 = pkin(5) + pkin(9);
	t357 = t337 * t345;
	t356 = t338 * t337;
	t355 = t346 * t340;
	t318 = t335 * t356 - t361;
	t311 = t318 * t332 + t336 * t363;
	t314 = t332 * t362 - t336 * t335;
	t344 = cos(qJ(3));
	t352 = -t311 * t341 - t314 * t344;
	t312 = t318 * t336 - t332 * t363;
	t319 = t332 * t335 + t336 * t362;
	t351 = t312 * t341 + t319 * t344;
	t325 = -t331 * pkin(2) + pkin(8) * t360;
	t328 = pkin(8) * t337 + qJ(2);
	t350 = t325 * t338 + t328 * t334;
	t321 = pkin(3) * t358 + t331 * t345;
	t349 = pkin(3) * t361 - t321 * t338;
	t348 = t344 * t364 + t366;
	t320 = -t331 * pkin(3) + t335 * t357;
	t347 = -t320 * t338 + t345 * t361;
	t342 = cos(qJ(6));
	t339 = sin(qJ(6));
	t324 = pkin(2) * t335 + pkin(8) * t365 + pkin(1);
	t323 = pkin(3) * t363 - t335 * t345;
	t322 = pkin(3) * t335 + t331 * t357;
	t316 = t334 * t360 - t356;
	t313 = t317 * t344 - t341 * t364;
	t308 = t312 * t344 - t319 * t341;
	t307 = t311 * t344 - t314 * t341;
	t306 = t316 * t343 + t348 * t340;
	t305 = t352 * t340 - t368;
	t304 = t351 * t340 + t367;
	t1 = [t305 * t339 + t307 * t342, t305 * t342 - t307 * t339, t309 * t340 + t352 * t343, (-t353 * t311 - t336 * t323 + t349 * t332) * t341 + (-t353 * t314 + t322 * t336 - t347 * t332) * t344 - qJ(5) * t368 + t309 * t355 + t350 * t332 + t324 * t336 + 0; t304 * t339 - t308 * t342, t304 * t342 + t308 * t339, -t310 * t340 + t351 * t343, (t353 * t312 - t332 * t323 - t349 * t336) * t341 + (t353 * t319 + t322 * t332 + t347 * t336) * t344 + qJ(5) * t367 - t310 * t355 - t350 * t336 + t324 * t332 + 0; t306 * t339 - t313 * t342, t306 * t342 + t313 * t339, -t340 * t316 + t348 * t343, qJ(1) + 0 + t353 * t366 + (t328 + (pkin(3) * t341 - t344 * t345) * t333) * t338 + (qJ(5) * t343 - t355) * t316 + (t321 * t341 + (t353 * t331 - t320) * t344 - t325) * t334; 0, 0, 0, 1;];
	Tc_mdh = t1;
end