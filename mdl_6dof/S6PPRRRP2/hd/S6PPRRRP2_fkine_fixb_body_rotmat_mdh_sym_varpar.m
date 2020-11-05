% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPRRRP2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:55
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PPRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:14
	% EndTime: 2020-11-04 20:55:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:14
	% EndTime: 2020-11-04 20:55:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t151 = cos(pkin(11));
	t150 = sin(pkin(11));
	t1 = [t151, -t150, 0, 0; t150, t151, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:14
	% EndTime: 2020-11-04 20:55:14
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t153 = sin(pkin(11));
	t154 = sin(pkin(6));
	t161 = t153 * t154;
	t157 = cos(pkin(6));
	t160 = t153 * t157;
	t156 = cos(pkin(11));
	t159 = t156 * t154;
	t158 = t156 * t157;
	t155 = cos(pkin(12));
	t152 = sin(pkin(12));
	t1 = [-t152 * t160 + t156 * t155, -t156 * t152 - t155 * t160, t161, t156 * pkin(1) + qJ(2) * t161 + 0; t152 * t158 + t153 * t155, -t153 * t152 + t155 * t158, -t159, t153 * pkin(1) - qJ(2) * t159 + 0; t154 * t152, t154 * t155, t157, t157 * qJ(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:14
	% EndTime: 2020-11-04 20:55:14
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (44->30), mult. (101->54), div. (0->0), fcn. (135->10), ass. (0->29)
	t173 = sin(pkin(7));
	t190 = pkin(8) * t173;
	t171 = sin(pkin(12));
	t176 = cos(pkin(11));
	t189 = t171 * t176;
	t172 = sin(pkin(11));
	t188 = t172 * t171;
	t178 = cos(pkin(6));
	t187 = t173 * t178;
	t174 = sin(pkin(6));
	t186 = t174 * t173;
	t175 = cos(pkin(12));
	t177 = cos(pkin(7));
	t185 = t175 * t177;
	t184 = t176 * t178;
	t183 = t177 * t174;
	t182 = t178 * t177;
	t168 = -t171 * pkin(2) + t175 * t190;
	t169 = t177 * pkin(8) + qJ(2);
	t181 = t168 * t178 + t174 * t169;
	t180 = cos(qJ(3));
	t179 = sin(qJ(3));
	t167 = t175 * pkin(2) + t171 * t190 + pkin(1);
	t166 = t176 * t175 - t178 * t188;
	t165 = t171 * t184 + t172 * t175;
	t164 = t175 * t182 - t186;
	t163 = -t164 * t172 - t177 * t189;
	t162 = t164 * t176 - t177 * t188;
	t1 = [t163 * t179 + t166 * t180, t163 * t180 - t166 * t179, (t175 * t187 + t183) * t172 + t173 * t189, t167 * t176 + t181 * t172 + 0; t162 * t179 + t165 * t180, t162 * t180 - t165 * t179, (-t175 * t184 + t188) * t173 - t176 * t183, t167 * t172 - t181 * t176 + 0; t179 * t187 + (t171 * t180 + t179 * t185) * t174, t180 * t187 + (-t171 * t179 + t180 * t185) * t174, -t175 * t186 + t182, -t168 * t174 + t169 * t178 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:14
	% EndTime: 2020-11-04 20:55:14
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (91->53), mult. (217->95), div. (0->0), fcn. (272->12), ass. (0->43)
	t211 = sin(pkin(12));
	t213 = sin(pkin(7));
	t236 = t211 * t213;
	t217 = cos(pkin(7));
	t235 = t211 * t217;
	t218 = cos(pkin(6));
	t234 = t211 * t218;
	t214 = sin(pkin(6));
	t233 = t214 * t213;
	t215 = cos(pkin(12));
	t232 = t215 * t217;
	t231 = t217 * t214;
	t230 = t218 * t213;
	t229 = t218 * t217;
	t198 = t215 * t229 - t233;
	t212 = sin(pkin(11));
	t216 = cos(pkin(11));
	t193 = t198 * t212 + t216 * t235;
	t200 = -t212 * t234 + t216 * t215;
	t220 = sin(qJ(3));
	t222 = cos(qJ(3));
	t228 = t193 * t220 - t200 * t222;
	t194 = -t198 * t216 + t212 * t235;
	t199 = t212 * t215 + t216 * t234;
	t227 = t194 * t220 - t199 * t222;
	t202 = t215 * t213 * pkin(8) - t211 * pkin(2);
	t208 = t217 * pkin(8) + qJ(2);
	t226 = t202 * t218 + t214 * t208;
	t206 = pkin(3) * t232 + t211 * pkin(9);
	t225 = pkin(3) * t233 - t206 * t218;
	t205 = -t211 * pkin(3) + pkin(9) * t232;
	t224 = pkin(9) * t233 - t205 * t218;
	t223 = t214 * t211 * t222 + (t215 * t231 + t230) * t220;
	t221 = cos(qJ(4));
	t219 = sin(qJ(4));
	t204 = pkin(3) * t235 - t215 * pkin(9);
	t203 = t215 * pkin(3) + pkin(9) * t235;
	t201 = t215 * pkin(2) + pkin(8) * t236 + pkin(1);
	t196 = t215 * t233 - t229;
	t195 = t215 * t230 + t231;
	t192 = t195 * t216 - t212 * t236;
	t191 = t195 * t212 + t216 * t236;
	t1 = [t191 * t219 - t228 * t221, t191 * t221 + t228 * t219, t193 * t222 + t200 * t220, (t216 * t203 - t224 * t212) * t222 + (-t216 * t204 + t225 * t212) * t220 + t226 * t212 + t201 * t216 + 0; -t219 * t192 - t227 * t221, -t221 * t192 + t227 * t219, t194 * t222 + t199 * t220, (t212 * t203 + t224 * t216) * t222 + (-t212 * t204 - t225 * t216) * t220 - t226 * t216 + t201 * t212 + 0; -t219 * t196 + t223 * t221, -t221 * t196 - t223 * t219, -t222 * t230 + (t211 * t220 - t222 * t232) * t214, (-pkin(9) * t230 - t205 * t214) * t222 + (pkin(3) * t230 + t206 * t214) * t220 - t202 * t214 + t208 * t218 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:14
	% EndTime: 2020-11-04 20:55:14
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (167->68), mult. (411->113), div. (0->0), fcn. (536->14), ass. (0->54)
	t266 = sin(pkin(12));
	t268 = sin(pkin(7));
	t293 = t266 * t268;
	t272 = cos(pkin(7));
	t292 = t266 * t272;
	t273 = cos(pkin(6));
	t291 = t266 * t273;
	t269 = sin(pkin(6));
	t290 = t269 * t268;
	t270 = cos(pkin(12));
	t289 = t270 * t272;
	t288 = t272 * t269;
	t287 = t273 * t268;
	t286 = t273 * t272;
	t253 = t270 * t286 - t290;
	t267 = sin(pkin(11));
	t271 = cos(pkin(11));
	t248 = t253 * t267 + t271 * t292;
	t255 = -t267 * t291 + t271 * t270;
	t276 = sin(qJ(3));
	t279 = cos(qJ(3));
	t285 = t248 * t276 - t255 * t279;
	t249 = -t253 * t271 + t267 * t292;
	t254 = t267 * t270 + t271 * t291;
	t284 = t249 * t276 - t254 * t279;
	t257 = t270 * t268 * pkin(8) - t266 * pkin(2);
	t263 = t272 * pkin(8) + qJ(2);
	t283 = t257 * t273 + t269 * t263;
	t261 = pkin(3) * t289 + t266 * pkin(9);
	t282 = pkin(3) * t290 - t261 * t273;
	t260 = -t266 * pkin(3) + pkin(9) * t289;
	t281 = pkin(9) * t290 - t260 * t273;
	t280 = t269 * t266 * t279 + (t270 * t288 + t287) * t276;
	t278 = cos(qJ(4));
	t277 = cos(qJ(5));
	t275 = sin(qJ(4));
	t274 = sin(qJ(5));
	t259 = pkin(3) * t292 - t270 * pkin(9);
	t258 = t270 * pkin(3) + pkin(9) * t292;
	t256 = t270 * pkin(2) + pkin(8) * t293 + pkin(1);
	t251 = t270 * t290 - t286;
	t250 = t270 * t287 + t288;
	t247 = t250 * t271 - t267 * t293;
	t246 = t250 * t267 + t271 * t293;
	t245 = -t279 * t287 + (t266 * t276 - t279 * t289) * t269;
	t244 = -t278 * t251 - t280 * t275;
	t243 = -t275 * t251 + t280 * t278;
	t242 = t249 * t279 + t254 * t276;
	t241 = t248 * t279 + t255 * t276;
	t240 = -t278 * t247 + t284 * t275;
	t239 = t246 * t275 - t285 * t278;
	t238 = -t275 * t247 - t284 * t278;
	t237 = t246 * t278 + t285 * t275;
	t1 = [t239 * t277 + t241 * t274, -t239 * t274 + t241 * t277, -t237, t239 * pkin(4) - t237 * pkin(10) + (t271 * t258 - t281 * t267) * t279 + (-t271 * t259 + t282 * t267) * t276 + t283 * t267 + t256 * t271 + 0; t238 * t277 + t242 * t274, -t238 * t274 + t242 * t277, -t240, t238 * pkin(4) - t240 * pkin(10) + (t267 * t258 + t281 * t271) * t279 + (-t267 * t259 - t282 * t271) * t276 - t283 * t271 + t256 * t267 + 0; t243 * t277 + t245 * t274, -t243 * t274 + t245 * t277, -t244, t243 * pkin(4) - t244 * pkin(10) + (-pkin(9) * t287 - t260 * t269) * t279 + (pkin(3) * t287 + t261 * t269) * t276 - t257 * t269 + t263 * t273 + 0 + qJ(1); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:55:14
	% EndTime: 2020-11-04 20:55:14
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (211->80), mult. (507->131), div. (0->0), fcn. (634->14), ass. (0->60)
	t326 = cos(pkin(12));
	t324 = sin(pkin(7));
	t329 = cos(pkin(6));
	t347 = t329 * t324;
	t325 = sin(pkin(6));
	t328 = cos(pkin(7));
	t348 = t328 * t325;
	t305 = t326 * t347 + t348;
	t323 = sin(pkin(11));
	t327 = cos(pkin(11));
	t322 = sin(pkin(12));
	t354 = t322 * t324;
	t299 = t305 * t323 + t327 * t354;
	t334 = cos(qJ(4));
	t356 = t299 * t334;
	t330 = sin(qJ(5));
	t333 = cos(qJ(5));
	t316 = pkin(5) * t333 + qJ(6) * t330 + pkin(4);
	t331 = sin(qJ(4));
	t355 = t316 * t331;
	t353 = t322 * t325;
	t352 = t322 * t328;
	t351 = t322 * t329;
	t350 = t325 * t324;
	t349 = t326 * t328;
	t346 = t329 * t328;
	t300 = t305 * t327 - t323 * t354;
	t345 = t331 * t300;
	t344 = t334 * t300;
	t343 = pkin(5) * t330 - qJ(6) * t333;
	t342 = pkin(10) * t331 + t316 * t334;
	t308 = t326 * t346 - t350;
	t301 = t308 * t323 + t327 * t352;
	t304 = t323 * t351 - t327 * t326;
	t332 = sin(qJ(3));
	t335 = cos(qJ(3));
	t341 = -t301 * t332 - t304 * t335;
	t302 = t308 * t327 - t323 * t352;
	t309 = t323 * t326 + t327 * t351;
	t340 = t302 * t332 + t309 * t335;
	t311 = pkin(8) * t324 * t326 - t322 * pkin(2);
	t319 = pkin(8) * t328 + qJ(2);
	t339 = t311 * t329 + t325 * t319;
	t315 = pkin(3) * t349 + pkin(9) * t322;
	t338 = pkin(3) * t350 - t315 * t329;
	t314 = -t322 * pkin(3) + pkin(9) * t349;
	t337 = pkin(9) * t350 - t314 * t329;
	t307 = t326 * t348 + t347;
	t336 = t307 * t332 + t335 * t353;
	t313 = pkin(3) * t352 - pkin(9) * t326;
	t312 = pkin(3) * t326 + pkin(9) * t352;
	t310 = pkin(2) * t326 + pkin(8) * t354 + pkin(1);
	t306 = t326 * t350 - t346;
	t303 = t307 * t335 - t332 * t353;
	t298 = t302 * t335 - t309 * t332;
	t297 = t301 * t335 - t304 * t332;
	t296 = -t331 * t306 + t334 * t336;
	t295 = t299 * t331 + t334 * t341;
	t294 = t334 * t340 - t345;
	t1 = [t295 * t333 + t297 * t330, t331 * t341 - t356, t295 * t330 - t297 * t333, (-t301 * t342 - t304 * t343 - t327 * t313 + t323 * t338) * t332 + (t301 * t343 - t304 * t342 + t327 * t312 - t323 * t337) * t335 - pkin(10) * t356 + t299 * t355 + t339 * t323 + t310 * t327 + 0; t294 * t333 - t298 * t330, t331 * t340 + t344, t294 * t330 + t298 * t333, (t302 * t342 + t309 * t343 - t323 * t313 - t327 * t338) * t332 + (-t302 * t343 + t309 * t342 + t323 * t312 + t327 * t337) * t335 + pkin(10) * t344 - t316 * t345 - t339 * t327 + t310 * t323 + 0; t296 * t333 - t303 * t330, t334 * t306 + t331 * t336, t296 * t330 + t303 * t333, qJ(1) + 0 + (t319 + (pkin(3) * t332 - pkin(9) * t335) * t324) * t329 + (t332 * t342 - t335 * t343) * t307 + (pkin(10) * t334 - t355) * t306 + (-t314 * t335 + t315 * t332 - t311 + (t332 * t343 + t335 * t342) * t322) * t325; 0, 0, 0, 1;];
	Tc_mdh = t1;
end