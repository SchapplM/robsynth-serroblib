% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR11 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:42
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:43
	% EndTime: 2020-11-04 21:42:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:43
	% EndTime: 2020-11-04 21:42:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t149 = cos(qJ(1));
	t148 = sin(qJ(1));
	t1 = [t149, -t148, 0, 0; t148, t149, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:43
	% EndTime: 2020-11-04 21:42:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t150 = sin(pkin(12));
	t154 = sin(qJ(1));
	t161 = t154 * t150;
	t151 = sin(pkin(6));
	t160 = t154 * t151;
	t152 = cos(pkin(12));
	t159 = t154 * t152;
	t155 = cos(qJ(1));
	t158 = t155 * t150;
	t157 = t155 * t151;
	t156 = t155 * t152;
	t153 = cos(pkin(6));
	t1 = [-t153 * t161 + t156, -t153 * t159 - t158, t160, t155 * pkin(1) + qJ(2) * t160 + 0; t153 * t158 + t159, t153 * t156 - t161, -t157, t154 * pkin(1) - qJ(2) * t157 + 0; t151 * t150, t151 * t152, t153, t153 * qJ(2) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:43
	% EndTime: 2020-11-04 21:42:43
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (44->32), mult. (101->58), div. (0->0), fcn. (135->10), ass. (0->31)
	t170 = sin(pkin(7));
	t193 = pkin(9) * t170;
	t169 = sin(pkin(12));
	t175 = sin(qJ(3));
	t192 = t169 * t175;
	t177 = cos(qJ(3));
	t191 = t169 * t177;
	t171 = sin(pkin(6));
	t190 = t171 * t170;
	t173 = cos(pkin(7));
	t189 = t171 * t173;
	t174 = cos(pkin(6));
	t188 = t174 * t173;
	t187 = t174 * t175;
	t186 = t174 * t177;
	t172 = cos(pkin(12));
	t185 = t175 * t172;
	t176 = sin(qJ(1));
	t184 = t176 * t169;
	t183 = t177 * t172;
	t178 = cos(qJ(1));
	t182 = t178 * t169;
	t181 = t178 * t172;
	t165 = -t169 * pkin(2) + t172 * t193;
	t167 = t173 * pkin(9) + qJ(2);
	t180 = t165 * t174 + t167 * t171;
	t162 = t172 * t188 - t190;
	t179 = t162 * t175 + t169 * t186;
	t164 = pkin(2) * t172 + t169 * t193 + pkin(1);
	t163 = t173 * t192 - t183;
	t1 = [-t178 * t163 - t176 * t179, (-t162 * t176 - t173 * t182) * t177 + t175 * (t174 * t184 - t181), (t172 * t174 * t176 + t182) * t170 + t176 * t189, t164 * t178 + t176 * t180 + 0; -t176 * t163 + t178 * t179, (t162 * t177 - t169 * t187) * t178 - t176 * (t173 * t191 + t185), -(t174 * t181 - t184) * t170 - t178 * t189, t164 * t176 - t178 * t180 + 0; t170 * t187 + (t173 * t185 + t191) * t171, t170 * t186 + (t173 * t183 - t192) * t171, -t172 * t190 + t188, -t165 * t171 + t167 * t174 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:43
	% EndTime: 2020-11-04 21:42:43
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (91->52), mult. (223->96), div. (0->0), fcn. (278->12), ass. (0->39)
	t210 = sin(pkin(12));
	t214 = cos(pkin(12));
	t215 = cos(pkin(7));
	t232 = t214 * t215;
	t203 = pkin(3) * t232 + qJ(4) * t210;
	t204 = -t210 * pkin(3) + qJ(4) * t232;
	t211 = sin(pkin(7));
	t205 = t214 * t211 * pkin(9) - t210 * pkin(2);
	t206 = t215 * pkin(9) + qJ(2);
	t212 = sin(pkin(6));
	t216 = cos(pkin(6));
	t217 = sin(qJ(3));
	t219 = cos(qJ(3));
	t233 = t212 * t211;
	t237 = (pkin(3) * t233 - t203 * t216) * t217 - (qJ(4) * t233 - t204 * t216) * t219 + t205 * t216 + t212 * t206;
	t236 = t210 * t211;
	t235 = t210 * t215;
	t234 = t210 * t217;
	t231 = t214 * t219;
	t230 = t215 * t217;
	t229 = t215 * t219;
	t228 = t216 * t211;
	t227 = t216 * t215;
	t226 = t216 * t219;
	t222 = t210 * t212 * t219 + (t212 * t232 + t228) * t217;
	t202 = t214 * t227 - t233;
	t221 = t202 * t217 + t210 * t226;
	t220 = cos(qJ(1));
	t218 = sin(qJ(1));
	t213 = cos(pkin(13));
	t209 = sin(pkin(13));
	t200 = t212 * t215 + t214 * t228;
	t199 = t214 * t233 - t227;
	t198 = t213 * t231 + (t209 * t211 - t213 * t230) * t210;
	t197 = -t213 * t236 + (-t210 * t230 + t231) * t209;
	t196 = (pkin(3) * t214 + qJ(4) * t235) * t219 + (-pkin(3) * t235 + qJ(4) * t214) * t217 + pkin(9) * t236 + t214 * pkin(2) + pkin(1);
	t195 = -t209 * t200 + t221 * t213;
	t194 = t213 * t200 + t221 * t209;
	t1 = [-t195 * t218 + t220 * t198, t194 * t218 - t220 * t197, (t202 * t218 + t220 * t235) * t219 - t217 * (t218 * t216 * t210 - t220 * t214), t196 * t220 + t237 * t218 + 0; t195 * t220 + t218 * t198, -t194 * t220 - t218 * t197, (-t202 * t219 + t216 * t234) * t220 + t218 * (t210 * t229 + t217 * t214), t196 * t218 - t237 * t220 + 0; -t209 * t199 + t222 * t213, -t213 * t199 - t222 * t209, -t211 * t226 + (-t214 * t229 + t234) * t212, (-qJ(4) * t228 - t204 * t212) * t219 + (pkin(3) * t228 + t203 * t212) * t217 - t205 * t212 + t206 * t216 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:43
	% EndTime: 2020-11-04 21:42:43
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (137->55), mult. (234->93), div. (0->0), fcn. (289->14), ass. (0->47)
	t263 = cos(pkin(12));
	t264 = cos(pkin(7));
	t265 = cos(pkin(6));
	t280 = t265 * t264;
	t261 = sin(pkin(7));
	t262 = sin(pkin(6));
	t283 = t262 * t261;
	t249 = t263 * t280 - t283;
	t260 = sin(pkin(12));
	t267 = sin(qJ(3));
	t269 = cos(qJ(3));
	t279 = t265 * t269;
	t241 = t249 * t267 + t260 * t279;
	t276 = t269 * t263;
	t288 = t260 * t267;
	t250 = t264 * t288 - t276;
	t268 = sin(qJ(1));
	t270 = cos(qJ(1));
	t293 = t241 * t268 + t270 * t250;
	t292 = t241 * t270 - t268 * t250;
	t255 = cos(pkin(13)) * pkin(4) + pkin(3);
	t282 = t263 * t264;
	t266 = qJ(4) + pkin(10);
	t289 = t260 * t266;
	t243 = t255 * t282 + t289;
	t278 = t266 * t263;
	t290 = t260 * t255;
	t244 = t264 * t278 - t290;
	t254 = sin(pkin(13)) * pkin(4) + pkin(9);
	t284 = t261 * t254;
	t245 = -t260 * pkin(2) + t263 * t284;
	t251 = t254 * t264 + qJ(2);
	t291 = (-t243 * t265 + t255 * t283) * t267 - (-t244 * t265 + t266 * t283) * t269 + t245 * t265 + t262 * t251;
	t287 = t260 * t268;
	t286 = t260 * t269;
	t285 = t260 * t270;
	t281 = t265 * t261;
	t271 = (t262 * t282 + t281) * t267 + t262 * t286;
	t259 = pkin(13) + qJ(5);
	t257 = cos(t259);
	t256 = sin(t259);
	t247 = t262 * t264 + t263 * t281;
	t246 = t263 * t283 - t280;
	t240 = t247 * t270 - t261 * t287;
	t239 = t247 * t268 + t261 * t285;
	t238 = (t255 * t263 + t264 * t289) * t269 + (-t264 * t290 + t278) * t267 + t263 * pkin(2) + t260 * t284 + pkin(1);
	t1 = [t239 * t256 - t293 * t257, t257 * t239 + t293 * t256, (t249 * t268 + t264 * t285) * t269 - t267 * (-t270 * t263 + t265 * t287), t238 * t270 + t291 * t268 + 0; -t240 * t256 + t292 * t257, -t257 * t240 - t292 * t256, (-t249 * t269 + t265 * t288) * t270 + t268 * (t267 * t263 + t264 * t286), t238 * t268 - t291 * t270 + 0; -t256 * t246 + t271 * t257, -t257 * t246 - t271 * t256, -t261 * t279 + (-t264 * t276 + t288) * t262, (-t244 * t262 - t266 * t281) * t269 + (t243 * t262 + t255 * t281) * t267 - t245 * t262 + t251 * t265 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:42:43
	% EndTime: 2020-11-04 21:42:44
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (231->70), mult. (428->111), div. (0->0), fcn. (553->16), ass. (0->58)
	t328 = cos(pkin(12));
	t329 = cos(pkin(7));
	t330 = cos(pkin(6));
	t347 = t330 * t329;
	t326 = sin(pkin(7));
	t327 = sin(pkin(6));
	t350 = t327 * t326;
	t314 = t328 * t347 - t350;
	t325 = sin(pkin(12));
	t333 = sin(qJ(3));
	t336 = cos(qJ(3));
	t346 = t330 * t336;
	t306 = t314 * t333 + t325 * t346;
	t343 = t336 * t328;
	t355 = t325 * t333;
	t315 = t329 * t355 - t343;
	t334 = sin(qJ(1));
	t337 = cos(qJ(1));
	t360 = t306 * t334 + t337 * t315;
	t359 = t306 * t337 - t334 * t315;
	t320 = cos(pkin(13)) * pkin(4) + pkin(3);
	t349 = t328 * t329;
	t331 = qJ(4) + pkin(10);
	t356 = t325 * t331;
	t308 = t320 * t349 + t356;
	t345 = t331 * t328;
	t357 = t325 * t320;
	t309 = t329 * t345 - t357;
	t319 = sin(pkin(13)) * pkin(4) + pkin(9);
	t351 = t326 * t319;
	t310 = -t325 * pkin(2) + t328 * t351;
	t316 = t319 * t329 + qJ(2);
	t358 = (-t308 * t330 + t320 * t350) * t333 - (-t309 * t330 + t331 * t350) * t336 + t310 * t330 + t327 * t316;
	t354 = t325 * t334;
	t353 = t325 * t336;
	t352 = t325 * t337;
	t348 = t330 * t326;
	t338 = (t327 * t349 + t348) * t333 + t327 * t353;
	t335 = cos(qJ(6));
	t332 = sin(qJ(6));
	t324 = pkin(13) + qJ(5);
	t322 = cos(t324);
	t321 = sin(t324);
	t312 = t327 * t329 + t328 * t348;
	t311 = t328 * t350 - t347;
	t305 = t312 * t337 - t326 * t354;
	t304 = t312 * t334 + t326 * t352;
	t303 = -t326 * t346 + (-t329 * t343 + t355) * t327;
	t302 = (-t314 * t336 + t330 * t355) * t337 + t334 * (t333 * t328 + t329 * t353);
	t301 = (t314 * t334 + t329 * t352) * t336 - t333 * (-t337 * t328 + t330 * t354);
	t300 = -t322 * t311 - t338 * t321;
	t299 = -t321 * t311 + t338 * t322;
	t298 = (t320 * t328 + t329 * t356) * t336 + (-t329 * t357 + t345) * t333 + t328 * pkin(2) + t325 * t351 + pkin(1);
	t297 = -t322 * t305 - t321 * t359;
	t296 = t304 * t321 - t322 * t360;
	t295 = -t305 * t321 + t322 * t359;
	t294 = t322 * t304 + t321 * t360;
	t1 = [t296 * t335 + t301 * t332, -t296 * t332 + t301 * t335, -t294, t296 * pkin(5) - t294 * pkin(11) + t298 * t337 + t358 * t334 + 0; t295 * t335 + t302 * t332, -t295 * t332 + t302 * t335, -t297, t295 * pkin(5) - t297 * pkin(11) + t298 * t334 - t358 * t337 + 0; t299 * t335 + t303 * t332, -t299 * t332 + t303 * t335, -t300, t299 * pkin(5) - t300 * pkin(11) + (-t309 * t327 - t331 * t348) * t336 + (t308 * t327 + t320 * t348) * t333 - t310 * t327 + t316 * t330 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
end