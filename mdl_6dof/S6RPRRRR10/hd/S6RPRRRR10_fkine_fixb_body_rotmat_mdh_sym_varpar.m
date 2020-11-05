% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRR10 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:57
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:57:15
	% EndTime: 2020-11-04 21:57:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:57:15
	% EndTime: 2020-11-04 21:57:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t142 = cos(qJ(1));
	t141 = sin(qJ(1));
	t1 = [t142, -t141, 0, 0; t141, t142, 0, 0; 0, 0, 1, pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:57:15
	% EndTime: 2020-11-04 21:57:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t143 = sin(pkin(13));
	t147 = sin(qJ(1));
	t154 = t147 * t143;
	t144 = sin(pkin(6));
	t153 = t147 * t144;
	t145 = cos(pkin(13));
	t152 = t147 * t145;
	t148 = cos(qJ(1));
	t151 = t148 * t143;
	t150 = t148 * t144;
	t149 = t148 * t145;
	t146 = cos(pkin(6));
	t1 = [-t146 * t154 + t149, -t146 * t152 - t151, t153, t148 * pkin(1) + qJ(2) * t153 + 0; t146 * t151 + t152, t146 * t149 - t154, -t150, t147 * pkin(1) - qJ(2) * t150 + 0; t144 * t143, t144 * t145, t146, t146 * qJ(2) + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:57:15
	% EndTime: 2020-11-04 21:57:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->32), mult. (101->58), div. (0->0), fcn. (135->10), ass. (0->31)
	t163 = sin(pkin(7));
	t186 = pkin(9) * t163;
	t162 = sin(pkin(13));
	t168 = sin(qJ(3));
	t185 = t162 * t168;
	t170 = cos(qJ(3));
	t184 = t162 * t170;
	t164 = sin(pkin(6));
	t183 = t164 * t163;
	t166 = cos(pkin(7));
	t182 = t164 * t166;
	t167 = cos(pkin(6));
	t181 = t167 * t166;
	t180 = t167 * t168;
	t179 = t167 * t170;
	t165 = cos(pkin(13));
	t178 = t168 * t165;
	t169 = sin(qJ(1));
	t177 = t169 * t162;
	t176 = t170 * t165;
	t171 = cos(qJ(1));
	t175 = t171 * t162;
	t174 = t171 * t165;
	t158 = -t162 * pkin(2) + t165 * t186;
	t160 = t166 * pkin(9) + qJ(2);
	t173 = t158 * t167 + t164 * t160;
	t155 = t165 * t181 - t183;
	t172 = t155 * t168 + t162 * t179;
	t157 = t165 * pkin(2) + t162 * t186 + pkin(1);
	t156 = t166 * t185 - t176;
	t1 = [-t171 * t156 - t172 * t169, (-t155 * t169 - t166 * t175) * t170 + t168 * (t167 * t177 - t174), (t169 * t167 * t165 + t175) * t163 + t169 * t182, t157 * t171 + t173 * t169 + 0; -t169 * t156 + t172 * t171, (t155 * t170 - t162 * t180) * t171 - t169 * (t166 * t184 + t178), -(t167 * t174 - t177) * t163 - t171 * t182, t157 * t169 - t173 * t171 + 0; t163 * t180 + (t166 * t178 + t184) * t164, t163 * t179 + (t166 * t176 - t185) * t164, -t165 * t183 + t181, -t158 * t164 + t160 * t167 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:57:15
	% EndTime: 2020-11-04 21:57:15
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (91->53), mult. (213->97), div. (0->0), fcn. (268->12), ass. (0->38)
	t206 = cos(pkin(13));
	t207 = cos(pkin(7));
	t208 = cos(pkin(6));
	t222 = t208 * t207;
	t204 = sin(pkin(7));
	t205 = sin(pkin(6));
	t225 = t205 * t204;
	t194 = t206 * t222 - t225;
	t203 = sin(pkin(13));
	t210 = sin(qJ(3));
	t213 = cos(qJ(3));
	t221 = t208 * t213;
	t188 = t194 * t210 + t203 * t221;
	t223 = t208 * t204;
	t192 = t205 * t207 + t206 * t223;
	t209 = sin(qJ(4));
	t212 = cos(qJ(4));
	t232 = t188 * t209 + t212 * t192;
	t196 = t206 * t204 * pkin(9) - t203 * pkin(2);
	t224 = t206 * t207;
	t197 = -t203 * pkin(3) + pkin(10) * t224;
	t198 = pkin(3) * t224 + t203 * pkin(10);
	t200 = t207 * pkin(9) + qJ(2);
	t231 = (pkin(3) * t225 - t198 * t208) * t210 - (pkin(10) * t225 - t197 * t208) * t213 + t196 * t208 + t205 * t200;
	t230 = t203 * t204;
	t229 = t203 * t207;
	t228 = t203 * t210;
	t227 = t203 * t213;
	t214 = cos(qJ(1));
	t226 = t203 * t214;
	t219 = t213 * t206;
	t215 = (t205 * t224 + t223) * t210 + t205 * t227;
	t211 = sin(qJ(1));
	t195 = -t207 * t228 + t219;
	t191 = t206 * t225 - t222;
	t190 = t195 * t209 - t212 * t230;
	t187 = (t206 * pkin(3) + pkin(10) * t229) * t213 + (-pkin(3) * t229 + t206 * pkin(10)) * t210 + pkin(9) * t230 + t206 * pkin(2) + pkin(1);
	t1 = [(-t188 * t211 + t214 * t195) * t212 + t209 * (t192 * t211 + t204 * t226), -t214 * t190 + t232 * t211, (t194 * t211 + t207 * t226) * t213 - t210 * (t211 * t208 * t203 - t214 * t206), t187 * t214 + t231 * t211 + 0; (t188 * t212 - t209 * t192) * t214 + (t195 * t212 + t209 * t230) * t211, -t190 * t211 - t232 * t214, (-t194 * t213 + t208 * t228) * t214 + t211 * (t210 * t206 + t207 * t227), t187 * t211 - t231 * t214 + 0; -t209 * t191 + t215 * t212, -t212 * t191 - t215 * t209, -t204 * t221 + (-t207 * t219 + t228) * t205, (-pkin(10) * t223 - t197 * t205) * t213 + (pkin(3) * t223 + t198 * t205) * t210 - t196 * t205 + t200 * t208 + 0 + pkin(8); 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:57:15
	% EndTime: 2020-11-04 21:57:15
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (117->52), mult. (245->81), div. (0->0), fcn. (325->14), ass. (0->46)
	t280 = pkin(4) * sin(qJ(4));
	t254 = sin(pkin(7));
	t279 = pkin(9) * t254;
	t253 = sin(pkin(13));
	t260 = sin(qJ(3));
	t278 = t253 * t260;
	t262 = cos(qJ(3));
	t277 = t253 * t262;
	t255 = sin(pkin(6));
	t276 = t255 * t254;
	t257 = cos(pkin(7));
	t275 = t255 * t257;
	t258 = cos(pkin(6));
	t274 = t258 * t257;
	t273 = t258 * t260;
	t272 = t258 * t262;
	t256 = cos(pkin(13));
	t271 = t260 * t256;
	t261 = sin(qJ(1));
	t270 = t261 * t253;
	t269 = t262 * t256;
	t263 = cos(qJ(1));
	t268 = t263 * t253;
	t267 = t263 * t256;
	t245 = -t253 * pkin(2) + t256 * t279;
	t247 = t257 * pkin(9) + qJ(2);
	t266 = t245 * t258 + t255 * t247;
	t242 = t256 * t274 - t276;
	t265 = t242 * t260 + t253 * t272;
	t264 = -pkin(11) - pkin(10);
	t252 = qJ(4) + qJ(5);
	t250 = cos(t252);
	t249 = sin(t252);
	t248 = cos(qJ(4)) * pkin(4) + pkin(3);
	t244 = t256 * pkin(2) + t253 * t279 + pkin(1);
	t243 = t257 * t278 - t269;
	t241 = -t256 * t276 + t274;
	t240 = -(t258 * t267 - t270) * t254 - t263 * t275;
	t239 = (t261 * t258 * t256 + t268) * t254 + t261 * t275;
	t238 = t254 * t273 + (t257 * t271 + t277) * t255;
	t237 = -t254 * t272 + (-t257 * t269 + t278) * t255;
	t236 = (-t242 * t261 - t257 * t268) * t262 + t260 * (t258 * t270 - t267);
	t235 = -t263 * t243 - t265 * t261;
	t234 = (t242 * t262 - t253 * t273) * t263 - t261 * (t257 * t277 + t271);
	t233 = -t261 * t243 + t265 * t263;
	t1 = [t235 * t250 + t239 * t249, -t235 * t249 + t239 * t250, -t236, t235 * t248 + t236 * t264 + t239 * t280 + t244 * t263 + t266 * t261 + 0; t233 * t250 + t240 * t249, -t233 * t249 + t240 * t250, -t234, t233 * t248 + t234 * t264 + t240 * t280 + t244 * t261 - t266 * t263 + 0; t238 * t250 + t241 * t249, -t238 * t249 + t241 * t250, t237, -t237 * t264 + t238 * t248 + t241 * t280 - t245 * t255 + t247 * t258 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:57:15
	% EndTime: 2020-11-04 21:57:15
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (211->62), mult. (442->99), div. (0->0), fcn. (592->16), ass. (0->54)
	t337 = pkin(4) * sin(qJ(4));
	t309 = sin(pkin(7));
	t336 = pkin(9) * t309;
	t308 = sin(pkin(13));
	t316 = sin(qJ(3));
	t335 = t308 * t316;
	t319 = cos(qJ(3));
	t334 = t308 * t319;
	t310 = sin(pkin(6));
	t333 = t310 * t309;
	t312 = cos(pkin(7));
	t332 = t310 * t312;
	t313 = cos(pkin(6));
	t331 = t313 * t312;
	t330 = t313 * t316;
	t329 = t313 * t319;
	t311 = cos(pkin(13));
	t328 = t316 * t311;
	t317 = sin(qJ(1));
	t327 = t317 * t308;
	t326 = t319 * t311;
	t320 = cos(qJ(1));
	t325 = t320 * t308;
	t324 = t320 * t311;
	t300 = -t308 * pkin(2) + t311 * t336;
	t302 = t312 * pkin(9) + qJ(2);
	t323 = t300 * t313 + t310 * t302;
	t297 = t311 * t331 - t333;
	t322 = t297 * t316 + t308 * t329;
	t321 = -pkin(11) - pkin(10);
	t318 = cos(qJ(6));
	t314 = sin(qJ(6));
	t307 = qJ(4) + qJ(5);
	t305 = cos(t307);
	t304 = sin(t307);
	t303 = cos(qJ(4)) * pkin(4) + pkin(3);
	t299 = t311 * pkin(2) + t308 * t336 + pkin(1);
	t298 = t312 * t335 - t326;
	t296 = -t311 * t333 + t331;
	t294 = -(t313 * t324 - t327) * t309 - t320 * t332;
	t293 = (t317 * t313 * t311 + t325) * t309 + t317 * t332;
	t292 = t309 * t330 + (t312 * t328 + t334) * t310;
	t291 = -t309 * t329 + (-t312 * t326 + t335) * t310;
	t290 = (-t297 * t317 - t312 * t325) * t319 + t316 * (t313 * t327 - t324);
	t289 = -t320 * t298 - t322 * t317;
	t288 = (t297 * t319 - t308 * t330) * t320 - t317 * (t312 * t334 + t328);
	t287 = -t317 * t298 + t322 * t320;
	t286 = t292 * t305 + t296 * t304;
	t285 = t292 * t304 - t296 * t305;
	t284 = t289 * t305 + t293 * t304;
	t283 = t289 * t304 - t293 * t305;
	t282 = t287 * t305 + t294 * t304;
	t281 = t287 * t304 - t294 * t305;
	t1 = [t284 * t318 - t290 * t314, -t284 * t314 - t290 * t318, t283, t284 * pkin(5) + t283 * pkin(12) + t289 * t303 + t290 * t321 + t293 * t337 + t299 * t320 + t323 * t317 + 0; t282 * t318 - t288 * t314, -t282 * t314 - t288 * t318, t281, t282 * pkin(5) + t281 * pkin(12) + t287 * t303 + t288 * t321 + t294 * t337 + t299 * t317 - t323 * t320 + 0; t286 * t318 + t291 * t314, -t286 * t314 + t291 * t318, t285, t286 * pkin(5) + t285 * pkin(12) - t291 * t321 + t292 * t303 + t296 * t337 - t300 * t310 + t302 * t313 + pkin(8) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end