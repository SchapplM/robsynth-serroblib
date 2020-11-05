% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRPR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:17
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRRRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [13x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:01
	% EndTime: 2020-11-04 21:17:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:01
	% EndTime: 2020-11-04 21:17:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t136 = cos(pkin(12));
	t135 = sin(pkin(12));
	t1 = [t136, -t135, 0, 0; t135, t136, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:01
	% EndTime: 2020-11-04 21:17:01
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t137 = sin(pkin(12));
	t138 = sin(pkin(6));
	t146 = t137 * t138;
	t139 = cos(pkin(12));
	t145 = t139 * t138;
	t140 = cos(pkin(6));
	t141 = sin(qJ(2));
	t144 = t140 * t141;
	t142 = cos(qJ(2));
	t143 = t140 * t142;
	t1 = [-t137 * t144 + t139 * t142, -t137 * t143 - t139 * t141, t146, t139 * pkin(1) + pkin(8) * t146 + 0; t137 * t142 + t139 * t144, -t137 * t141 + t139 * t143, -t145, t137 * pkin(1) - pkin(8) * t145 + 0; t138 * t141, t138 * t142, t140, t140 * pkin(8) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:01
	% EndTime: 2020-11-04 21:17:01
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (44->35), mult. (112->66), div. (0->0), fcn. (146->10), ass. (0->30)
	t153 = sin(pkin(7));
	t175 = pkin(9) * t153;
	t152 = sin(pkin(12));
	t174 = t152 * pkin(2);
	t155 = cos(pkin(12));
	t173 = t155 * pkin(2);
	t157 = cos(pkin(6));
	t172 = t153 * t157;
	t161 = cos(qJ(2));
	t171 = t153 * t161;
	t156 = cos(pkin(7));
	t151 = t156 * pkin(9) + pkin(8);
	t154 = sin(pkin(6));
	t170 = t154 * t151;
	t169 = t154 * t156;
	t159 = sin(qJ(2));
	t168 = t156 * t159;
	t167 = t156 * t161;
	t166 = t157 * t159;
	t165 = t157 * t161;
	t164 = t152 * t175;
	t163 = t155 * t175;
	t162 = -t153 * t154 + t156 * t165;
	t160 = cos(qJ(3));
	t158 = sin(qJ(3));
	t150 = t152 * t161 + t155 * t166;
	t149 = t152 * t166 - t155 * t161;
	t148 = -t152 * t168 + t162 * t155;
	t147 = -t162 * t152 - t155 * t168;
	t1 = [t147 * t158 - t160 * t149, t147 * t160 + t158 * t149, (t152 * t165 + t155 * t159) * t153 + t152 * t169, (t157 * t164 + t173) * t161 + (-t157 * t174 + t163) * t159 + t152 * t170 + t155 * pkin(1) + 0; t148 * t158 + t150 * t160, t148 * t160 - t150 * t158, -(-t152 * t159 + t155 * t165) * t153 - t155 * t169, (-t157 * t163 + t174) * t161 + (t157 * t173 + t164) * t159 - t155 * t170 + t152 * pkin(1) + 0; t158 * t172 + (t158 * t167 + t159 * t160) * t154, t160 * t172 + (-t158 * t159 + t160 * t167) * t154, -t154 * t171 + t157 * t156, t151 * t157 + qJ(1) + 0 + (pkin(2) * t159 - pkin(9) * t171) * t154; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:01
	% EndTime: 2020-11-04 21:17:02
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (91->59), mult. (247->109), div. (0->0), fcn. (302->12), ass. (0->42)
	t187 = sin(pkin(6));
	t186 = sin(pkin(7));
	t189 = cos(pkin(7));
	t192 = sin(qJ(3));
	t195 = cos(qJ(3));
	t204 = pkin(3) * t192 - pkin(10) * t195;
	t201 = t189 * pkin(9) + t204 * t186 + pkin(8);
	t222 = t201 * t187;
	t193 = sin(qJ(2));
	t196 = cos(qJ(2));
	t220 = t186 * pkin(9);
	t200 = -t204 * t189 + t220;
	t203 = pkin(3) * t195 + pkin(10) * t192 + pkin(2);
	t221 = t203 * t193 - t200 * t196;
	t185 = sin(pkin(12));
	t219 = t185 * t189;
	t218 = t186 * t187;
	t217 = t187 * t189;
	t188 = cos(pkin(12));
	t190 = cos(pkin(6));
	t216 = t188 * t190;
	t215 = t189 * t192;
	t214 = t189 * t193;
	t213 = t189 * t196;
	t212 = t190 * t189;
	t211 = t190 * t192;
	t210 = t190 * t193;
	t209 = t190 * t195;
	t208 = t190 * t196;
	t207 = t188 * t212;
	t206 = t189 * t211;
	t205 = t192 * t218;
	t202 = t189 * t208 - t218;
	t199 = -(t185 * t206 - t188 * t195) * t196 - (t185 * t209 + t188 * t215) * t193 + t185 * t205;
	t198 = -(t185 * t195 + t188 * t206) * t196 - (-t185 * t215 + t188 * t209) * t193 + t188 * t205;
	t194 = cos(qJ(4));
	t191 = sin(qJ(4));
	t183 = t196 * t218 - t212;
	t178 = t186 * t211 + (t192 * t213 + t195 * t193) * t187;
	t177 = t188 * t217 + (-t185 * t193 + t188 * t208) * t186;
	t176 = t186 * t193 * t188 + (t186 * t208 + t217) * t185;
	t1 = [t191 * t176 + t199 * t194, t194 * t176 - t199 * t191, (t202 * t185 + t188 * t214) * t195 - t192 * (t185 * t210 - t188 * t196), 0 + (t200 * t193 + t203 * t196 + pkin(1)) * t188 + (-t221 * t190 + t222) * t185; -t191 * t177 - t198 * t194, -t194 * t177 + t198 * t191, (t185 * t214 - t202 * t188) * t195 + (t185 * t196 + t188 * t210) * t192, ((t185 * pkin(3) - pkin(10) * t207) * t195 + (pkin(3) * t207 + t185 * pkin(10)) * t192 - t216 * t220 + t185 * pkin(2)) * t196 + ((pkin(3) * t216 + pkin(10) * t219) * t195 + (-pkin(3) * t219 + pkin(10) * t216) * t192 + pkin(2) * t216 + t185 * t220) * t193 + t185 * pkin(1) + 0 - t188 * t222; t178 * t194 - t191 * t183, -t178 * t191 - t194 * t183, -t186 * t209 + (t192 * t193 - t195 * t213) * t187, t221 * t187 + t201 * t190 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:02
	% EndTime: 2020-11-04 21:17:02
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (167->74), mult. (464->127), div. (0->0), fcn. (589->14), ass. (0->53)
	t244 = sin(pkin(6));
	t243 = sin(pkin(7));
	t247 = cos(pkin(7));
	t250 = sin(qJ(3));
	t253 = cos(qJ(3));
	t262 = pkin(3) * t250 - pkin(10) * t253;
	t259 = t247 * pkin(9) + t262 * t243 + pkin(8);
	t280 = t259 * t244;
	t251 = sin(qJ(2));
	t254 = cos(qJ(2));
	t278 = t243 * pkin(9);
	t258 = -t262 * t247 + t278;
	t261 = pkin(3) * t253 + pkin(10) * t250 + pkin(2);
	t279 = t261 * t251 - t258 * t254;
	t242 = sin(pkin(12));
	t277 = t242 * t247;
	t276 = t243 * t244;
	t275 = t244 * t247;
	t246 = cos(pkin(12));
	t248 = cos(pkin(6));
	t274 = t246 * t248;
	t273 = t247 * t250;
	t272 = t247 * t251;
	t271 = t247 * t254;
	t270 = t248 * t247;
	t269 = t248 * t250;
	t268 = t248 * t251;
	t267 = t248 * t253;
	t266 = t248 * t254;
	t265 = t246 * t270;
	t264 = t247 * t269;
	t263 = t250 * t276;
	t260 = t247 * t266 - t276;
	t257 = -(t242 * t264 - t246 * t253) * t254 - (t242 * t267 + t246 * t273) * t251 + t242 * t263;
	t256 = -(t242 * t253 + t246 * t264) * t254 - (-t242 * t273 + t246 * t267) * t251 + t246 * t263;
	t252 = cos(qJ(4));
	t249 = sin(qJ(4));
	t245 = cos(pkin(13));
	t241 = sin(pkin(13));
	t239 = t254 * t276 - t270;
	t234 = t243 * t269 + (t250 * t271 + t253 * t251) * t244;
	t233 = -t243 * t267 + (t250 * t251 - t253 * t271) * t244;
	t232 = t246 * t275 + (-t242 * t251 + t246 * t266) * t243;
	t231 = t243 * t251 * t246 + (t243 * t266 + t275) * t242;
	t230 = -t234 * t249 - t252 * t239;
	t229 = t234 * t252 - t249 * t239;
	t228 = (t242 * t272 - t260 * t246) * t253 + (t242 * t254 + t246 * t268) * t250;
	t227 = (t260 * t242 + t246 * t272) * t253 - t250 * (t242 * t268 - t246 * t254);
	t226 = -t252 * t232 + t256 * t249;
	t225 = -t249 * t232 - t256 * t252;
	t224 = t249 * t231 + t257 * t252;
	t223 = t252 * t231 - t257 * t249;
	t1 = [t224 * t245 + t227 * t241, -t224 * t241 + t227 * t245, -t223, t224 * pkin(4) - t223 * qJ(5) + 0 + (t258 * t251 + t261 * t254 + pkin(1)) * t246 + (-t279 * t248 + t280) * t242; t225 * t245 + t228 * t241, -t225 * t241 + t228 * t245, -t226, t225 * pkin(4) - t226 * qJ(5) + ((t242 * pkin(3) - pkin(10) * t265) * t253 + (pkin(3) * t265 + t242 * pkin(10)) * t250 - t274 * t278 + t242 * pkin(2)) * t254 + ((pkin(3) * t274 + pkin(10) * t277) * t253 + (-pkin(3) * t277 + pkin(10) * t274) * t250 + pkin(2) * t274 + t242 * t278) * t251 + t242 * pkin(1) + 0 - t246 * t280; t229 * t245 + t233 * t241, -t229 * t241 + t233 * t245, -t230, t229 * pkin(4) - t230 * qJ(5) + t279 * t244 + t259 * t248 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:17:02
	% EndTime: 2020-11-04 21:17:02
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (198->80), mult. (504->132), div. (0->0), fcn. (639->16), ass. (0->57)
	t306 = sin(pkin(6));
	t305 = sin(pkin(7));
	t308 = cos(pkin(7));
	t312 = sin(qJ(3));
	t315 = cos(qJ(3));
	t324 = pkin(3) * t312 - pkin(10) * t315;
	t321 = t308 * pkin(9) + t324 * t305 + pkin(8);
	t343 = t321 * t306;
	t313 = sin(qJ(2));
	t316 = cos(qJ(2));
	t340 = t305 * pkin(9);
	t320 = -t324 * t308 + t340;
	t323 = pkin(3) * t315 + pkin(10) * t312 + pkin(2);
	t342 = t323 * t313 - t320 * t316;
	t341 = sin(pkin(13)) * pkin(5);
	t304 = sin(pkin(12));
	t339 = t304 * t308;
	t338 = t305 * t306;
	t337 = t306 * t308;
	t307 = cos(pkin(12));
	t309 = cos(pkin(6));
	t336 = t307 * t309;
	t335 = t308 * t312;
	t334 = t308 * t313;
	t333 = t308 * t316;
	t332 = t309 * t308;
	t331 = t309 * t312;
	t330 = t309 * t313;
	t329 = t309 * t315;
	t328 = t309 * t316;
	t327 = t307 * t332;
	t326 = t308 * t331;
	t325 = t312 * t338;
	t322 = t308 * t328 - t338;
	t319 = -(t304 * t326 - t307 * t315) * t316 - (t304 * t329 + t307 * t335) * t313 + t304 * t325;
	t318 = -(t304 * t315 + t307 * t326) * t316 - (-t304 * t335 + t307 * t329) * t313 + t307 * t325;
	t314 = cos(qJ(4));
	t311 = sin(qJ(4));
	t310 = -pkin(11) - qJ(5);
	t302 = pkin(13) + qJ(6);
	t301 = cos(t302);
	t300 = sin(t302);
	t298 = cos(pkin(13)) * pkin(5) + pkin(4);
	t297 = t316 * t338 - t332;
	t292 = t305 * t331 + (t312 * t333 + t315 * t313) * t306;
	t291 = -t305 * t329 + (t312 * t313 - t315 * t333) * t306;
	t290 = t307 * t337 + (-t304 * t313 + t307 * t328) * t305;
	t289 = t305 * t313 * t307 + (t305 * t328 + t337) * t304;
	t288 = -t292 * t311 - t314 * t297;
	t287 = t292 * t314 - t311 * t297;
	t286 = (t304 * t334 - t322 * t307) * t315 + (t304 * t316 + t307 * t330) * t312;
	t285 = (t322 * t304 + t307 * t334) * t315 - t312 * (t304 * t330 - t307 * t316);
	t284 = -t314 * t290 + t318 * t311;
	t283 = -t311 * t290 - t318 * t314;
	t282 = t311 * t289 + t319 * t314;
	t281 = t314 * t289 - t319 * t311;
	t1 = [t282 * t301 + t285 * t300, -t282 * t300 + t285 * t301, -t281, t285 * t341 + t281 * t310 + t282 * t298 + 0 + (t320 * t313 + t323 * t316 + pkin(1)) * t307 + (-t342 * t309 + t343) * t304; t283 * t301 + t286 * t300, -t283 * t300 + t286 * t301, -t284, t283 * t298 + t284 * t310 + t286 * t341 + ((t304 * pkin(3) - pkin(10) * t327) * t315 + (pkin(3) * t327 + t304 * pkin(10)) * t312 - t336 * t340 + t304 * pkin(2)) * t316 + ((pkin(3) * t336 + pkin(10) * t339) * t315 + (-pkin(3) * t339 + pkin(10) * t336) * t312 + pkin(2) * t336 + t304 * t340) * t313 + t304 * pkin(1) + 0 - t307 * t343; t287 * t301 + t291 * t300, -t287 * t300 + t291 * t301, -t288, t287 * t298 + t288 * t310 + t291 * t341 + t342 * t306 + t321 * t309 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end