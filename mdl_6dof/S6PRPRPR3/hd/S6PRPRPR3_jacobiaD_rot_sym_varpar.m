% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PRPRPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(10));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(10)) * cos(pkin(6));
	t37 = t39 * t42 - t41 * t45;
	t34 = 0.1e1 / t37 ^ 2;
	t49 = qJD(2) * t34;
	t36 = t39 * t41 + t42 * t45;
	t33 = t36 ^ 2;
	t30 = t33 * t34 + 0.1e1;
	t46 = t37 * t49;
	t47 = t36 / t37 * t49;
	t48 = (t33 * t47 + t36 * t46) / t30 ^ 2;
	t28 = 0.1e1 / t30;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (153->13), mult. (451->28), div. (25->4), fcn. (552->7), ass. (0->20)
	t63 = sin(pkin(11));
	t64 = cos(pkin(11));
	t65 = cos(pkin(10));
	t66 = sin(qJ(2));
	t67 = cos(qJ(2));
	t78 = cos(pkin(6));
	t74 = t67 * t78;
	t75 = t66 * t78;
	t77 = sin(pkin(10));
	t54 = (t63 * t75 - t64 * t74) * t77 + t65 * (-t67 * t63 - t66 * t64);
	t48 = t54 * qJD(2);
	t72 = (-t63 * t74 - t64 * t75) * t77 - t65 * (t66 * t63 - t67 * t64);
	t51 = 0.1e1 / t72 ^ 2;
	t84 = t51 * t54 ^ 2;
	t50 = 0.1e1 / t72;
	t83 = t50 * t84;
	t82 = t51 * t72;
	t76 = t82 * t48;
	t46 = 0.1e1 + t84;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0.2e1 * (-t50 * t72 - t84) / t46 ^ 2 * (-t48 * t83 - t76) + (-0.2e1 * t76 + (t50 - t82 - 0.2e1 * t83) * t48) / t46, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:41
	% DurationCPUTime: 0.73s
	% Computational Cost: add. (1747->58), mult. (5333->135), div. (281->12), fcn. (6885->13), ass. (0->71)
	t185 = sin(qJ(4));
	t187 = cos(qJ(4));
	t184 = cos(pkin(6));
	t179 = sin(pkin(11));
	t182 = cos(pkin(11));
	t186 = sin(qJ(2));
	t188 = cos(qJ(2));
	t200 = t188 * t179 + t186 * t182;
	t171 = t200 * t184;
	t174 = t186 * t179 - t188 * t182;
	t180 = sin(pkin(10));
	t183 = cos(pkin(10));
	t201 = -t180 * t171 - t183 * t174;
	t181 = sin(pkin(6));
	t205 = t180 * t181;
	t197 = -t185 * t201 + t187 * t205;
	t214 = t197 * qJD(4);
	t172 = t174 * qJD(2);
	t196 = t174 * t184;
	t156 = -t180 * t200 - t183 * t196;
	t169 = t174 * t181;
	t142 = atan2(t156, t169);
	t137 = sin(t142);
	t138 = cos(t142);
	t131 = t137 * t156 + t138 * t169;
	t128 = 0.1e1 / t131;
	t148 = t185 * t205 + t187 * t201;
	t144 = 0.1e1 / t148;
	t166 = 0.1e1 / t169;
	t129 = 0.1e1 / t131 ^ 2;
	t145 = 0.1e1 / t148 ^ 2;
	t167 = 0.1e1 / t169 ^ 2;
	t153 = t156 ^ 2;
	t141 = t153 * t167 + 0.1e1;
	t139 = 0.1e1 / t141;
	t195 = qJD(2) * t171;
	t149 = t180 * t172 - t183 * t195;
	t170 = t200 * t181;
	t163 = qJD(2) * t170;
	t208 = t156 * t167;
	t122 = (t149 * t166 - t163 * t208) * t139;
	t202 = -t137 * t169 + t138 * t156;
	t119 = t202 * t122 + t137 * t149 + t138 * t163;
	t213 = t119 * t128 * t129;
	t143 = t197 ^ 2;
	t134 = t143 * t145 + 0.1e1;
	t165 = t184 * t172;
	t173 = t200 * qJD(2);
	t152 = t180 * t165 - t183 * t173;
	t135 = t148 * qJD(4) + t152 * t185;
	t209 = t145 * t197;
	t136 = t152 * t187 + t214;
	t210 = t136 * t144 * t145;
	t212 = (-t135 * t209 - t143 * t210) / t134 ^ 2;
	t158 = t180 * t196 - t183 * t200;
	t211 = t129 * t158;
	t207 = t156 * t170;
	t206 = t163 * t166 * t167;
	t199 = -t144 * t185 - t187 * t209;
	t155 = -t183 * t171 + t180 * t174;
	t198 = -t155 * t166 + t167 * t207;
	t164 = t181 * t172;
	t154 = t158 ^ 2;
	t151 = t183 * t172 + t180 * t195;
	t150 = t183 * t165 + t180 * t173;
	t132 = 0.1e1 / t134;
	t126 = t154 * t129 + 0.1e1;
	t123 = t198 * t139;
	t120 = -t202 * t123 + t137 * t155 + t138 * t170;
	t118 = 0.2e1 * t198 / t141 ^ 2 * (t149 * t208 - t153 * t206) + (0.2e1 * t206 * t207 + t150 * t166 + (-t149 * t170 - t155 * t163 + t156 * t164) * t167) * t139;
	t1 = [0, t118, 0, 0, 0, 0; 0, 0.2e1 * (-t120 * t211 - t128 * t201) / t126 ^ 2 * (t151 * t211 - t154 * t213) + (t152 * t128 + (-t119 * t201 + t120 * t151) * t129 + (-0.2e1 * t120 * t213 + ((t118 * t156 - t123 * t149 - t164 + (t123 * t169 + t155) * t122) * t138 + (-t118 * t169 + t123 * t163 + t150 + (t123 * t156 - t170) * t122) * t137) * t129) * t158) / t126, 0, 0, 0, 0; 0, 0.2e1 * t199 * t158 * t212 + (-t199 * t151 + ((qJD(4) * t144 - 0.2e1 * t197 * t210) * t187 + (-t135 * t187 + (-t136 - t214) * t185) * t145) * t158) * t132, 0, -0.2e1 * t212 - 0.2e1 * (t132 * t135 * t145 - (-t132 * t210 - t145 * t212) * t197) * t197, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:41
	% EndTime: 2019-10-09 21:33:42
	% DurationCPUTime: 1.54s
	% Computational Cost: add. (4524->94), mult. (13242->201), div. (524->12), fcn. (17302->13), ass. (0->92)
	t213 = sin(pkin(11));
	t216 = cos(pkin(11));
	t220 = sin(qJ(2));
	t222 = cos(qJ(2));
	t206 = t220 * t213 - t222 * t216;
	t218 = cos(pkin(6));
	t230 = t206 * t218;
	t200 = qJD(2) * t230;
	t234 = t222 * t213 + t220 * t216;
	t205 = t234 * qJD(2);
	t214 = sin(pkin(10));
	t217 = cos(pkin(10));
	t179 = -t217 * t200 - t214 * t205;
	t203 = t234 * t218;
	t186 = t217 * t203 - t214 * t206;
	t219 = sin(qJ(4));
	t215 = sin(pkin(6));
	t246 = t215 * t219;
	t238 = t217 * t246;
	t221 = cos(qJ(4));
	t241 = qJD(4) * t221;
	t155 = -qJD(4) * t238 + t179 * t219 + t186 * t241;
	t245 = t215 * t221;
	t173 = t186 * t219 + t217 * t245;
	t170 = t173 ^ 2;
	t202 = t234 * t215;
	t194 = t202 * t219 - t218 * t221;
	t192 = 0.1e1 / t194 ^ 2;
	t166 = t170 * t192 + 0.1e1;
	t164 = 0.1e1 / t166;
	t195 = t202 * t221 + t218 * t219;
	t201 = t206 * t215;
	t199 = qJD(2) * t201;
	t168 = t195 * qJD(4) - t199 * t219;
	t191 = 0.1e1 / t194;
	t249 = t173 * t192;
	t143 = (-t155 * t191 + t168 * t249) * t164;
	t167 = atan2(-t173, t194);
	t162 = sin(t167);
	t163 = cos(t167);
	t236 = -t162 * t194 - t163 * t173;
	t140 = t236 * t143 - t162 * t155 + t163 * t168;
	t154 = -t162 * t173 + t163 * t194;
	t151 = 0.1e1 / t154;
	t152 = 0.1e1 / t154 ^ 2;
	t259 = t140 * t151 * t152;
	t181 = t214 * t200 - t217 * t205;
	t235 = -t214 * t203 - t217 * t206;
	t231 = t214 * t245 - t219 * t235;
	t158 = t231 * qJD(4) + t181 * t221;
	t177 = t214 * t246 + t221 * t235;
	t172 = t177 ^ 2;
	t188 = t214 * t230 - t217 * t234;
	t183 = 0.1e1 / t188 ^ 2;
	t161 = t172 * t183 + 0.1e1;
	t204 = t206 * qJD(2);
	t229 = t218 * t205;
	t180 = t217 * t204 + t214 * t229;
	t182 = 0.1e1 / t188;
	t184 = t182 * t183;
	t258 = 0.2e1 * (t177 * t183 * t158 - t172 * t184 * t180) / t161 ^ 2;
	t257 = -0.2e1 * t231 * t259;
	t185 = -t214 * t234 - t217 * t230;
	t232 = -t185 * t191 - t201 * t249;
	t256 = t219 * t232;
	t250 = t168 * t191 * t192;
	t255 = -0.2e1 * (t155 * t249 - t170 * t250) / t166 ^ 2;
	t253 = t152 * t231;
	t252 = t162 * t231;
	t251 = t163 * t231;
	t248 = t177 * t235;
	t247 = t188 * t219;
	t171 = t231 ^ 2;
	t150 = t171 * t152 + 0.1e1;
	t157 = t177 * qJD(4) + t181 * t219;
	t240 = 0.2e1 * (-t157 * t253 - t171 * t259) / t150 ^ 2;
	t237 = -0.2e1 * t173 * t250;
	t175 = t186 * t221 - t238;
	t233 = -t175 * t191 + t195 * t249;
	t198 = t215 * t205;
	t178 = t214 * t204 - t217 * t229;
	t169 = -t194 * qJD(4) - t199 * t221;
	t159 = 0.1e1 / t161;
	t156 = -t173 * qJD(4) + t179 * t221;
	t148 = 0.1e1 / t150;
	t145 = t164 * t256;
	t144 = t233 * t164;
	t142 = (-t162 * t185 - t163 * t201) * t219 + t236 * t145;
	t141 = t236 * t144 - t162 * t175 + t163 * t195;
	t139 = t233 * t255 + (t195 * t237 - t156 * t191 + (t155 * t195 + t168 * t175 + t169 * t173) * t192) * t164;
	t137 = t255 * t256 + (t232 * t241 + (-t201 * t237 - t178 * t191 + (-t155 * t201 + t168 * t185 - t173 * t198) * t192) * t219) * t164;
	t1 = [0, t137, 0, t139, 0, 0; 0, (-t142 * t253 - t151 * t247) * t240 + ((t180 * t219 + t188 * t241) * t151 + t142 * t257 + (-t142 * t157 - t247 * t140 + (-t201 * t241 - t137 * t173 - t145 * t155 - t198 * t219 + (-t145 * t194 - t185 * t219) * t143) * t251 + (-t185 * t241 - t137 * t194 - t145 * t168 - t178 * t219 + (t145 * t173 + t201 * t219) * t143) * t252) * t152) * t148, 0, (-t141 * t253 - t151 * t177) * t240 + (t141 * t257 + t158 * t151 + (-t177 * t140 - t141 * t157 + (-t139 * t173 - t144 * t155 + t169 + (-t144 * t194 - t175) * t143) * t251 + (-t139 * t194 - t144 * t168 - t156 + (t144 * t173 - t195) * t143) * t252) * t152) * t148, 0, 0; 0, (t182 * t188 * t221 + t183 * t248) * t258 + (qJD(4) * t182 * t247 + (-t158 * t235 - t181 * t177) * t183 + (0.2e1 * t184 * t248 + (t183 * t188 - t182) * t221) * t180) * t159, 0, t231 * t182 * t258 + (t180 * t183 * t231 + t157 * t182) * t159, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:41
	% EndTime: 2019-10-09 21:33:43
	% DurationCPUTime: 1.80s
	% Computational Cost: add. (5642->115), mult. (16325->236), div. (559->12), fcn. (21250->15), ass. (0->112)
	t278 = sin(qJ(4));
	t281 = cos(qJ(4));
	t276 = cos(pkin(6));
	t273 = sin(pkin(11));
	t275 = cos(pkin(11));
	t279 = sin(qJ(2));
	t282 = cos(qJ(2));
	t302 = t282 * t273 + t279 * t275;
	t263 = t302 * t276;
	t266 = t279 * t273 - t282 * t275;
	t274 = sin(pkin(10));
	t334 = cos(pkin(10));
	t297 = t263 * t334 - t274 * t266;
	t333 = sin(pkin(6));
	t305 = t334 * t333;
	t236 = t297 * t278 + t281 * t305;
	t294 = t266 * t276;
	t260 = qJD(2) * t294;
	t265 = t302 * qJD(2);
	t242 = -t260 * t334 - t274 * t265;
	t214 = -t236 * qJD(4) + t242 * t281;
	t237 = -t278 * t305 + t281 * t297;
	t234 = t237 ^ 2;
	t309 = t282 * t333;
	t310 = t279 * t333;
	t289 = -t273 * t309 - t275 * t310;
	t255 = t276 * t278 - t281 * t289;
	t252 = 0.1e1 / t255 ^ 2;
	t229 = t234 * t252 + 0.1e1;
	t227 = 0.1e1 / t229;
	t254 = t276 * t281 + t278 * t289;
	t261 = -t273 * t310 + t275 * t309;
	t259 = t261 * qJD(2);
	t232 = qJD(4) * t254 + t259 * t281;
	t251 = 0.1e1 / t255;
	t322 = t237 * t252;
	t197 = (-t214 * t251 + t232 * t322) * t227;
	t230 = atan2(-t237, t255);
	t225 = sin(t230);
	t226 = cos(t230);
	t304 = -t225 * t255 - t226 * t237;
	t193 = t197 * t304 - t225 * t214 + t226 * t232;
	t209 = -t225 * t237 + t226 * t255;
	t206 = 0.1e1 / t209;
	t207 = 0.1e1 / t209 ^ 2;
	t339 = t193 * t206 * t207;
	t249 = t274 * t294 - t302 * t334;
	t277 = sin(qJ(6));
	t280 = cos(qJ(6));
	t295 = -t274 * t263 - t266 * t334;
	t311 = t274 * t333;
	t293 = -t278 * t295 + t281 * t311;
	t303 = t249 * t277 - t280 * t293;
	t338 = qJD(6) * t303;
	t240 = t278 * t311 + t281 * t295;
	t337 = 0.2e1 * t240 * t339;
	t247 = -t274 * t302 - t294 * t334;
	t298 = -t247 * t251 + t261 * t322;
	t336 = t281 * t298;
	t323 = t232 * t251 * t252;
	t335 = -0.2e1 * (t214 * t322 - t234 * t323) / t229 ^ 2;
	t320 = t249 * t280;
	t222 = -t277 * t293 - t320;
	t218 = 0.1e1 / t222;
	t219 = 0.1e1 / t222 ^ 2;
	t296 = t274 * t260 - t265 * t334;
	t215 = qJD(4) * t240 + t278 * t296;
	t264 = t266 * qJD(2);
	t292 = t276 * t265;
	t243 = t264 * t334 + t274 * t292;
	t204 = qJD(6) * t222 - t215 * t280 - t243 * t277;
	t217 = t303 ^ 2;
	t212 = t217 * t219 + 0.1e1;
	t327 = t219 * t303;
	t205 = t215 * t277 - t243 * t280 + t338;
	t331 = t205 * t218 * t219;
	t332 = (-t204 * t327 - t217 * t331) / t212 ^ 2;
	t330 = t207 * t240;
	t216 = qJD(4) * t293 + t281 * t296;
	t329 = t216 * t207;
	t328 = t218 * t280;
	t326 = t303 * t277;
	t325 = t225 * t240;
	t324 = t226 * t240;
	t321 = t249 * t278;
	t319 = t249 * t281;
	t315 = qJD(4) * t278;
	t235 = t240 ^ 2;
	t203 = t235 * t207 + 0.1e1;
	t314 = 0.2e1 * (-t235 * t339 + t240 * t329) / t203 ^ 2;
	t313 = 0.2e1 * t332;
	t308 = -0.2e1 * t303 * t331;
	t307 = -0.2e1 * t237 * t323;
	t306 = qJD(6) * t321 + t296;
	t300 = -t219 * t326 + t328;
	t299 = t236 * t251 + t254 * t322;
	t290 = qJD(4) * t319 - qJD(6) * t295 + t243 * t278;
	t258 = t289 * qJD(2);
	t241 = t274 * t264 - t292 * t334;
	t231 = -qJD(4) * t255 - t259 * t278;
	t224 = t277 * t321 + t280 * t295;
	t223 = t277 * t295 - t278 * t320;
	t213 = qJD(4) * t237 + t242 * t278;
	t210 = 0.1e1 / t212;
	t201 = 0.1e1 / t203;
	t199 = t227 * t336;
	t198 = t299 * t227;
	t195 = (-t225 * t247 + t226 * t261) * t281 + t304 * t199;
	t194 = t198 * t304 + t225 * t236 + t226 * t254;
	t192 = t299 * t335 + (t254 * t307 + t213 * t251 + (t214 * t254 + t231 * t237 - t232 * t236) * t252) * t227;
	t190 = t335 * t336 + (-t298 * t315 + (t261 * t307 - t241 * t251 + (t214 * t261 + t232 * t247 + t237 * t258) * t252) * t281) * t227;
	t1 = [0, t190, 0, t192, 0, 0; 0, (t195 * t330 - t206 * t319) * t314 + ((t243 * t281 - t249 * t315) * t206 + (-t329 + t337) * t195 + (-t319 * t193 - (-t261 * t315 - t190 * t237 - t199 * t214 + t258 * t281 + (-t199 * t255 - t247 * t281) * t197) * t324 - (t247 * t315 - t190 * t255 - t199 * t232 - t241 * t281 + (t199 * t237 - t261 * t281) * t197) * t325) * t207) * t201, 0, (t194 * t330 - t206 * t293) * t314 + (t194 * t337 - t215 * t206 + (-t293 * t193 - t194 * t216 - (-t192 * t237 - t198 * t214 + t231 + (-t198 * t255 + t236) * t197) * t324 - (-t192 * t255 - t198 * t232 + t213 + (t198 * t237 - t254) * t197) * t325) * t207) * t201, 0, 0; 0, (-t218 * t223 - t224 * t327) * t313 + (t224 * t308 + t306 * t218 * t277 - t290 * t328 + (t280 * t303 * t306 - t224 * t204 - t223 * t205 + t290 * t326) * t219) * t210, 0, t300 * t240 * t313 + (-t300 * t216 + ((qJD(6) * t218 + t308) * t277 + (-t204 * t277 + (t205 + t338) * t280) * t219) * t240) * t210, 0, -0.2e1 * t332 - 0.2e1 * (t204 * t219 * t210 - (-t210 * t331 - t219 * t332) * t303) * t303;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end