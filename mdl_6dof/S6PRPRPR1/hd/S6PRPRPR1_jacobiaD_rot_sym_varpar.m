% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR1
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
%   Wie in S6PRPRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
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
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:30:00
	% DurationCPUTime: 0.71s
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
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:30:00
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (1932->59), mult. (5333->135), div. (281->12), fcn. (6885->13), ass. (0->72)
	t190 = qJ(4) + pkin(12);
	t188 = sin(t190);
	t189 = cos(t190);
	t196 = cos(pkin(6));
	t191 = sin(pkin(11));
	t194 = cos(pkin(11));
	t197 = sin(qJ(2));
	t198 = cos(qJ(2));
	t210 = t191 * t198 + t197 * t194;
	t180 = t210 * t196;
	t183 = t197 * t191 - t194 * t198;
	t192 = sin(pkin(10));
	t195 = cos(pkin(10));
	t211 = -t180 * t192 - t183 * t195;
	t193 = sin(pkin(6));
	t214 = t192 * t193;
	t207 = -t188 * t211 + t189 * t214;
	t224 = t207 * qJD(4);
	t181 = t183 * qJD(2);
	t206 = t183 * t196;
	t165 = -t192 * t210 - t195 * t206;
	t178 = t183 * t193;
	t151 = atan2(t165, t178);
	t146 = sin(t151);
	t147 = cos(t151);
	t140 = t146 * t165 + t147 * t178;
	t137 = 0.1e1 / t140;
	t157 = t188 * t214 + t189 * t211;
	t153 = 0.1e1 / t157;
	t175 = 0.1e1 / t178;
	t138 = 0.1e1 / t140 ^ 2;
	t154 = 0.1e1 / t157 ^ 2;
	t176 = 0.1e1 / t178 ^ 2;
	t162 = t165 ^ 2;
	t150 = t162 * t176 + 0.1e1;
	t148 = 0.1e1 / t150;
	t205 = qJD(2) * t180;
	t158 = t192 * t181 - t195 * t205;
	t179 = t210 * t193;
	t172 = qJD(2) * t179;
	t218 = t165 * t176;
	t131 = (t158 * t175 - t172 * t218) * t148;
	t212 = -t146 * t178 + t147 * t165;
	t128 = t212 * t131 + t146 * t158 + t147 * t172;
	t223 = t128 * t137 * t138;
	t152 = t207 ^ 2;
	t143 = t152 * t154 + 0.1e1;
	t174 = t196 * t181;
	t182 = t210 * qJD(2);
	t161 = t174 * t192 - t182 * t195;
	t144 = t157 * qJD(4) + t161 * t188;
	t219 = t154 * t207;
	t145 = t161 * t189 + t224;
	t220 = t145 * t153 * t154;
	t222 = (-t144 * t219 - t152 * t220) / t143 ^ 2;
	t167 = t192 * t206 - t195 * t210;
	t221 = t138 * t167;
	t217 = t165 * t179;
	t216 = t172 * t175 * t176;
	t209 = -t153 * t188 - t189 * t219;
	t164 = -t180 * t195 + t183 * t192;
	t208 = -t164 * t175 + t176 * t217;
	t173 = t193 * t181;
	t163 = t167 ^ 2;
	t160 = t181 * t195 + t192 * t205;
	t159 = t174 * t195 + t182 * t192;
	t141 = 0.1e1 / t143;
	t135 = t138 * t163 + 0.1e1;
	t132 = t208 * t148;
	t129 = -t212 * t132 + t146 * t164 + t147 * t179;
	t127 = 0.2e1 * t208 / t150 ^ 2 * (t158 * t218 - t162 * t216) + (0.2e1 * t216 * t217 + t159 * t175 + (-t158 * t179 - t164 * t172 + t165 * t173) * t176) * t148;
	t1 = [0, t127, 0, 0, 0, 0; 0, 0.2e1 * (-t129 * t221 - t137 * t211) / t135 ^ 2 * (t160 * t221 - t163 * t223) + (t161 * t137 + (-t128 * t211 + t129 * t160) * t138 + (-0.2e1 * t129 * t223 + ((t127 * t165 - t132 * t158 - t173 + (t132 * t178 + t164) * t131) * t147 + (-t127 * t178 + t132 * t172 + t159 + (t132 * t165 - t179) * t131) * t146) * t138) * t167) / t135, 0, 0, 0, 0; 0, 0.2e1 * t209 * t167 * t222 + (-t209 * t160 + ((qJD(4) * t153 - 0.2e1 * t207 * t220) * t189 + (-t144 * t189 + (-t145 - t224) * t188) * t154) * t167) * t141, 0, -0.2e1 * t222 - 0.2e1 * (t141 * t144 * t154 - (-t141 * t220 - t154 * t222) * t207) * t207, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:30:01
	% DurationCPUTime: 1.94s
	% Computational Cost: add. (8444->115), mult. (16325->235), div. (559->12), fcn. (21250->15), ass. (0->112)
	t284 = sin(pkin(11));
	t287 = cos(pkin(11));
	t291 = sin(qJ(2));
	t293 = cos(qJ(2));
	t274 = t291 * t284 - t293 * t287;
	t289 = cos(pkin(6));
	t302 = t274 * t289;
	t268 = qJD(2) * t302;
	t307 = t293 * t284 + t291 * t287;
	t273 = t307 * qJD(2);
	t285 = sin(pkin(10));
	t288 = cos(pkin(10));
	t249 = -t288 * t268 - t285 * t273;
	t271 = t307 * t289;
	t255 = t288 * t271 - t285 * t274;
	t283 = qJ(4) + pkin(12);
	t281 = sin(t283);
	t286 = sin(pkin(6));
	t323 = t286 * t288;
	t313 = t281 * t323;
	t282 = cos(t283);
	t319 = qJD(4) * t282;
	t221 = -qJD(4) * t313 + t249 * t281 + t255 * t319;
	t243 = t255 * t281 + t282 * t323;
	t241 = t243 ^ 2;
	t270 = t307 * t286;
	t262 = t270 * t281 - t289 * t282;
	t260 = 0.1e1 / t262 ^ 2;
	t235 = t241 * t260 + 0.1e1;
	t233 = 0.1e1 / t235;
	t263 = t270 * t282 + t289 * t281;
	t269 = t274 * t286;
	t267 = qJD(2) * t269;
	t239 = t263 * qJD(4) - t267 * t281;
	t259 = 0.1e1 / t262;
	t328 = t243 * t260;
	t205 = (-t221 * t259 + t239 * t328) * t233;
	t236 = atan2(-t243, t262);
	t231 = sin(t236);
	t232 = cos(t236);
	t310 = -t231 * t262 - t232 * t243;
	t201 = t310 * t205 - t231 * t221 + t232 * t239;
	t215 = -t231 * t243 + t232 * t262;
	t212 = 0.1e1 / t215;
	t213 = 0.1e1 / t215 ^ 2;
	t342 = t201 * t212 * t213;
	t308 = -t285 * t271 - t288 * t274;
	t324 = t285 * t286;
	t303 = -t281 * t308 + t282 * t324;
	t341 = -0.2e1 * t303 * t342;
	t254 = -t285 * t307 - t288 * t302;
	t304 = -t254 * t259 - t269 * t328;
	t340 = t281 * t304;
	t329 = t239 * t259 * t260;
	t339 = -0.2e1 * (t221 * t328 - t241 * t329) / t235 ^ 2;
	t247 = t281 * t324 + t282 * t308;
	t292 = cos(qJ(6));
	t257 = t285 * t302 - t288 * t307;
	t290 = sin(qJ(6));
	t326 = t257 * t290;
	t230 = t247 * t292 - t326;
	t226 = 0.1e1 / t230;
	t227 = 0.1e1 / t230 ^ 2;
	t309 = t285 * t268 - t288 * t273;
	t224 = t303 * qJD(4) + t282 * t309;
	t272 = t274 * qJD(2);
	t301 = t289 * t273;
	t250 = t288 * t272 + t285 * t301;
	t216 = t230 * qJD(6) + t224 * t290 + t250 * t292;
	t325 = t257 * t292;
	t229 = t247 * t290 + t325;
	t225 = t229 ^ 2;
	t220 = t225 * t227 + 0.1e1;
	t333 = t227 * t229;
	t318 = qJD(6) * t229;
	t217 = t224 * t292 - t250 * t290 - t318;
	t336 = t217 * t226 * t227;
	t338 = (t216 * t333 - t225 * t336) / t220 ^ 2;
	t337 = t213 * t303;
	t223 = t247 * qJD(4) + t281 * t309;
	t335 = t223 * t213;
	t334 = t226 * t290;
	t332 = t229 * t292;
	t331 = t231 * t303;
	t330 = t232 * t303;
	t327 = t257 * t281;
	t242 = t303 ^ 2;
	t211 = t242 * t213 + 0.1e1;
	t317 = 0.2e1 * (-t242 * t342 - t303 * t335) / t211 ^ 2;
	t316 = -0.2e1 * t338;
	t314 = t229 * t336;
	t312 = -0.2e1 * t243 * t329;
	t311 = qJD(6) * t257 * t282 - t309;
	t306 = t227 * t332 - t334;
	t245 = t255 * t282 - t313;
	t305 = -t245 * t259 + t263 * t328;
	t300 = -qJD(4) * t327 + qJD(6) * t308 + t250 * t282;
	t266 = t286 * t273;
	t248 = t285 * t272 - t288 * t301;
	t240 = -t262 * qJD(4) - t267 * t282;
	t238 = t282 * t325 + t290 * t308;
	t237 = t282 * t326 - t292 * t308;
	t222 = -t243 * qJD(4) + t249 * t282;
	t218 = 0.1e1 / t220;
	t209 = 0.1e1 / t211;
	t207 = t233 * t340;
	t206 = t305 * t233;
	t203 = (-t231 * t254 - t232 * t269) * t281 + t310 * t207;
	t202 = t310 * t206 - t231 * t245 + t232 * t263;
	t199 = t305 * t339 + (t263 * t312 - t222 * t259 + (t221 * t263 + t239 * t245 + t240 * t243) * t260) * t233;
	t198 = t339 * t340 + (t304 * t319 + (-t269 * t312 - t248 * t259 + (-t221 * t269 + t239 * t254 - t243 * t266) * t260) * t281) * t233;
	t1 = [0, t198, 0, t199, 0, 0; 0, (-t203 * t337 - t212 * t327) * t317 + ((t250 * t281 + t257 * t319) * t212 + (-t335 + t341) * t203 + (-t327 * t201 + (-t269 * t319 - t198 * t243 - t207 * t221 - t266 * t281 + (-t207 * t262 - t254 * t281) * t205) * t330 + (-t254 * t319 - t198 * t262 - t207 * t239 - t248 * t281 + (t207 * t243 + t269 * t281) * t205) * t331) * t213) * t209, 0, (-t202 * t337 - t212 * t247) * t317 + (t202 * t341 + t224 * t212 + (-t247 * t201 - t202 * t223 + (-t199 * t243 - t206 * t221 + t240 + (-t206 * t262 - t245) * t205) * t330 + (-t199 * t262 - t206 * t239 - t222 + (t206 * t243 - t263) * t205) * t331) * t213) * t209, 0, 0; 0, 0.2e1 * (-t226 * t237 + t238 * t333) * t338 + (0.2e1 * t238 * t314 + t311 * t226 * t292 + t300 * t334 + (t311 * t229 * t290 - t238 * t216 - t237 * t217 - t300 * t332) * t227) * t218, 0, -t306 * t303 * t316 + (t306 * t223 - ((-qJD(6) * t226 - 0.2e1 * t314) * t292 + (t216 * t292 + (t217 - t318) * t290) * t227) * t303) * t218, 0, t316 + 0.2e1 * (t216 * t227 * t218 + (-t218 * t336 - t227 * t338) * t229) * t229;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end