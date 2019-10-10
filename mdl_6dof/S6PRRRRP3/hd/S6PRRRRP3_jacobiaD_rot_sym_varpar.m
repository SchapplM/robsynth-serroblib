% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRP3
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
%   Wie in S6PRRRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:05
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:05:49
	% EndTime: 2019-10-09 23:05:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:05:49
	% EndTime: 2019-10-09 23:05:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:05:49
	% EndTime: 2019-10-09 23:05:49
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(11));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(11)) * cos(pkin(6));
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
	% StartTime: 2019-10-09 23:05:49
	% EndTime: 2019-10-09 23:05:50
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (756->55), mult. (2271->133), div. (423->14), fcn. (2956->11), ass. (0->65)
	t139 = sin(qJ(2));
	t141 = cos(qJ(2));
	t136 = sin(pkin(11));
	t166 = cos(pkin(6));
	t155 = t136 * t166;
	t165 = cos(pkin(11));
	t127 = -t139 * t155 + t165 * t141;
	t138 = sin(qJ(3));
	t140 = cos(qJ(3));
	t137 = sin(pkin(6));
	t160 = t136 * t137;
	t149 = -t127 * t138 + t140 * t160;
	t170 = t149 * qJD(3);
	t152 = t166 * t165;
	t123 = t136 * t139 - t141 * t152;
	t159 = t137 * t141;
	t113 = atan2(-t123, -t159);
	t111 = sin(t113);
	t112 = cos(t113);
	t98 = -t111 * t123 - t112 * t159;
	t95 = 0.1e1 / t98;
	t110 = t127 * t140 + t138 * t160;
	t106 = 0.1e1 / t110;
	t133 = 0.1e1 / t141;
	t107 = 0.1e1 / t110 ^ 2;
	t134 = 0.1e1 / t141 ^ 2;
	t96 = 0.1e1 / t98 ^ 2;
	t105 = t149 ^ 2;
	t102 = t105 * t107 + 0.1e1;
	t148 = -t165 * t139 - t141 * t155;
	t119 = t148 * qJD(2);
	t103 = t110 * qJD(3) + t119 * t138;
	t163 = t107 * t149;
	t104 = t119 * t140 + t170;
	t164 = t104 * t106 * t107;
	t169 = 0.1e1 / t102 ^ 2 * (-t103 * t163 - t105 * t164);
	t125 = t136 * t141 + t139 * t152;
	t161 = t134 * t139;
	t156 = t123 * t161;
	t150 = t125 * t133 + t156;
	t121 = t123 ^ 2;
	t132 = 0.1e1 / t137 ^ 2;
	t116 = t121 * t132 * t134 + 0.1e1;
	t114 = 0.1e1 / t116;
	t131 = 0.1e1 / t137;
	t162 = t114 * t131;
	t91 = t150 * t162;
	t168 = t123 * t91;
	t167 = t148 * t96;
	t158 = qJD(2) * t139;
	t157 = -0.2e1 * t169;
	t151 = -t106 * t138 - t140 * t163;
	t135 = t133 * t134;
	t122 = t148 ^ 2;
	t120 = t127 * qJD(2);
	t118 = t125 * qJD(2);
	t117 = qJD(2) * t123;
	t100 = 0.1e1 / t102;
	t97 = t95 * t96;
	t94 = t122 * t96 + 0.1e1;
	t90 = (qJD(2) * t156 + t118 * t133) * t162;
	t88 = (t137 * t139 - t168) * t112 + (t91 * t159 - t125) * t111;
	t87 = (-t123 * t90 + t137 * t158) * t112 + (t90 * t159 - t118) * t111;
	t86 = (-0.2e1 * t150 * (t118 * t123 * t134 + t121 * t135 * t158) * t132 / t116 ^ 2 + (t118 * t161 - t117 * t133 + (t125 * t161 + (0.2e1 * t135 * t139 ^ 2 + t133) * t123) * qJD(2)) * t114) * t131;
	t1 = [0, t86, 0, 0, 0, 0; 0, 0.2e1 * (-t127 * t95 - t88 * t167) / t94 ^ 2 * (-t122 * t97 * t87 - t120 * t167) + (-t88 * t120 * t96 + t119 * t95 + (-0.2e1 * t148 * t88 * t97 - t127 * t96) * t87 - (-(-t118 * t91 - t123 * t86 - t125 * t90 + (t90 * t91 + qJD(2)) * t159) * t112 - (t90 * t168 + t117 + (t141 * t86 + (-qJD(2) * t91 - t90) * t139) * t137) * t111) * t167) / t94, 0, 0, 0, 0; 0, -t151 * t148 * t157 + (t151 * t120 - ((-qJD(3) * t106 + 0.2e1 * t149 * t164) * t140 + (t103 * t140 + (t104 + t170) * t138) * t107) * t148) * t100, t157 - 0.2e1 * (t100 * t103 * t107 - (-t100 * t164 - t107 * t169) * t149) * t149, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:05:49
	% EndTime: 2019-10-09 23:05:51
	% DurationCPUTime: 1.38s
	% Computational Cost: add. (3002->109), mult. (9085->226), div. (559->12), fcn. (11668->13), ass. (0->103)
	t207 = sin(pkin(11));
	t209 = cos(pkin(11));
	t213 = sin(qJ(2));
	t210 = cos(pkin(6));
	t216 = cos(qJ(2));
	t242 = t210 * t216;
	t197 = -t207 * t213 + t209 * t242;
	t190 = t197 * qJD(2);
	t243 = t210 * t213;
	t198 = t207 * t216 + t209 * t243;
	t212 = sin(qJ(3));
	t208 = sin(pkin(6));
	t246 = t208 * t212;
	t232 = t209 * t246;
	t215 = cos(qJ(3));
	t239 = qJD(3) * t215;
	t162 = -qJD(3) * t232 + t190 * t212 + t198 * t239;
	t245 = t208 * t215;
	t182 = t198 * t212 + t209 * t245;
	t180 = t182 ^ 2;
	t201 = -t210 * t215 + t213 * t246;
	t195 = 0.1e1 / t201 ^ 2;
	t176 = t180 * t195 + 0.1e1;
	t174 = 0.1e1 / t176;
	t202 = t210 * t212 + t213 * t245;
	t240 = qJD(2) * t216;
	t231 = t208 * t240;
	t187 = t202 * qJD(3) + t212 * t231;
	t194 = 0.1e1 / t201;
	t250 = t182 * t195;
	t146 = (-t162 * t194 + t187 * t250) * t174;
	t177 = atan2(-t182, t201);
	t172 = sin(t177);
	t173 = cos(t177);
	t228 = -t172 * t201 - t173 * t182;
	t142 = t228 * t146 - t172 * t162 + t173 * t187;
	t156 = -t172 * t182 + t173 * t201;
	t153 = 0.1e1 / t156;
	t154 = 0.1e1 / t156 ^ 2;
	t263 = t142 * t153 * t154;
	t233 = t207 * t243;
	t200 = t209 * t216 - t233;
	t225 = -t200 * t212 + t207 * t245;
	t262 = -0.2e1 * t225 * t263;
	t244 = t208 * t216;
	t224 = -t194 * t197 + t244 * t250;
	t261 = t212 * t224;
	t249 = t187 * t194 * t195;
	t260 = -0.2e1 * (t162 * t250 - t180 * t249) / t176 ^ 2;
	t186 = t200 * t215 + t207 * t246;
	t199 = t207 * t242 + t209 * t213;
	t211 = sin(qJ(4));
	t214 = cos(qJ(4));
	t171 = t186 * t214 + t199 * t211;
	t167 = 0.1e1 / t171;
	t168 = 0.1e1 / t171 ^ 2;
	t192 = t199 * qJD(2);
	t165 = t225 * qJD(3) - t192 * t215;
	t193 = -qJD(2) * t233 + t209 * t240;
	t157 = t171 * qJD(4) + t165 * t211 - t193 * t214;
	t170 = t186 * t211 - t199 * t214;
	t166 = t170 ^ 2;
	t161 = t166 * t168 + 0.1e1;
	t254 = t168 * t170;
	t238 = qJD(4) * t170;
	t158 = t165 * t214 + t193 * t211 - t238;
	t257 = t158 * t167 * t168;
	t259 = (t157 * t254 - t166 * t257) / t161 ^ 2;
	t258 = t154 * t225;
	t164 = t186 * qJD(3) - t192 * t212;
	t256 = t164 * t154;
	t255 = t167 * t211;
	t253 = t170 * t214;
	t252 = t172 * t225;
	t251 = t173 * t225;
	t248 = t199 * t212;
	t247 = t199 * t215;
	t241 = qJD(2) * t213;
	t181 = t225 ^ 2;
	t152 = t181 * t154 + 0.1e1;
	t237 = 0.2e1 * (-t181 * t263 - t225 * t256) / t152 ^ 2;
	t236 = -0.2e1 * t259;
	t234 = t170 * t257;
	t230 = -0.2e1 * t182 * t249;
	t229 = qJD(4) * t247 - t192;
	t227 = t168 * t253 - t255;
	t184 = t198 * t215 - t232;
	t226 = -t184 * t194 + t202 * t250;
	t223 = qJD(3) * t248 + qJD(4) * t200 - t193 * t215;
	t191 = t198 * qJD(2);
	t188 = -t201 * qJD(3) + t215 * t231;
	t179 = t200 * t211 - t214 * t247;
	t178 = -t200 * t214 - t211 * t247;
	t163 = -t182 * qJD(3) + t190 * t215;
	t159 = 0.1e1 / t161;
	t149 = 0.1e1 / t152;
	t148 = t174 * t261;
	t147 = t226 * t174;
	t144 = (-t172 * t197 + t173 * t244) * t212 + t228 * t148;
	t143 = t228 * t147 - t172 * t184 + t173 * t202;
	t141 = t226 * t260 + (t202 * t230 - t163 * t194 + (t162 * t202 + t182 * t188 + t184 * t187) * t195) * t174;
	t139 = t260 * t261 + (t224 * t239 + (t230 * t244 + t191 * t194 + (t187 * t197 + (t162 * t216 - t182 * t241) * t208) * t195) * t212) * t174;
	t1 = [0, t139, t141, 0, 0, 0; 0, (-t144 * t258 + t153 * t248) * t237 + ((-t193 * t212 - t199 * t239) * t153 + (-t256 + t262) * t144 + (t248 * t142 + (-t139 * t182 - t148 * t162 + (-t212 * t241 + t216 * t239) * t208 + (-t148 * t201 - t197 * t212) * t146) * t251 + (-t197 * t239 - t139 * t201 - t148 * t187 + t191 * t212 + (t148 * t182 - t212 * t244) * t146) * t252) * t154) * t149, (-t143 * t258 - t153 * t186) * t237 + (t143 * t262 + t165 * t153 + (-t186 * t142 - t143 * t164 + (-t141 * t182 - t147 * t162 + t188 + (-t147 * t201 - t184) * t146) * t251 + (-t141 * t201 - t147 * t187 - t163 + (t147 * t182 - t202) * t146) * t252) * t154) * t149, 0, 0, 0; 0, 0.2e1 * (-t167 * t178 + t179 * t254) * t259 + (0.2e1 * t179 * t234 - t229 * t167 * t214 + t223 * t255 + (-t229 * t170 * t211 - t179 * t157 - t178 * t158 - t223 * t253) * t168) * t159, -t227 * t225 * t236 + (t227 * t164 - ((-qJD(4) * t167 - 0.2e1 * t234) * t214 + (t157 * t214 + (t158 - t238) * t211) * t168) * t225) * t159, t236 + 0.2e1 * (t157 * t168 * t159 + (-t159 * t257 - t168 * t259) * t170) * t170, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:05:50
	% EndTime: 2019-10-09 23:05:51
	% DurationCPUTime: 1.38s
	% Computational Cost: add. (3659->111), mult. (9671->227), div. (577->12), fcn. (12382->13), ass. (0->107)
	t241 = sin(pkin(11));
	t243 = cos(pkin(11));
	t246 = sin(qJ(2));
	t244 = cos(pkin(6));
	t248 = cos(qJ(2));
	t275 = t244 * t248;
	t227 = -t241 * t246 + t243 * t275;
	t220 = t227 * qJD(2);
	t276 = t244 * t246;
	t228 = t241 * t248 + t243 * t276;
	t245 = sin(qJ(3));
	t242 = sin(pkin(6));
	t279 = t242 * t245;
	t266 = t243 * t279;
	t247 = cos(qJ(3));
	t272 = qJD(3) * t247;
	t198 = -qJD(3) * t266 + t220 * t245 + t228 * t272;
	t278 = t242 * t247;
	t212 = t228 * t245 + t243 * t278;
	t210 = t212 ^ 2;
	t231 = -t244 * t247 + t246 * t279;
	t225 = 0.1e1 / t231 ^ 2;
	t206 = t210 * t225 + 0.1e1;
	t204 = 0.1e1 / t206;
	t232 = t244 * t245 + t246 * t278;
	t273 = qJD(2) * t248;
	t265 = t242 * t273;
	t217 = t232 * qJD(3) + t245 * t265;
	t224 = 0.1e1 / t231;
	t283 = t212 * t225;
	t176 = (-t198 * t224 + t217 * t283) * t204;
	t207 = atan2(-t212, t231);
	t202 = sin(t207);
	t203 = cos(t207);
	t260 = -t202 * t231 - t203 * t212;
	t172 = t260 * t176 - t198 * t202 + t203 * t217;
	t188 = -t202 * t212 + t203 * t231;
	t185 = 0.1e1 / t188;
	t186 = 0.1e1 / t188 ^ 2;
	t296 = t172 * t185 * t186;
	t267 = t241 * t276;
	t230 = t243 * t248 - t267;
	t257 = -t230 * t245 + t241 * t278;
	t295 = -0.2e1 * t257 * t296;
	t277 = t242 * t248;
	t256 = -t224 * t227 + t277 * t283;
	t294 = t245 * t256;
	t282 = t217 * t224 * t225;
	t293 = -0.2e1 * (t198 * t283 - t210 * t282) / t206 ^ 2;
	t216 = t230 * t247 + t241 * t279;
	t229 = t241 * t275 + t243 * t246;
	t240 = qJ(4) + qJ(5);
	t237 = sin(t240);
	t238 = cos(t240);
	t197 = t216 * t238 + t229 * t237;
	t193 = 0.1e1 / t197;
	t194 = 0.1e1 / t197 ^ 2;
	t223 = -qJD(2) * t267 + t243 * t273;
	t239 = qJD(4) + qJD(5);
	t262 = t216 * t239 - t223;
	t222 = t229 * qJD(2);
	t201 = t257 * qJD(3) - t222 * t247;
	t263 = t229 * t239 + t201;
	t183 = t263 * t237 + t262 * t238;
	t196 = t216 * t237 - t229 * t238;
	t192 = t196 ^ 2;
	t191 = t192 * t194 + 0.1e1;
	t288 = t194 * t196;
	t184 = -t262 * t237 + t263 * t238;
	t291 = t184 * t193 * t194;
	t292 = (t183 * t288 - t192 * t291) / t191 ^ 2;
	t290 = t186 * t257;
	t289 = t193 * t237;
	t287 = t196 * t238;
	t200 = t216 * qJD(3) - t222 * t245;
	t286 = t200 * t186;
	t285 = t202 * t257;
	t284 = t203 * t257;
	t281 = t229 * t245;
	t280 = t229 * t247;
	t274 = qJD(2) * t246;
	t211 = t257 ^ 2;
	t182 = t186 * t211 + 0.1e1;
	t271 = 0.2e1 * (-t211 * t296 - t257 * t286) / t182 ^ 2;
	t270 = -0.2e1 * t292;
	t268 = t196 * t291;
	t264 = -0.2e1 * t212 * t282;
	t261 = t239 * t280 - t222;
	t259 = t194 * t287 - t289;
	t214 = t228 * t247 - t266;
	t258 = -t214 * t224 + t232 * t283;
	t255 = qJD(3) * t281 - t223 * t247 + t230 * t239;
	t221 = t228 * qJD(2);
	t218 = -t231 * qJD(3) + t247 * t265;
	t209 = t230 * t237 - t238 * t280;
	t208 = -t230 * t238 - t237 * t280;
	t199 = -t212 * qJD(3) + t220 * t247;
	t189 = 0.1e1 / t191;
	t179 = 0.1e1 / t182;
	t178 = t204 * t294;
	t177 = t258 * t204;
	t174 = (-t202 * t227 + t203 * t277) * t245 + t260 * t178;
	t173 = t260 * t177 - t202 * t214 + t203 * t232;
	t171 = t258 * t293 + (t232 * t264 - t199 * t224 + (t198 * t232 + t212 * t218 + t214 * t217) * t225) * t204;
	t169 = t293 * t294 + (t256 * t272 + (t264 * t277 + t221 * t224 + (t217 * t227 + (t198 * t248 - t212 * t274) * t242) * t225) * t245) * t204;
	t168 = t270 + 0.2e1 * (t183 * t189 * t194 + (-t189 * t291 - t194 * t292) * t196) * t196;
	t1 = [0, t169, t171, 0, 0, 0; 0, (-t174 * t290 + t185 * t281) * t271 + ((-t223 * t245 - t229 * t272) * t185 + (-t286 + t295) * t174 + (t281 * t172 + (-t169 * t212 - t178 * t198 + (-t245 * t274 + t248 * t272) * t242 + (-t178 * t231 - t227 * t245) * t176) * t284 + (-t227 * t272 - t169 * t231 - t178 * t217 + t221 * t245 + (t178 * t212 - t245 * t277) * t176) * t285) * t186) * t179, (-t173 * t290 - t185 * t216) * t271 + (t173 * t295 + t201 * t185 + (-t216 * t172 - t173 * t200 + (-t171 * t212 - t177 * t198 + t218 + (-t177 * t231 - t214) * t176) * t284 + (-t171 * t231 - t177 * t217 - t199 + (t177 * t212 - t232) * t176) * t285) * t186) * t179, 0, 0, 0; 0, 0.2e1 * (-t193 * t208 + t209 * t288) * t292 + (0.2e1 * t209 * t268 - t261 * t193 * t238 + t255 * t289 + (-t261 * t196 * t237 - t209 * t183 - t208 * t184 - t255 * t287) * t194) * t189, -t259 * t257 * t270 + (t259 * t200 - ((-t193 * t239 - 0.2e1 * t268) * t238 + (t183 * t238 + (-t196 * t239 + t184) * t237) * t194) * t257) * t189, t168, t168, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:05:50
	% EndTime: 2019-10-09 23:05:51
	% DurationCPUTime: 1.42s
	% Computational Cost: add. (3659->111), mult. (9671->227), div. (577->12), fcn. (12382->13), ass. (0->107)
	t248 = sin(pkin(11));
	t250 = cos(pkin(11));
	t253 = sin(qJ(2));
	t251 = cos(pkin(6));
	t255 = cos(qJ(2));
	t282 = t251 * t255;
	t234 = -t248 * t253 + t250 * t282;
	t227 = t234 * qJD(2);
	t283 = t251 * t253;
	t235 = t248 * t255 + t250 * t283;
	t252 = sin(qJ(3));
	t249 = sin(pkin(6));
	t286 = t249 * t252;
	t273 = t250 * t286;
	t254 = cos(qJ(3));
	t279 = qJD(3) * t254;
	t205 = -qJD(3) * t273 + t227 * t252 + t235 * t279;
	t285 = t249 * t254;
	t219 = t235 * t252 + t250 * t285;
	t217 = t219 ^ 2;
	t238 = -t251 * t254 + t253 * t286;
	t232 = 0.1e1 / t238 ^ 2;
	t213 = t217 * t232 + 0.1e1;
	t211 = 0.1e1 / t213;
	t239 = t251 * t252 + t253 * t285;
	t280 = qJD(2) * t255;
	t272 = t249 * t280;
	t224 = qJD(3) * t239 + t252 * t272;
	t231 = 0.1e1 / t238;
	t290 = t219 * t232;
	t183 = (-t205 * t231 + t224 * t290) * t211;
	t214 = atan2(-t219, t238);
	t209 = sin(t214);
	t210 = cos(t214);
	t267 = -t209 * t238 - t210 * t219;
	t179 = t183 * t267 - t205 * t209 + t210 * t224;
	t195 = -t209 * t219 + t210 * t238;
	t192 = 0.1e1 / t195;
	t193 = 0.1e1 / t195 ^ 2;
	t303 = t179 * t192 * t193;
	t284 = t249 * t255;
	t263 = -t231 * t234 + t284 * t290;
	t302 = t252 * t263;
	t274 = t248 * t283;
	t237 = t250 * t255 - t274;
	t264 = -t237 * t252 + t248 * t285;
	t301 = -0.2e1 * t264 * t303;
	t289 = t224 * t231 * t232;
	t300 = -0.2e1 * (t205 * t290 - t217 * t289) / t213 ^ 2;
	t223 = t237 * t254 + t248 * t286;
	t236 = t248 * t282 + t250 * t253;
	t247 = qJ(4) + qJ(5);
	t244 = sin(t247);
	t245 = cos(t247);
	t204 = t223 * t245 + t236 * t244;
	t200 = 0.1e1 / t204;
	t201 = 0.1e1 / t204 ^ 2;
	t230 = -qJD(2) * t274 + t250 * t280;
	t246 = qJD(4) + qJD(5);
	t269 = t223 * t246 - t230;
	t229 = t236 * qJD(2);
	t208 = qJD(3) * t264 - t229 * t254;
	t270 = t236 * t246 + t208;
	t190 = t244 * t270 + t245 * t269;
	t203 = t223 * t244 - t236 * t245;
	t199 = t203 ^ 2;
	t198 = t199 * t201 + 0.1e1;
	t295 = t201 * t203;
	t191 = -t244 * t269 + t245 * t270;
	t298 = t191 * t200 * t201;
	t299 = (t190 * t295 - t199 * t298) / t198 ^ 2;
	t297 = t193 * t264;
	t296 = t200 * t244;
	t294 = t203 * t245;
	t207 = qJD(3) * t223 - t229 * t252;
	t293 = t207 * t193;
	t292 = t209 * t264;
	t291 = t210 * t264;
	t288 = t236 * t252;
	t287 = t236 * t254;
	t281 = qJD(2) * t253;
	t218 = t264 ^ 2;
	t189 = t193 * t218 + 0.1e1;
	t278 = 0.2e1 * (-t218 * t303 - t264 * t293) / t189 ^ 2;
	t277 = -0.2e1 * t299;
	t275 = t203 * t298;
	t271 = -0.2e1 * t219 * t289;
	t268 = t246 * t287 - t229;
	t266 = t201 * t294 - t296;
	t221 = t235 * t254 - t273;
	t265 = -t221 * t231 + t239 * t290;
	t262 = qJD(3) * t288 - t230 * t254 + t237 * t246;
	t228 = t235 * qJD(2);
	t225 = -qJD(3) * t238 + t254 * t272;
	t216 = t237 * t244 - t245 * t287;
	t215 = -t237 * t245 - t244 * t287;
	t206 = -qJD(3) * t219 + t227 * t254;
	t196 = 0.1e1 / t198;
	t186 = 0.1e1 / t189;
	t185 = t211 * t302;
	t184 = t265 * t211;
	t181 = (-t209 * t234 + t210 * t284) * t252 + t267 * t185;
	t180 = t184 * t267 - t209 * t221 + t210 * t239;
	t178 = t265 * t300 + (t239 * t271 - t206 * t231 + (t205 * t239 + t219 * t225 + t221 * t224) * t232) * t211;
	t176 = t300 * t302 + (t263 * t279 + (t271 * t284 + t228 * t231 + (t224 * t234 + (t205 * t255 - t219 * t281) * t249) * t232) * t252) * t211;
	t175 = t277 + 0.2e1 * (t190 * t196 * t201 + (-t196 * t298 - t201 * t299) * t203) * t203;
	t1 = [0, t176, t178, 0, 0, 0; 0, (-t181 * t297 + t192 * t288) * t278 + ((-t230 * t252 - t236 * t279) * t192 + (-t293 + t301) * t181 + (t288 * t179 + (-t176 * t219 - t185 * t205 + (-t252 * t281 + t255 * t279) * t249 + (-t185 * t238 - t234 * t252) * t183) * t291 + (-t234 * t279 - t176 * t238 - t185 * t224 + t228 * t252 + (t185 * t219 - t252 * t284) * t183) * t292) * t193) * t186, (-t180 * t297 - t192 * t223) * t278 + (t180 * t301 + t208 * t192 + (-t223 * t179 - t180 * t207 + (-t178 * t219 - t184 * t205 + t225 + (-t184 * t238 - t221) * t183) * t291 + (-t178 * t238 - t184 * t224 - t206 + (t184 * t219 - t239) * t183) * t292) * t193) * t186, 0, 0, 0; 0, 0.2e1 * (-t200 * t215 + t216 * t295) * t299 + (0.2e1 * t216 * t275 - t268 * t200 * t245 + t262 * t296 + (-t203 * t244 * t268 - t216 * t190 - t215 * t191 - t262 * t294) * t201) * t196, -t266 * t264 * t277 + (t266 * t207 - ((-t200 * t246 - 0.2e1 * t275) * t245 + (t190 * t245 + (-t203 * t246 + t191) * t244) * t201) * t264) * t196, t175, t175, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end