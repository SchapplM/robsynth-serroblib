% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPP1
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
%   Wie in S6PRRRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
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
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:09
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (756->55), mult. (2271->133), div. (423->14), fcn. (2956->11), ass. (0->65)
	t139 = sin(qJ(2));
	t141 = cos(qJ(2));
	t136 = sin(pkin(10));
	t166 = cos(pkin(6));
	t155 = t136 * t166;
	t165 = cos(pkin(10));
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
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:10
	% DurationCPUTime: 1.35s
	% Computational Cost: add. (3002->109), mult. (9085->226), div. (559->12), fcn. (11668->13), ass. (0->103)
	t207 = sin(pkin(10));
	t209 = cos(pkin(10));
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
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:10
	% DurationCPUTime: 1.41s
	% Computational Cost: add. (3313->110), mult. (9085->226), div. (559->12), fcn. (11668->13), ass. (0->103)
	t219 = sin(pkin(10));
	t221 = cos(pkin(10));
	t224 = sin(qJ(2));
	t222 = cos(pkin(6));
	t226 = cos(qJ(2));
	t252 = t222 * t226;
	t206 = -t219 * t224 + t221 * t252;
	t199 = t206 * qJD(2);
	t253 = t222 * t224;
	t207 = t219 * t226 + t221 * t253;
	t223 = sin(qJ(3));
	t220 = sin(pkin(6));
	t256 = t220 * t223;
	t242 = t221 * t256;
	t225 = cos(qJ(3));
	t249 = qJD(3) * t225;
	t177 = -qJD(3) * t242 + t199 * t223 + t207 * t249;
	t255 = t220 * t225;
	t191 = t207 * t223 + t221 * t255;
	t189 = t191 ^ 2;
	t210 = -t222 * t225 + t224 * t256;
	t204 = 0.1e1 / t210 ^ 2;
	t185 = t189 * t204 + 0.1e1;
	t183 = 0.1e1 / t185;
	t211 = t222 * t223 + t224 * t255;
	t250 = qJD(2) * t226;
	t241 = t220 * t250;
	t196 = t211 * qJD(3) + t223 * t241;
	t203 = 0.1e1 / t210;
	t260 = t191 * t204;
	t155 = (-t177 * t203 + t196 * t260) * t183;
	t186 = atan2(-t191, t210);
	t181 = sin(t186);
	t182 = cos(t186);
	t238 = -t181 * t210 - t182 * t191;
	t151 = t238 * t155 - t181 * t177 + t182 * t196;
	t167 = -t181 * t191 + t182 * t210;
	t164 = 0.1e1 / t167;
	t165 = 0.1e1 / t167 ^ 2;
	t272 = t151 * t164 * t165;
	t243 = t219 * t253;
	t209 = t221 * t226 - t243;
	t235 = -t209 * t223 + t219 * t255;
	t271 = -0.2e1 * t235 * t272;
	t254 = t220 * t226;
	t234 = -t203 * t206 + t254 * t260;
	t270 = t223 * t234;
	t259 = t196 * t203 * t204;
	t269 = -0.2e1 * (t177 * t260 - t189 * t259) / t185 ^ 2;
	t195 = t209 * t225 + t219 * t256;
	t208 = t219 * t252 + t221 * t224;
	t218 = qJ(4) + pkin(11);
	t216 = sin(t218);
	t217 = cos(t218);
	t176 = t195 * t217 + t208 * t216;
	t172 = 0.1e1 / t176;
	t173 = 0.1e1 / t176 ^ 2;
	t201 = t208 * qJD(2);
	t180 = t235 * qJD(3) - t201 * t225;
	t202 = -qJD(2) * t243 + t221 * t250;
	t162 = t176 * qJD(4) + t180 * t216 - t202 * t217;
	t175 = t195 * t216 - t208 * t217;
	t171 = t175 ^ 2;
	t170 = t171 * t173 + 0.1e1;
	t264 = t173 * t175;
	t248 = qJD(4) * t175;
	t163 = t180 * t217 + t202 * t216 - t248;
	t267 = t163 * t172 * t173;
	t268 = (t162 * t264 - t171 * t267) / t170 ^ 2;
	t266 = t165 * t235;
	t265 = t172 * t216;
	t263 = t175 * t217;
	t262 = t181 * t235;
	t261 = t182 * t235;
	t258 = t208 * t223;
	t257 = t208 * t225;
	t251 = qJD(2) * t224;
	t190 = t235 ^ 2;
	t161 = t165 * t190 + 0.1e1;
	t179 = t195 * qJD(3) - t201 * t223;
	t247 = 0.2e1 * (-t179 * t266 - t190 * t272) / t161 ^ 2;
	t246 = -0.2e1 * t268;
	t244 = t175 * t267;
	t240 = -0.2e1 * t191 * t259;
	t239 = qJD(4) * t257 - t201;
	t237 = t173 * t263 - t265;
	t193 = t207 * t225 - t242;
	t236 = -t193 * t203 + t211 * t260;
	t233 = qJD(3) * t258 + qJD(4) * t209 - t202 * t225;
	t200 = t207 * qJD(2);
	t197 = -t210 * qJD(3) + t225 * t241;
	t188 = t209 * t216 - t217 * t257;
	t187 = -t209 * t217 - t216 * t257;
	t178 = -t191 * qJD(3) + t199 * t225;
	t168 = 0.1e1 / t170;
	t158 = 0.1e1 / t161;
	t157 = t183 * t270;
	t156 = t236 * t183;
	t153 = (-t181 * t206 + t182 * t254) * t223 + t238 * t157;
	t152 = t238 * t156 - t181 * t193 + t182 * t211;
	t150 = t236 * t269 + (t211 * t240 - t178 * t203 + (t177 * t211 + t191 * t197 + t193 * t196) * t204) * t183;
	t148 = t269 * t270 + (t234 * t249 + (t240 * t254 + t200 * t203 + (t196 * t206 + (t177 * t226 - t191 * t251) * t220) * t204) * t223) * t183;
	t1 = [0, t148, t150, 0, 0, 0; 0, (-t153 * t266 + t164 * t258) * t247 + ((-t202 * t223 - t208 * t249) * t164 + t153 * t271 + (-t153 * t179 + t258 * t151 + (-t148 * t191 - t157 * t177 + (-t223 * t251 + t226 * t249) * t220 + (-t157 * t210 - t206 * t223) * t155) * t261 + (-t206 * t249 - t148 * t210 - t157 * t196 + t200 * t223 + (t157 * t191 - t223 * t254) * t155) * t262) * t165) * t158, (-t152 * t266 - t164 * t195) * t247 + (t152 * t271 + t180 * t164 + (-t195 * t151 - t152 * t179 + (-t150 * t191 - t156 * t177 + t197 + (-t156 * t210 - t193) * t155) * t261 + (-t150 * t210 - t156 * t196 - t178 + (t156 * t191 - t211) * t155) * t262) * t165) * t158, 0, 0, 0; 0, 0.2e1 * (-t172 * t187 + t188 * t264) * t268 + (0.2e1 * t188 * t244 - t239 * t172 * t217 + t233 * t265 + (-t239 * t175 * t216 - t188 * t162 - t187 * t163 - t233 * t263) * t173) * t168, -t237 * t235 * t246 + (t237 * t179 - ((-qJD(4) * t172 - 0.2e1 * t244) * t217 + (t162 * t217 + (t163 - t248) * t216) * t173) * t235) * t168, t246 + 0.2e1 * (t162 * t168 * t173 + (-t168 * t267 - t173 * t268) * t175) * t175, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:11
	% DurationCPUTime: 2.66s
	% Computational Cost: add. (11295->157), mult. (21125->308), div. (784->12), fcn. (27225->13), ass. (0->124)
	t261 = qJ(4) + pkin(11);
	t259 = sin(t261);
	t260 = cos(t261);
	t263 = sin(qJ(3));
	t265 = cos(qJ(3));
	t264 = sin(qJ(2));
	t266 = cos(qJ(2));
	t323 = cos(pkin(10));
	t324 = cos(pkin(6));
	t287 = t324 * t323;
	t322 = sin(pkin(10));
	t275 = -t264 * t287 - t322 * t266;
	t262 = sin(pkin(6));
	t292 = t262 * t323;
	t277 = t263 * t292 + t265 * t275;
	t285 = t322 * t264 - t266 * t287;
	t222 = t285 * t259 - t260 * t277;
	t238 = t263 * t275 - t265 * t292;
	t249 = t285 * qJD(2);
	t226 = t238 * qJD(3) - t249 * t265;
	t273 = t275 * qJD(2);
	t198 = t222 * qJD(4) + t226 * t259 + t260 * t273;
	t278 = t285 * t260;
	t220 = -t259 * t277 - t278;
	t215 = t220 ^ 2;
	t307 = t262 * t264;
	t256 = t324 * t263 + t265 * t307;
	t306 = t262 * t266;
	t235 = t256 * t259 + t260 * t306;
	t233 = 0.1e1 / t235 ^ 2;
	t207 = t215 * t233 + 0.1e1;
	t205 = 0.1e1 / t207;
	t236 = t256 * t260 - t259 * t306;
	t255 = -t263 * t307 + t324 * t265;
	t305 = qJD(2) * t262;
	t293 = t266 * t305;
	t243 = t255 * qJD(3) + t265 * t293;
	t294 = t264 * t305;
	t212 = t236 * qJD(4) + t243 * t259 - t260 * t294;
	t232 = 0.1e1 / t235;
	t313 = t220 * t233;
	t185 = (-t198 * t232 + t212 * t313) * t205;
	t208 = atan2(-t220, t235);
	t203 = sin(t208);
	t204 = cos(t208);
	t284 = -t203 * t235 - t204 * t220;
	t181 = t284 * t185 - t203 * t198 + t204 * t212;
	t197 = -t203 * t220 + t204 * t235;
	t194 = 0.1e1 / t197;
	t195 = 0.1e1 / t197 ^ 2;
	t328 = t181 * t194 * t195;
	t286 = t324 * t322;
	t274 = -t323 * t264 - t266 * t286;
	t254 = -t264 * t286 + t323 * t266;
	t291 = t262 * t322;
	t276 = -t254 * t265 - t263 * t291;
	t223 = -t259 * t276 + t260 * t274;
	t327 = -0.2e1 * t223;
	t289 = 0.2e1 * t223 * t328;
	t280 = -t232 * t238 + t255 * t313;
	t326 = t259 * t280;
	t315 = t212 * t232 * t233;
	t325 = -0.2e1 * (t198 * t313 - t215 * t315) / t207 ^ 2;
	t224 = -t259 * t274 - t260 * t276;
	t217 = 0.1e1 / t224;
	t218 = 0.1e1 / t224 ^ 2;
	t240 = -t254 * t263 + t265 * t291;
	t237 = t240 ^ 2;
	t311 = t237 * t218;
	t211 = 0.1e1 + t311;
	t250 = t274 * qJD(2);
	t228 = t240 * qJD(3) + t250 * t265;
	t251 = t254 * qJD(2);
	t201 = -t223 * qJD(4) + t228 * t260 + t251 * t259;
	t318 = t201 * t217 * t218;
	t295 = t237 * t318;
	t227 = t276 * qJD(3) - t250 * t263;
	t312 = t227 * t240;
	t321 = (t218 * t312 - t295) / t211 ^ 2;
	t320 = t195 * t223;
	t200 = t224 * qJD(4) + t228 * t259 - t251 * t260;
	t319 = t200 * t195;
	t317 = t203 * t223;
	t316 = t204 * t223;
	t314 = t218 * t240;
	t310 = t240 * t259;
	t309 = t259 * t265;
	t308 = t260 * t265;
	t304 = qJD(2) * t265;
	t303 = qJD(3) * t263;
	t302 = qJD(4) * t259;
	t301 = qJD(4) * t260;
	t300 = t265 * qJD(4);
	t216 = t223 ^ 2;
	t193 = t216 * t195 + 0.1e1;
	t299 = 0.2e1 * (-t216 * t328 + t223 * t319) / t193 ^ 2;
	t298 = 0.2e1 * t321;
	t296 = t240 * t318;
	t288 = -0.2e1 * t220 * t315;
	t282 = -t222 * t232 + t236 * t313;
	t279 = t265 * t285;
	t229 = -t259 * t279 + t260 * t275;
	t244 = (-t260 * t264 + t266 * t309) * t262;
	t281 = -t229 * t232 + t244 * t313;
	t242 = -t256 * qJD(3) - t263 * t293;
	t231 = t254 * t259 + t274 * t308;
	t230 = -t254 * t260 + t274 * t309;
	t225 = t277 * qJD(3) + t249 * t263;
	t214 = ((-qJD(2) + t300) * t266 * t260 + (-t266 * t303 + (qJD(4) - t304) * t264) * t259) * t262;
	t213 = -t235 * qJD(4) + t243 * t260 + t259 * t294;
	t209 = 0.1e1 / t211;
	t202 = t249 * t260 - t275 * t302 - t279 * t301 + (t275 * t304 + t285 * t303) * t259;
	t199 = qJD(4) * t278 + t226 * t260 - t259 * t273 + t277 * t302;
	t191 = 0.1e1 / t193;
	t190 = t205 * t326;
	t188 = t281 * t205;
	t187 = t282 * t205;
	t184 = (-t203 * t238 + t204 * t255) * t259 + t284 * t190;
	t183 = t284 * t188 - t203 * t229 + t204 * t244;
	t182 = t284 * t187 - t203 * t222 + t204 * t236;
	t180 = t281 * t325 + (t244 * t288 - t202 * t232 + (t198 * t244 + t212 * t229 + t214 * t220) * t233) * t205;
	t178 = t282 * t325 + (t236 * t288 - t199 * t232 + (t198 * t236 + t212 * t222 + t213 * t220) * t233) * t205;
	t177 = t325 * t326 + (t280 * t301 + (t255 * t288 - t225 * t232 + (t198 * t255 + t212 * t238 + t220 * t242) * t233) * t259) * t205;
	t1 = [0, t180, t177, t178, 0, 0; 0, (t183 * t320 - t194 * t230) * t299 + (t183 * t289 + (-t230 * t181 - t183 * t200 - (-t180 * t220 - t188 * t198 + t214 + (-t188 * t235 - t229) * t185) * t316 - (-t180 * t235 - t188 * t212 - t202 + (t188 * t220 - t244) * t185) * t317) * t195 + ((t274 * t300 - t250) * t260 + (qJD(4) * t254 - t251 * t265 - t274 * t303) * t259) * t194) * t191, (t184 * t320 - t194 * t310) * t299 + ((t227 * t259 + t240 * t301) * t194 + (-t319 + t289) * t184 + (-t310 * t181 - (t255 * t301 - t177 * t220 - t190 * t198 + t242 * t259 + (-t190 * t235 - t238 * t259) * t185) * t316 - (-t238 * t301 - t177 * t235 - t190 * t212 - t225 * t259 + (t190 * t220 - t255 * t259) * t185) * t317) * t195) * t191, (t182 * t320 - t194 * t224) * t299 + (t182 * t289 + t201 * t194 + (-t224 * t181 - t182 * t200 - (-t178 * t220 - t187 * t198 + t213 + (-t187 * t235 - t222) * t185) * t316 - (-t178 * t235 - t187 * t212 - t199 + (t187 * t220 - t236) * t185) * t317) * t195) * t191, 0, 0; 0, (t217 * t263 * t274 + t231 * t314) * t298 + (0.2e1 * t231 * t296 + (-qJD(3) * t265 * t274 + t251 * t263) * t217 + (-(t250 * t259 - t251 * t308 + t254 * t301) * t240 - t231 * t227 - (-t263 * t201 - (t259 * t300 + t260 * t303) * t240) * t274) * t218) * t209, (-t217 * t276 + t260 * t311) * t298 + (0.2e1 * t260 * t295 - t217 * t228 + (-t201 * t276 + t237 * t302 - 0.2e1 * t260 * t312) * t218) * t209, t314 * t321 * t327 + (t296 * t327 + (t200 * t240 + t223 * t227) * t218) * t209, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end