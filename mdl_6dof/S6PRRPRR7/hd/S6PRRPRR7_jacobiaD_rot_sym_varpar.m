% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR7
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
%   Wie in S6PRRPRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
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
	% StartTime: 2019-10-09 22:37:19
	% EndTime: 2019-10-09 22:37:19
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
	% StartTime: 2019-10-09 22:37:19
	% EndTime: 2019-10-09 22:37:20
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (2404->89), mult. (7372->192), div. (524->12), fcn. (9540->11), ass. (0->86)
	t161 = sin(pkin(11));
	t163 = cos(pkin(11));
	t166 = sin(qJ(2));
	t164 = cos(pkin(6));
	t168 = cos(qJ(2));
	t189 = t164 * t168;
	t151 = -t161 * t166 + t163 * t189;
	t141 = t151 * qJD(2);
	t190 = t164 * t166;
	t152 = t161 * t168 + t163 * t190;
	t165 = sin(qJ(3));
	t162 = sin(pkin(6));
	t193 = t162 * t165;
	t181 = t163 * t193;
	t167 = cos(qJ(3));
	t186 = qJD(3) * t167;
	t118 = -qJD(3) * t181 + t141 * t165 + t152 * t186;
	t192 = t162 * t167;
	t134 = t152 * t165 + t163 * t192;
	t131 = t134 ^ 2;
	t155 = -t164 * t167 + t166 * t193;
	t149 = 0.1e1 / t155 ^ 2;
	t129 = t131 * t149 + 0.1e1;
	t126 = 0.1e1 / t129;
	t156 = t164 * t165 + t166 * t192;
	t187 = qJD(2) * t168;
	t180 = t162 * t187;
	t139 = t156 * qJD(3) + t165 * t180;
	t148 = 0.1e1 / t155;
	t197 = t134 * t149;
	t106 = (-t118 * t148 + t139 * t197) * t126;
	t130 = atan2(-t134, t155);
	t122 = sin(t130);
	t123 = cos(t130);
	t178 = -t122 * t155 - t123 * t134;
	t103 = t178 * t106 - t122 * t118 + t123 * t139;
	t117 = -t122 * t134 + t123 * t155;
	t114 = 0.1e1 / t117;
	t115 = 0.1e1 / t117 ^ 2;
	t205 = t103 * t114 * t115;
	t182 = t161 * t190;
	t154 = t163 * t168 - t182;
	t176 = -t154 * t165 + t161 * t192;
	t204 = -0.2e1 * t176 * t205;
	t191 = t162 * t168;
	t175 = -t148 * t151 + t191 * t197;
	t203 = t165 * t175;
	t195 = t139 * t148 * t149;
	t202 = -0.2e1 * (t118 * t197 - t131 * t195) / t129 ^ 2;
	t153 = t161 * t189 + t163 * t166;
	t145 = 0.1e1 / t153;
	t146 = 0.1e1 / t153 ^ 2;
	t201 = t115 * t176;
	t138 = t154 * t167 + t161 * t193;
	t143 = t153 * qJD(2);
	t120 = t138 * qJD(3) - t143 * t165;
	t200 = t120 * t115;
	t199 = t122 * t176;
	t198 = t123 * t176;
	t196 = t138 * t154;
	t194 = t153 * t165;
	t188 = qJD(2) * t166;
	t132 = t176 ^ 2;
	t112 = t132 * t115 + 0.1e1;
	t185 = 0.2e1 * (-t132 * t205 - t176 * t200) / t112 ^ 2;
	t121 = t176 * qJD(3) - t143 * t167;
	t133 = t138 ^ 2;
	t128 = t133 * t146 + 0.1e1;
	t144 = -qJD(2) * t182 + t163 * t187;
	t147 = t145 * t146;
	t184 = 0.2e1 * (t138 * t146 * t121 - t133 * t147 * t144) / t128 ^ 2;
	t179 = -0.2e1 * t134 * t195;
	t136 = t152 * t167 - t181;
	t177 = -t136 * t148 + t156 * t197;
	t142 = t152 * qJD(2);
	t140 = -t155 * qJD(3) + t167 * t180;
	t124 = 0.1e1 / t128;
	t119 = -t134 * qJD(3) + t141 * t167;
	t109 = 0.1e1 / t112;
	t108 = t126 * t203;
	t107 = t177 * t126;
	t105 = (-t122 * t151 + t123 * t191) * t165 + t178 * t108;
	t104 = t178 * t107 - t122 * t136 + t123 * t156;
	t102 = t177 * t202 + (t156 * t179 - t119 * t148 + (t118 * t156 + t134 * t140 + t136 * t139) * t149) * t126;
	t100 = t202 * t203 + (t175 * t186 + (t179 * t191 + t142 * t148 + (t139 * t151 + (t118 * t168 - t134 * t188) * t162) * t149) * t165) * t126;
	t1 = [0, t100, t102, 0, 0, 0; 0, (-t105 * t201 + t114 * t194) * t185 + ((-t144 * t165 - t153 * t186) * t114 + (-t200 + t204) * t105 + (t194 * t103 + (-t100 * t134 - t108 * t118 + (-t165 * t188 + t168 * t186) * t162 + (-t108 * t155 - t151 * t165) * t106) * t198 + (-t151 * t186 - t100 * t155 - t108 * t139 + t142 * t165 + (t108 * t134 - t165 * t191) * t106) * t199) * t115) * t109, (-t104 * t201 - t114 * t138) * t185 + (t104 * t204 + t121 * t114 + (-t138 * t103 - t104 * t120 + (-t102 * t134 - t107 * t118 + t140 + (-t107 * t155 - t136) * t106) * t198 + (-t102 * t155 - t107 * t139 - t119 + (t107 * t134 - t156) * t106) * t199) * t115) * t109, 0, 0, 0; 0, (t145 * t153 * t167 + t146 * t196) * t184 + (qJD(3) * t145 * t194 + (-t121 * t154 + t138 * t143) * t146 + (0.2e1 * t147 * t196 + (t146 * t153 - t145) * t167) * t144) * t124, -t176 * t145 * t184 + (-t144 * t146 * t176 - t120 * t145) * t124, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:19
	% EndTime: 2019-10-09 22:37:20
	% DurationCPUTime: 1.41s
	% Computational Cost: add. (3002->107), mult. (9085->220), div. (559->12), fcn. (11668->13), ass. (0->104)
	t212 = cos(pkin(6));
	t218 = cos(qJ(2));
	t268 = cos(pkin(11));
	t239 = t268 * t218;
	t210 = sin(pkin(11));
	t215 = sin(qJ(2));
	t253 = t210 * t215;
	t201 = t212 * t239 - t253;
	t194 = t201 * qJD(2);
	t217 = cos(qJ(3));
	t214 = sin(qJ(3));
	t240 = t268 * t215;
	t252 = t210 * t218;
	t227 = -t212 * t240 - t252;
	t211 = sin(pkin(6));
	t241 = t211 * t268;
	t270 = t227 * t214 - t217 * t241;
	t166 = t270 * qJD(3) + t194 * t217;
	t187 = -t214 * t241 - t227 * t217;
	t184 = t187 ^ 2;
	t250 = t211 * t217;
	t205 = t212 * t214 + t215 * t250;
	t199 = 0.1e1 / t205 ^ 2;
	t179 = t184 * t199 + 0.1e1;
	t177 = 0.1e1 / t179;
	t251 = t211 * t214;
	t204 = t212 * t217 - t215 * t251;
	t249 = t211 * t218;
	t242 = qJD(2) * t249;
	t192 = t204 * qJD(3) + t217 * t242;
	t198 = 0.1e1 / t205;
	t258 = t187 * t199;
	t149 = (-t166 * t198 + t192 * t258) * t177;
	t180 = atan2(-t187, t205);
	t175 = sin(t180);
	t176 = cos(t180);
	t234 = -t175 * t205 - t176 * t187;
	t145 = t234 * t149 - t175 * t166 + t176 * t192;
	t159 = -t175 * t187 + t176 * t205;
	t156 = 0.1e1 / t159;
	t157 = 0.1e1 / t159 ^ 2;
	t274 = t145 * t156 * t157;
	t213 = sin(qJ(5));
	t216 = cos(qJ(5));
	t228 = -t212 * t252 - t240;
	t203 = -t212 * t253 + t239;
	t230 = -t203 * t214 + t210 * t250;
	t233 = t213 * t228 - t216 * t230;
	t273 = t233 * qJD(5);
	t190 = t203 * t217 + t210 * t251;
	t272 = 0.2e1 * t190 * t274;
	t229 = -t198 * t201 + t249 * t258;
	t271 = t217 * t229;
	t257 = t192 * t198 * t199;
	t269 = -0.2e1 * (t166 * t258 - t184 * t257) / t179 ^ 2;
	t255 = t228 * t216;
	t174 = -t213 * t230 - t255;
	t170 = 0.1e1 / t174;
	t171 = 0.1e1 / t174 ^ 2;
	t196 = t228 * qJD(2);
	t167 = t190 * qJD(3) + t196 * t214;
	t197 = t203 * qJD(2);
	t160 = t174 * qJD(5) - t167 * t216 + t197 * t213;
	t169 = t233 ^ 2;
	t164 = t169 * t171 + 0.1e1;
	t262 = t171 * t233;
	t161 = t167 * t213 + t197 * t216 + t273;
	t265 = t161 * t170 * t171;
	t267 = (-t160 * t262 - t169 * t265) / t164 ^ 2;
	t266 = t157 * t190;
	t168 = t230 * qJD(3) + t196 * t217;
	t264 = t168 * t157;
	t263 = t170 * t216;
	t261 = t233 * t213;
	t260 = t175 * t190;
	t259 = t176 * t190;
	t256 = t228 * t214;
	t254 = t228 * t217;
	t248 = qJD(2) * t215;
	t247 = qJD(3) * t214;
	t185 = t190 ^ 2;
	t155 = t185 * t157 + 0.1e1;
	t246 = 0.2e1 * (-t185 * t274 + t190 * t264) / t155 ^ 2;
	t245 = 0.2e1 * t267;
	t238 = -0.2e1 * t233 * t265;
	t237 = -0.2e1 * t187 * t257;
	t235 = qJD(5) * t256 + t196;
	t232 = -t171 * t261 + t263;
	t231 = -t198 * t270 + t204 * t258;
	t226 = -qJD(3) * t254 + qJD(5) * t203 + t197 * t214;
	t195 = t227 * qJD(2);
	t191 = -t205 * qJD(3) - t214 * t242;
	t182 = t203 * t216 + t213 * t256;
	t181 = t203 * t213 - t214 * t255;
	t165 = t187 * qJD(3) + t194 * t214;
	t162 = 0.1e1 / t164;
	t152 = 0.1e1 / t155;
	t151 = t177 * t271;
	t150 = t231 * t177;
	t147 = (-t175 * t201 + t176 * t249) * t217 + t234 * t151;
	t146 = t234 * t150 - t175 * t270 + t176 * t204;
	t144 = t231 * t269 + (t204 * t237 + t165 * t198 + (t166 * t204 + t187 * t191 + t192 * t270) * t199) * t177;
	t142 = t269 * t271 + (-t229 * t247 + (t237 * t249 - t195 * t198 + (t192 * t201 + (t166 * t218 - t187 * t248) * t211) * t199) * t217) * t177;
	t1 = [0, t142, t144, 0, 0, 0; 0, (t147 * t266 - t156 * t254) * t246 + ((-t197 * t217 - t228 * t247) * t156 + (-t264 + t272) * t147 + (-t254 * t145 - (-t142 * t187 - t151 * t166 + (-t217 * t248 - t218 * t247) * t211 + (-t151 * t205 - t201 * t217) * t149) * t259 - (t201 * t247 - t142 * t205 - t151 * t192 - t195 * t217 + (t151 * t187 - t217 * t249) * t149) * t260) * t157) * t152, (t146 * t266 - t156 * t230) * t246 + (t146 * t272 - t167 * t156 + (-t230 * t145 - t146 * t168 - (-t144 * t187 - t150 * t166 + t191 + (-t150 * t205 - t270) * t149) * t259 - (-t144 * t205 - t150 * t192 + t165 + (t150 * t187 - t204) * t149) * t260) * t157) * t152, 0, 0, 0; 0, (-t170 * t181 - t182 * t262) * t245 + (t182 * t238 + t235 * t170 * t213 + t226 * t263 + (t216 * t233 * t235 - t182 * t160 - t181 * t161 - t226 * t261) * t171) * t162, t232 * t190 * t245 + (-t232 * t168 + ((qJD(5) * t170 + t238) * t213 + (-t160 * t213 + (t161 + t273) * t216) * t171) * t190) * t162, 0, -0.2e1 * t267 - 0.2e1 * (t160 * t171 * t162 - (-t162 * t265 - t171 * t267) * t233) * t233, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:19
	% EndTime: 2019-10-09 22:37:20
	% DurationCPUTime: 1.40s
	% Computational Cost: add. (3659->109), mult. (9671->221), div. (577->12), fcn. (12382->13), ass. (0->107)
	t245 = cos(pkin(6));
	t249 = cos(qJ(2));
	t299 = cos(pkin(11));
	t271 = t299 * t249;
	t243 = sin(pkin(11));
	t247 = sin(qJ(2));
	t285 = t243 * t247;
	t230 = t245 * t271 - t285;
	t223 = t230 * qJD(2);
	t248 = cos(qJ(3));
	t246 = sin(qJ(3));
	t272 = t299 * t247;
	t284 = t243 * t249;
	t258 = -t245 * t272 - t284;
	t244 = sin(pkin(6));
	t273 = t244 * t299;
	t301 = t258 * t246 - t248 * t273;
	t201 = t301 * qJD(3) + t223 * t248;
	t216 = -t246 * t273 - t258 * t248;
	t213 = t216 ^ 2;
	t282 = t244 * t248;
	t234 = t245 * t246 + t247 * t282;
	t228 = 0.1e1 / t234 ^ 2;
	t208 = t213 * t228 + 0.1e1;
	t206 = 0.1e1 / t208;
	t283 = t244 * t246;
	t233 = t245 * t248 - t247 * t283;
	t281 = t244 * t249;
	t274 = qJD(2) * t281;
	t221 = t233 * qJD(3) + t248 * t274;
	t227 = 0.1e1 / t234;
	t289 = t216 * t228;
	t178 = (-t201 * t227 + t221 * t289) * t206;
	t209 = atan2(-t216, t234);
	t204 = sin(t209);
	t205 = cos(t209);
	t264 = -t204 * t234 - t205 * t216;
	t174 = t264 * t178 - t201 * t204 + t205 * t221;
	t190 = -t204 * t216 + t205 * t234;
	t187 = 0.1e1 / t190;
	t188 = 0.1e1 / t190 ^ 2;
	t304 = t174 * t187 * t188;
	t232 = -t245 * t285 + t271;
	t219 = t232 * t248 + t243 * t283;
	t303 = 0.2e1 * t219 * t304;
	t260 = -t227 * t230 + t281 * t289;
	t302 = t248 * t260;
	t288 = t221 * t227 * t228;
	t300 = -0.2e1 * (t201 * t289 - t213 * t288) / t208 ^ 2;
	t242 = qJ(5) + qJ(6);
	t239 = sin(t242);
	t240 = cos(t242);
	t259 = -t245 * t284 - t272;
	t261 = -t232 * t246 + t243 * t282;
	t199 = -t239 * t261 - t240 * t259;
	t195 = 0.1e1 / t199;
	t196 = 0.1e1 / t199 ^ 2;
	t226 = t232 * qJD(2);
	t241 = qJD(5) + qJD(6);
	t267 = -t241 * t261 + t226;
	t225 = t259 * qJD(2);
	t202 = t219 * qJD(3) + t225 * t246;
	t268 = -t241 * t259 - t202;
	t185 = t267 * t239 + t268 * t240;
	t198 = -t239 * t259 + t240 * t261;
	t194 = t198 ^ 2;
	t193 = t194 * t196 + 0.1e1;
	t294 = t196 * t198;
	t186 = -t268 * t239 + t267 * t240;
	t297 = t186 * t195 * t196;
	t298 = (t185 * t294 - t194 * t297) / t193 ^ 2;
	t296 = t188 * t219;
	t295 = t195 * t240;
	t293 = t198 * t239;
	t203 = t261 * qJD(3) + t225 * t248;
	t292 = t203 * t188;
	t291 = t204 * t219;
	t290 = t205 * t219;
	t287 = t259 * t246;
	t286 = t259 * t248;
	t280 = qJD(2) * t247;
	t279 = qJD(3) * t246;
	t214 = t219 ^ 2;
	t184 = t188 * t214 + 0.1e1;
	t278 = 0.2e1 * (-t214 * t304 + t219 * t292) / t184 ^ 2;
	t277 = 0.2e1 * t298;
	t270 = 0.2e1 * t198 * t297;
	t269 = -0.2e1 * t216 * t288;
	t265 = t241 * t287 + t225;
	t263 = t196 * t293 + t295;
	t262 = -t227 * t301 + t233 * t289;
	t257 = -qJD(3) * t286 + t226 * t246 + t232 * t241;
	t224 = t258 * qJD(2);
	t220 = -t234 * qJD(3) - t246 * t274;
	t211 = t232 * t240 + t239 * t287;
	t210 = t232 * t239 - t240 * t287;
	t200 = t216 * qJD(3) + t223 * t246;
	t191 = 0.1e1 / t193;
	t181 = 0.1e1 / t184;
	t180 = t206 * t302;
	t179 = t262 * t206;
	t176 = (-t204 * t230 + t205 * t281) * t248 + t264 * t180;
	t175 = t264 * t179 - t204 * t301 + t205 * t233;
	t173 = t262 * t300 + (t233 * t269 + t200 * t227 + (t201 * t233 + t216 * t220 + t221 * t301) * t228) * t206;
	t171 = t300 * t302 + (-t260 * t279 + (t269 * t281 - t224 * t227 + (t221 * t230 + (t201 * t249 - t216 * t280) * t244) * t228) * t248) * t206;
	t170 = -0.2e1 * t298 + 0.2e1 * (t185 * t191 * t196 + (-t191 * t297 - t196 * t298) * t198) * t198;
	t1 = [0, t171, t173, 0, 0, 0; 0, (t176 * t296 - t187 * t286) * t278 + ((-t226 * t248 - t259 * t279) * t187 + (-t292 + t303) * t176 + (-t286 * t174 - (-t171 * t216 - t180 * t201 + (-t248 * t280 - t249 * t279) * t244 + (-t180 * t234 - t230 * t248) * t178) * t290 - (t230 * t279 - t171 * t234 - t180 * t221 - t224 * t248 + (t180 * t216 - t248 * t281) * t178) * t291) * t188) * t181, (t175 * t296 - t187 * t261) * t278 + (t175 * t303 - t202 * t187 + (-t261 * t174 - t175 * t203 - (-t173 * t216 - t179 * t201 + t220 + (-t179 * t234 - t301) * t178) * t290 - (-t173 * t234 - t179 * t221 + t200 + (t179 * t216 - t233) * t178) * t291) * t188) * t181, 0, 0, 0; 0, (-t195 * t210 + t211 * t294) * t277 + (t211 * t270 + t265 * t195 * t239 + t257 * t295 + (-t265 * t198 * t240 - t211 * t185 - t210 * t186 + t257 * t293) * t196) * t191, t263 * t219 * t277 + (-t263 * t203 + ((t195 * t241 + t270) * t239 + (-t185 * t239 + (-t198 * t241 + t186) * t240) * t196) * t219) * t191, 0, t170, t170;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end