% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR3
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
%   Wie in S6PRPRRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:57
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:35
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
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (610->48), mult. (1820->126), div. (396->14), fcn. (2421->11), ass. (0->64)
	t116 = sin(qJ(2));
	t117 = cos(qJ(2));
	t113 = sin(pkin(11));
	t136 = cos(pkin(6));
	t129 = t113 * t136;
	t135 = cos(pkin(11));
	t104 = -t116 * t129 + t135 * t117;
	t126 = t136 * t135;
	t100 = t113 * t116 - t117 * t126;
	t114 = sin(pkin(6));
	t132 = t114 * t117;
	t90 = atan2(-t100, -t132);
	t88 = sin(t90);
	t89 = cos(t90);
	t77 = -t88 * t100 - t89 * t132;
	t74 = 0.1e1 / t77;
	t112 = sin(pkin(12));
	t115 = cos(pkin(12));
	t133 = t113 * t114;
	t87 = t104 * t115 + t112 * t133;
	t83 = 0.1e1 / t87;
	t109 = 0.1e1 / t117;
	t110 = 0.1e1 / t117 ^ 2;
	t75 = 0.1e1 / t77 ^ 2;
	t84 = 0.1e1 / t87 ^ 2;
	t131 = qJD(2) * t116;
	t137 = t117 * t88;
	t143 = t100 * t89;
	t134 = t110 * t116;
	t130 = t100 * t134;
	t107 = 0.1e1 / t114;
	t108 = 0.1e1 / t114 ^ 2;
	t98 = t100 ^ 2;
	t93 = t98 * t108 * t110 + 0.1e1;
	t91 = 0.1e1 / t93;
	t141 = t107 * t91;
	t102 = t113 * t117 + t116 * t126;
	t95 = t102 * qJD(2);
	t69 = (qJD(2) * t130 + t109 * t95) * t141;
	t66 = -t69 * t143 - t88 * t95 + (t89 * t131 + t69 * t137) * t114;
	t146 = t66 * t74 * t75;
	t139 = t112 * t84;
	t86 = t104 * t112 - t115 * t133;
	t82 = t86 ^ 2;
	t81 = t82 * t84 + 0.1e1;
	t85 = t83 * t84;
	t124 = -t135 * t116 - t117 * t129;
	t96 = t124 * qJD(2);
	t145 = (-t115 * t82 * t85 + t86 * t139) * t96 / t81 ^ 2;
	t125 = t102 * t109 + t130;
	t70 = t125 * t141;
	t144 = t69 * t70;
	t142 = t124 * t75;
	t140 = t112 * t83;
	t138 = t115 * t86;
	t111 = t109 * t110;
	t99 = t124 ^ 2;
	t97 = t104 * qJD(2);
	t94 = qJD(2) * t100;
	t79 = 0.1e1 / t81;
	t73 = t99 * t75 + 0.1e1;
	t67 = -t70 * t143 - t88 * t102 + (t116 * t89 + t70 * t137) * t114;
	t65 = (-0.2e1 * t125 / t93 ^ 2 * (t100 * t110 * t95 + t111 * t98 * t131) * t108 + (t95 * t134 - t109 * t94 + (t102 * t134 + (0.2e1 * t111 * t116 ^ 2 + t109) * t100) * qJD(2)) * t91) * t107;
	t1 = [0, t65, 0, 0, 0, 0; 0, 0.2e1 * (-t104 * t74 - t67 * t142) / t73 ^ 2 * (-t97 * t142 - t99 * t146) + (-0.2e1 * t67 * t124 * t146 + t96 * t74 + (-t104 * t66 - t67 * t97 - (-(t100 * t144 + t94) * t88 - (-t100 * t65 - t102 * t69 - t70 * t95) * t89) * t124) * t75 + ((-qJD(2) * t70 - t69) * t88 * t116 + (t65 * t88 + (qJD(2) + t144) * t89) * t117) * t114 * t142) / t73, 0, 0, 0, 0; 0, (t84 * t138 - t140) * t97 * t79 - 0.2e1 * (t140 * t145 + (-t84 * t86 * t145 + (-t85 * t138 + t139) * t96 * t79) * t115) * t124, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (941->56), mult. (2271->131), div. (423->14), fcn. (2956->11), ass. (0->66)
	t146 = sin(qJ(2));
	t147 = cos(qJ(2));
	t144 = sin(pkin(11));
	t173 = cos(pkin(6));
	t161 = t144 * t173;
	t172 = cos(pkin(11));
	t132 = -t146 * t161 + t172 * t147;
	t140 = pkin(12) + qJ(4);
	t136 = sin(t140);
	t137 = cos(t140);
	t145 = sin(pkin(6));
	t166 = t144 * t145;
	t155 = -t132 * t136 + t137 * t166;
	t177 = t155 * qJD(4);
	t158 = t173 * t172;
	t128 = t144 * t146 - t147 * t158;
	t165 = t145 * t147;
	t118 = atan2(-t128, -t165);
	t116 = sin(t118);
	t117 = cos(t118);
	t103 = -t116 * t128 - t117 * t165;
	t100 = 0.1e1 / t103;
	t115 = t132 * t137 + t136 * t166;
	t111 = 0.1e1 / t115;
	t141 = 0.1e1 / t147;
	t101 = 0.1e1 / t103 ^ 2;
	t112 = 0.1e1 / t115 ^ 2;
	t142 = 0.1e1 / t147 ^ 2;
	t130 = t144 * t147 + t146 * t158;
	t123 = t130 * qJD(2);
	t164 = qJD(2) * t146;
	t167 = t142 * t146;
	t162 = t128 * t167;
	t126 = t128 ^ 2;
	t139 = 0.1e1 / t145 ^ 2;
	t121 = t126 * t139 * t142 + 0.1e1;
	t119 = 0.1e1 / t121;
	t138 = 0.1e1 / t145;
	t168 = t119 * t138;
	t95 = (qJD(2) * t162 + t123 * t141) * t168;
	t92 = (-t128 * t95 + t145 * t164) * t117 + (t95 * t165 - t123) * t116;
	t176 = t100 * t101 * t92;
	t110 = t155 ^ 2;
	t106 = t110 * t112 + 0.1e1;
	t154 = -t172 * t146 - t147 * t161;
	t124 = t154 * qJD(2);
	t108 = t115 * qJD(4) + t124 * t136;
	t169 = t112 * t155;
	t109 = t124 * t137 + t177;
	t170 = t109 * t111 * t112;
	t175 = 0.1e1 / t106 ^ 2 * (-t108 * t169 - t110 * t170);
	t156 = t130 * t141 + t162;
	t96 = t156 * t168;
	t174 = t128 * t96;
	t171 = t101 * t154;
	t163 = -0.2e1 * t175;
	t157 = -t111 * t136 - t137 * t169;
	t143 = t141 * t142;
	t127 = t154 ^ 2;
	t125 = t132 * qJD(2);
	t122 = qJD(2) * t128;
	t104 = 0.1e1 / t106;
	t99 = t127 * t101 + 0.1e1;
	t93 = (t145 * t146 - t174) * t117 + (t96 * t165 - t130) * t116;
	t91 = (-0.2e1 * t156 / t121 ^ 2 * (t123 * t128 * t142 + t126 * t143 * t164) * t139 + (t123 * t167 - t122 * t141 + (t130 * t167 + (0.2e1 * t143 * t146 ^ 2 + t141) * t128) * qJD(2)) * t119) * t138;
	t1 = [0, t91, 0, 0, 0, 0; 0, 0.2e1 * (-t100 * t132 - t93 * t171) / t99 ^ 2 * (-t125 * t171 - t127 * t176) + (t124 * t100 + (-t93 * t125 - t132 * t92) * t101 - (0.2e1 * t93 * t176 + (-(-t123 * t96 - t128 * t91 - t130 * t95 + (t95 * t96 + qJD(2)) * t165) * t117 - (t95 * t174 + t122 + (t147 * t91 + (-qJD(2) * t96 - t95) * t146) * t145) * t116) * t101) * t154) / t99, 0, 0, 0, 0; 0, -t157 * t154 * t163 + (t157 * t125 - ((-qJD(4) * t111 + 0.2e1 * t155 * t170) * t137 + (t108 * t137 + (t109 + t177) * t136) * t112) * t154) * t104, 0, t163 - 0.2e1 * (t104 * t108 * t112 - (-t104 * t170 - t112 * t175) * t155) * t155, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (1419->57), mult. (2577->133), div. (441->14), fcn. (3313->11), ass. (0->69)
	t173 = sin(qJ(2));
	t174 = cos(qJ(2));
	t171 = sin(pkin(11));
	t204 = cos(pkin(6));
	t189 = t171 * t204;
	t203 = cos(pkin(11));
	t158 = -t173 * t189 + t203 * t174;
	t185 = t204 * t203;
	t154 = t171 * t173 - t174 * t185;
	t172 = sin(pkin(6));
	t193 = t172 * t174;
	t144 = atan2(-t154, -t193);
	t142 = sin(t144);
	t143 = cos(t144);
	t129 = -t142 * t154 - t143 * t193;
	t126 = 0.1e1 / t129;
	t164 = pkin(12) + qJ(4) + qJ(5);
	t162 = sin(t164);
	t163 = cos(t164);
	t194 = t171 * t172;
	t141 = t158 * t163 + t162 * t194;
	t137 = 0.1e1 / t141;
	t168 = 0.1e1 / t174;
	t127 = 0.1e1 / t129 ^ 2;
	t138 = 0.1e1 / t141 ^ 2;
	t169 = 0.1e1 / t174 ^ 2;
	t156 = t171 * t174 + t173 * t185;
	t149 = t156 * qJD(2);
	t195 = t169 * t173;
	t190 = t154 * t195;
	t152 = t154 ^ 2;
	t166 = 0.1e1 / t172 ^ 2;
	t147 = t152 * t166 * t169 + 0.1e1;
	t145 = 0.1e1 / t147;
	t165 = 0.1e1 / t172;
	t197 = t145 * t165;
	t121 = (qJD(2) * t190 + t149 * t168) * t197;
	t183 = t142 * t193 - t143 * t154;
	t191 = t143 * t172 * t173;
	t118 = qJD(2) * t191 + t183 * t121 - t142 * t149;
	t202 = t118 * t126 * t127;
	t140 = t158 * t162 - t163 * t194;
	t136 = t140 ^ 2;
	t132 = t136 * t138 + 0.1e1;
	t181 = -t203 * t173 - t174 * t189;
	t150 = t181 * qJD(2);
	t167 = qJD(4) + qJD(5);
	t186 = t167 * t194 + t150;
	t196 = t158 * t167;
	t133 = t186 * t162 + t163 * t196;
	t198 = t138 * t140;
	t134 = -t162 * t196 + t186 * t163;
	t199 = t134 * t137 * t138;
	t201 = (t133 * t198 - t136 * t199) / t132 ^ 2;
	t200 = t127 * t181;
	t192 = -0.2e1 * t201;
	t184 = -t137 * t162 + t163 * t198;
	t182 = t156 * t168 + t190;
	t170 = t168 * t169;
	t153 = t181 ^ 2;
	t151 = t158 * qJD(2);
	t148 = qJD(2) * t154;
	t130 = 0.1e1 / t132;
	t125 = t153 * t127 + 0.1e1;
	t122 = t182 * t197;
	t119 = t183 * t122 - t142 * t156 + t191;
	t117 = (-0.2e1 * t182 / t147 ^ 2 * (qJD(2) * t152 * t170 * t173 + t149 * t154 * t169) * t166 + (t149 * t195 - t148 * t168 + (t156 * t195 + (0.2e1 * t170 * t173 ^ 2 + t168) * t154) * qJD(2)) * t145) * t165;
	t115 = t192 + 0.2e1 * (t130 * t133 * t138 + (-t130 * t199 - t138 * t201) * t140) * t140;
	t1 = [0, t117, 0, 0, 0, 0; 0, 0.2e1 * (-t119 * t200 - t126 * t158) / t125 ^ 2 * (-t151 * t200 - t153 * t202) + (t150 * t126 + (-t158 * t118 - t119 * t151) * t127 - (0.2e1 * t119 * t202 + (-(qJD(2) * t193 - t117 * t154 - t122 * t149 + (t122 * t193 - t156) * t121) * t143 - (t121 * t122 * t154 + t148 + (t117 * t174 + (-qJD(2) * t122 - t121) * t173) * t172) * t142) * t127) * t181) / t125, 0, 0, 0, 0; 0, -t184 * t181 * t192 + (t184 * t151 - ((-t137 * t167 - 0.2e1 * t140 * t199) * t163 + (t133 * t163 + (-t140 * t167 + t134) * t162) * t138) * t181) * t130, 0, t115, t115, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:37
	% DurationCPUTime: 1.76s
	% Computational Cost: add. (13183->114), mult. (13312->231), div. (822->12), fcn. (17117->13), ass. (0->113)
	t260 = sin(pkin(11));
	t262 = cos(pkin(11));
	t265 = sin(qJ(2));
	t263 = cos(pkin(6));
	t267 = cos(qJ(2));
	t294 = t263 * t267;
	t248 = -t260 * t265 + t262 * t294;
	t244 = t248 * qJD(2);
	t295 = t263 * t265;
	t249 = t260 * t267 + t262 * t295;
	t258 = pkin(12) + qJ(4) + qJ(5);
	t256 = sin(t258);
	t259 = qJD(4) + qJD(5);
	t261 = sin(pkin(6));
	t298 = t261 * t262;
	t283 = t256 * t298;
	t257 = cos(t258);
	t300 = t257 * t259;
	t211 = t244 * t256 + t249 * t300 - t259 * t283;
	t233 = t249 * t256 + t257 * t298;
	t231 = t233 ^ 2;
	t297 = t261 * t265;
	t285 = t256 * t297;
	t241 = -t263 * t257 + t285;
	t239 = 0.1e1 / t241 ^ 2;
	t219 = t231 * t239 + 0.1e1;
	t217 = 0.1e1 / t219;
	t292 = qJD(2) * t267;
	t276 = t259 * t263 + t261 * t292;
	t284 = t257 * t297;
	t229 = t276 * t256 + t259 * t284;
	t238 = 0.1e1 / t241;
	t305 = t233 * t239;
	t195 = (-t211 * t238 + t229 * t305) * t217;
	t220 = atan2(-t233, t241);
	t215 = sin(t220);
	t216 = cos(t220);
	t279 = -t215 * t241 - t216 * t233;
	t191 = t279 * t195 - t215 * t211 + t216 * t229;
	t205 = -t215 * t233 + t216 * t241;
	t202 = 0.1e1 / t205;
	t203 = 0.1e1 / t205 ^ 2;
	t319 = t191 * t202 * t203;
	t286 = t260 * t295;
	t251 = t262 * t267 - t286;
	t299 = t260 * t261;
	t236 = t251 * t256 - t257 * t299;
	t318 = 0.2e1 * t236 * t319;
	t296 = t261 * t267;
	t275 = -t238 * t248 + t296 * t305;
	t317 = t256 * t275;
	t306 = t229 * t238 * t239;
	t316 = -0.2e1 * (t211 * t305 - t231 * t306) / t219 ^ 2;
	t237 = t251 * t257 + t256 * t299;
	t266 = cos(qJ(6));
	t250 = t260 * t294 + t262 * t265;
	t264 = sin(qJ(6));
	t303 = t250 * t264;
	t226 = t237 * t266 + t303;
	t222 = 0.1e1 / t226;
	t223 = 0.1e1 / t226 ^ 2;
	t246 = t250 * qJD(2);
	t281 = t259 * t299 - t246;
	t301 = t256 * t259;
	t214 = -t251 * t301 + t281 * t257;
	t247 = -qJD(2) * t286 + t262 * t292;
	t206 = t226 * qJD(6) + t214 * t264 - t247 * t266;
	t302 = t250 * t266;
	t225 = t237 * t264 - t302;
	t221 = t225 ^ 2;
	t210 = t221 * t223 + 0.1e1;
	t308 = t223 * t225;
	t291 = qJD(6) * t225;
	t207 = t214 * t266 + t247 * t264 - t291;
	t313 = t207 * t222 * t223;
	t315 = (t206 * t308 - t221 * t313) / t210 ^ 2;
	t314 = t203 * t236;
	t213 = t251 * t300 + t281 * t256;
	t312 = t213 * t203;
	t311 = t215 * t236;
	t310 = t216 * t236;
	t309 = t222 * t264;
	t307 = t225 * t266;
	t304 = t250 * t256;
	t293 = qJD(2) * t265;
	t232 = t236 ^ 2;
	t201 = t232 * t203 + 0.1e1;
	t290 = 0.2e1 * (-t232 * t319 + t236 * t312) / t201 ^ 2;
	t289 = -0.2e1 * t315;
	t287 = t225 * t313;
	t282 = -0.2e1 * t233 * t306;
	t280 = qJD(6) * t250 * t257 - t246;
	t278 = t223 * t307 - t309;
	t235 = t249 * t257 - t283;
	t242 = t263 * t256 + t284;
	t277 = -t235 * t238 + t242 * t305;
	t274 = qJD(6) * t251 - t247 * t257 + t250 * t301;
	t245 = t249 * qJD(2);
	t230 = t276 * t257 - t259 * t285;
	t228 = t251 * t264 - t257 * t302;
	t227 = -t251 * t266 - t257 * t303;
	t212 = -t249 * t301 + (-t259 * t298 + t244) * t257;
	t208 = 0.1e1 / t210;
	t199 = 0.1e1 / t201;
	t197 = t217 * t317;
	t196 = t277 * t217;
	t193 = (-t215 * t248 + t216 * t296) * t256 + t279 * t197;
	t192 = t279 * t196 - t215 * t235 + t216 * t242;
	t189 = t277 * t316 + (t242 * t282 - t212 * t238 + (t211 * t242 + t229 * t235 + t230 * t233) * t239) * t217;
	t188 = t316 * t317 + (t275 * t300 + (t282 * t296 + t238 * t245 + (t229 * t248 + (t211 * t267 - t233 * t293) * t261) * t239) * t256) * t217;
	t187 = t278 * t236 * t289 + (t278 * t213 + ((-qJD(6) * t222 - 0.2e1 * t287) * t266 + (t206 * t266 + (t207 - t291) * t264) * t223) * t236) * t208;
	t186 = (t192 * t314 - t202 * t237) * t290 + (t192 * t318 + t214 * t202 + (-t237 * t191 - t192 * t213 - (-t189 * t233 - t196 * t211 + t230 + (-t196 * t241 - t235) * t195) * t310 - (-t189 * t241 - t196 * t229 - t212 + (t196 * t233 - t242) * t195) * t311) * t203) * t199;
	t1 = [0, t188, 0, t189, t189, 0; 0, (t193 * t314 + t202 * t304) * t290 + ((-t247 * t256 - t250 * t300) * t202 + (-t312 + t318) * t193 + (t304 * t191 - (-t188 * t233 - t197 * t211 + (-t256 * t293 + t267 * t300) * t261 + (-t197 * t241 - t248 * t256) * t195) * t310 - (-t248 * t300 - t188 * t241 - t197 * t229 + t245 * t256 + (t197 * t233 - t256 * t296) * t195) * t311) * t203) * t199, 0, t186, t186, 0; 0, 0.2e1 * (-t222 * t227 + t228 * t308) * t315 + (0.2e1 * t228 * t287 - t280 * t222 * t266 + t274 * t309 + (-t280 * t225 * t264 - t228 * t206 - t227 * t207 - t274 * t307) * t223) * t208, 0, t187, t187, t289 + 0.2e1 * (t206 * t223 * t208 + (-t208 * t313 - t223 * t315) * t225) * t225;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end