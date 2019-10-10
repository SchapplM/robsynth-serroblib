% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR5
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
%   Wie in S6PRPRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
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
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:22
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (610->48), mult. (1820->126), div. (396->14), fcn. (2421->11), ass. (0->64)
	t116 = sin(qJ(2));
	t117 = cos(qJ(2));
	t113 = sin(pkin(10));
	t136 = cos(pkin(6));
	t129 = t113 * t136;
	t135 = cos(pkin(10));
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
	t112 = sin(pkin(11));
	t115 = cos(pkin(11));
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
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:22
	% DurationCPUTime: 0.54s
	% Computational Cost: add. (941->56), mult. (2271->131), div. (423->14), fcn. (2956->11), ass. (0->66)
	t146 = sin(qJ(2));
	t147 = cos(qJ(2));
	t144 = sin(pkin(10));
	t173 = cos(pkin(6));
	t161 = t144 * t173;
	t172 = cos(pkin(10));
	t132 = -t146 * t161 + t172 * t147;
	t140 = pkin(11) + qJ(4);
	t136 = sin(t140);
	t137 = cos(t140);
	t145 = sin(pkin(6));
	t166 = t144 * t145;
	t155 = -t132 * t136 + t137 * t166;
	t177 = qJD(4) * t155;
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
	t108 = qJD(4) * t115 + t124 * t136;
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
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:23
	% DurationCPUTime: 1.17s
	% Computational Cost: add. (4935->90), mult. (7372->193), div. (524->12), fcn. (9540->11), ass. (0->87)
	t173 = sin(pkin(10));
	t175 = cos(pkin(10));
	t177 = sin(qJ(2));
	t176 = cos(pkin(6));
	t178 = cos(qJ(2));
	t199 = t176 * t178;
	t162 = -t173 * t177 + t175 * t199;
	t155 = t162 * qJD(2);
	t200 = t176 * t177;
	t163 = t173 * t178 + t175 * t200;
	t172 = pkin(11) + qJ(4);
	t170 = sin(t172);
	t174 = sin(pkin(6));
	t203 = t174 * t175;
	t191 = t170 * t203;
	t171 = cos(t172);
	t196 = qJD(4) * t171;
	t127 = -qJD(4) * t191 + t155 * t170 + t163 * t196;
	t145 = t163 * t170 + t171 * t203;
	t142 = t145 ^ 2;
	t202 = t174 * t177;
	t153 = t170 * t202 - t171 * t176;
	t151 = 0.1e1 / t153 ^ 2;
	t135 = t142 * t151 + 0.1e1;
	t133 = 0.1e1 / t135;
	t154 = t170 * t176 + t171 * t202;
	t197 = qJD(2) * t178;
	t190 = t174 * t197;
	t140 = t154 * qJD(4) + t170 * t190;
	t150 = 0.1e1 / t153;
	t207 = t145 * t151;
	t115 = (-t127 * t150 + t140 * t207) * t133;
	t136 = atan2(-t145, t153);
	t131 = sin(t136);
	t132 = cos(t136);
	t188 = -t131 * t153 - t132 * t145;
	t112 = t188 * t115 - t131 * t127 + t132 * t140;
	t126 = -t131 * t145 + t132 * t153;
	t123 = 0.1e1 / t126;
	t124 = 0.1e1 / t126 ^ 2;
	t215 = t112 * t123 * t124;
	t192 = t173 * t200;
	t165 = t175 * t178 - t192;
	t204 = t173 * t174;
	t186 = -t165 * t170 + t171 * t204;
	t214 = -0.2e1 * t186 * t215;
	t201 = t174 * t178;
	t185 = -t150 * t162 + t201 * t207;
	t213 = t170 * t185;
	t208 = t140 * t150 * t151;
	t212 = -0.2e1 * (t127 * t207 - t142 * t208) / t135 ^ 2;
	t164 = t173 * t199 + t175 * t177;
	t159 = 0.1e1 / t164;
	t160 = 0.1e1 / t164 ^ 2;
	t211 = t124 * t186;
	t210 = t131 * t186;
	t209 = t132 * t186;
	t149 = t165 * t171 + t170 * t204;
	t206 = t149 * t165;
	t205 = t164 * t170;
	t198 = qJD(2) * t177;
	t143 = t186 ^ 2;
	t121 = t124 * t143 + 0.1e1;
	t157 = t164 * qJD(2);
	t129 = t149 * qJD(4) - t157 * t170;
	t195 = 0.2e1 * (-t129 * t211 - t143 * t215) / t121 ^ 2;
	t130 = t186 * qJD(4) - t157 * t171;
	t144 = t149 ^ 2;
	t139 = t144 * t160 + 0.1e1;
	t158 = -qJD(2) * t192 + t175 * t197;
	t161 = t159 * t160;
	t194 = 0.2e1 * (t130 * t149 * t160 - t144 * t158 * t161) / t139 ^ 2;
	t189 = -0.2e1 * t145 * t208;
	t147 = t163 * t171 - t191;
	t187 = -t147 * t150 + t154 * t207;
	t156 = t163 * qJD(2);
	t141 = -t153 * qJD(4) + t171 * t190;
	t137 = 0.1e1 / t139;
	t128 = -t145 * qJD(4) + t155 * t171;
	t118 = 0.1e1 / t121;
	t117 = t133 * t213;
	t116 = t187 * t133;
	t114 = (-t131 * t162 + t132 * t201) * t170 + t188 * t117;
	t113 = t188 * t116 - t131 * t147 + t132 * t154;
	t111 = t187 * t212 + (t154 * t189 - t128 * t150 + (t127 * t154 + t140 * t147 + t141 * t145) * t151) * t133;
	t109 = t212 * t213 + (t185 * t196 + (t189 * t201 + t150 * t156 + (t140 * t162 + (t127 * t178 - t145 * t198) * t174) * t151) * t170) * t133;
	t1 = [0, t109, 0, t111, 0, 0; 0, (-t114 * t211 + t123 * t205) * t195 + ((-t158 * t170 - t164 * t196) * t123 + t114 * t214 + (-t114 * t129 + t205 * t112 + (-t109 * t145 - t117 * t127 + (-t170 * t198 + t178 * t196) * t174 + (-t117 * t153 - t162 * t170) * t115) * t209 + (-t162 * t196 - t109 * t153 - t117 * t140 + t156 * t170 + (t117 * t145 - t170 * t201) * t115) * t210) * t124) * t118, 0, (-t113 * t211 - t123 * t149) * t195 + (t113 * t214 + t130 * t123 + (-t149 * t112 - t113 * t129 + (-t111 * t145 - t116 * t127 + t141 + (-t116 * t153 - t147) * t115) * t209 + (-t111 * t153 - t116 * t140 - t128 + (t116 * t145 - t154) * t115) * t210) * t124) * t118, 0, 0; 0, (t159 * t164 * t171 + t160 * t206) * t194 + (qJD(4) * t159 * t205 + (-t130 * t165 + t149 * t157) * t160 + (0.2e1 * t161 * t206 + (t160 * t164 - t159) * t171) * t158) * t137, 0, -t186 * t159 * t194 + (-t158 * t160 * t186 - t129 * t159) * t137, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:22
	% EndTime: 2019-10-09 21:37:23
	% DurationCPUTime: 1.52s
	% Computational Cost: add. (5804->108), mult. (9085->220), div. (559->12), fcn. (11668->13), ass. (0->104)
	t224 = cos(pkin(6));
	t228 = cos(qJ(2));
	t277 = cos(pkin(10));
	t249 = t277 * t228;
	t222 = sin(pkin(10));
	t226 = sin(qJ(2));
	t262 = t222 * t226;
	t212 = t224 * t249 - t262;
	t208 = t212 * qJD(2);
	t221 = pkin(11) + qJ(4);
	t220 = cos(t221);
	t219 = sin(t221);
	t250 = t277 * t226;
	t261 = t222 * t228;
	t237 = -t224 * t250 - t261;
	t223 = sin(pkin(6));
	t251 = t223 * t277;
	t279 = t237 * t219 - t220 * t251;
	t175 = t279 * qJD(4) + t208 * t220;
	t198 = -t219 * t251 - t237 * t220;
	t195 = t198 ^ 2;
	t260 = t223 * t226;
	t207 = t219 * t224 + t220 * t260;
	t204 = 0.1e1 / t207 ^ 2;
	t188 = t195 * t204 + 0.1e1;
	t186 = 0.1e1 / t188;
	t206 = -t219 * t260 + t220 * t224;
	t259 = t223 * t228;
	t252 = qJD(2) * t259;
	t194 = t206 * qJD(4) + t220 * t252;
	t203 = 0.1e1 / t207;
	t267 = t198 * t204;
	t158 = (-t175 * t203 + t194 * t267) * t186;
	t189 = atan2(-t198, t207);
	t182 = sin(t189);
	t183 = cos(t189);
	t244 = -t182 * t207 - t183 * t198;
	t154 = t244 * t158 - t182 * t175 + t183 * t194;
	t168 = -t182 * t198 + t183 * t207;
	t165 = 0.1e1 / t168;
	t166 = 0.1e1 / t168 ^ 2;
	t283 = t154 * t165 * t166;
	t227 = cos(qJ(6));
	t214 = -t224 * t262 + t249;
	t263 = t222 * t223;
	t240 = -t214 * t219 + t220 * t263;
	t225 = sin(qJ(6));
	t238 = -t224 * t261 - t250;
	t265 = t238 * t225;
	t243 = -t227 * t240 + t265;
	t282 = t243 * qJD(6);
	t201 = t214 * t220 + t219 * t263;
	t281 = 0.2e1 * t201 * t283;
	t239 = -t203 * t212 + t259 * t267;
	t280 = t220 * t239;
	t268 = t194 * t203 * t204;
	t278 = -0.2e1 * (t175 * t267 - t195 * t268) / t188 ^ 2;
	t264 = t238 * t227;
	t185 = -t225 * t240 - t264;
	t179 = 0.1e1 / t185;
	t180 = 0.1e1 / t185 ^ 2;
	t276 = t166 * t201;
	t210 = t238 * qJD(2);
	t176 = t201 * qJD(4) + t210 * t219;
	t211 = t214 * qJD(2);
	t170 = t176 * t225 + t211 * t227 + t282;
	t275 = t170 * t179 * t180;
	t169 = t185 * qJD(6) - t176 * t227 + t211 * t225;
	t178 = t243 ^ 2;
	t173 = t178 * t180 + 0.1e1;
	t272 = t180 * t243;
	t274 = 0.1e1 / t173 ^ 2 * (-t169 * t272 - t178 * t275);
	t273 = t179 * t227;
	t271 = t182 * t201;
	t270 = t183 * t201;
	t269 = t243 * t225;
	t266 = t238 * t220;
	t258 = qJD(2) * t226;
	t257 = qJD(4) * t219;
	t196 = t201 ^ 2;
	t164 = t166 * t196 + 0.1e1;
	t177 = t240 * qJD(4) + t210 * t220;
	t256 = 0.2e1 * (t177 * t276 - t196 * t283) / t164 ^ 2;
	t254 = 0.2e1 * t274;
	t248 = -0.2e1 * t243 * t275;
	t247 = -0.2e1 * t198 * t268;
	t245 = qJD(6) * t219 * t238 + t210;
	t242 = -t180 * t269 + t273;
	t241 = -t203 * t279 + t206 * t267;
	t236 = -qJD(4) * t266 + qJD(6) * t214 + t211 * t219;
	t209 = t237 * qJD(2);
	t193 = -t207 * qJD(4) - t219 * t252;
	t191 = t214 * t227 + t219 * t265;
	t190 = t214 * t225 - t219 * t264;
	t174 = t198 * qJD(4) + t208 * t219;
	t171 = 0.1e1 / t173;
	t161 = 0.1e1 / t164;
	t160 = t186 * t280;
	t159 = t241 * t186;
	t156 = (-t182 * t212 + t183 * t259) * t220 + t244 * t160;
	t155 = t244 * t159 - t182 * t279 + t183 * t206;
	t153 = t241 * t278 + (t206 * t247 + t174 * t203 + (t175 * t206 + t193 * t198 + t194 * t279) * t204) * t186;
	t151 = t278 * t280 + (-t239 * t257 + (t247 * t259 - t203 * t209 + (t194 * t212 + (t175 * t228 - t198 * t258) * t223) * t204) * t220) * t186;
	t1 = [0, t151, 0, t153, 0, 0; 0, (t156 * t276 - t165 * t266) * t256 + ((-t211 * t220 - t238 * t257) * t165 + t156 * t281 + (-t156 * t177 - t266 * t154 - (-t151 * t198 - t160 * t175 + (-t220 * t258 - t228 * t257) * t223 + (-t160 * t207 - t212 * t220) * t158) * t270 - (t212 * t257 - t151 * t207 - t160 * t194 - t209 * t220 + (t160 * t198 - t220 * t259) * t158) * t271) * t166) * t161, 0, (t155 * t276 - t165 * t240) * t256 + (t155 * t281 - t176 * t165 + (-t240 * t154 - t155 * t177 - (-t153 * t198 - t159 * t175 + t193 + (-t159 * t207 - t279) * t158) * t270 - (-t153 * t207 - t159 * t194 + t174 + (t159 * t198 - t206) * t158) * t271) * t166) * t161, 0, 0; 0, (-t179 * t190 - t191 * t272) * t254 + (t191 * t248 + t245 * t179 * t225 + t236 * t273 + (t227 * t243 * t245 - t191 * t169 - t190 * t170 - t236 * t269) * t180) * t171, 0, t242 * t201 * t254 + (-t242 * t177 + ((qJD(6) * t179 + t248) * t225 + (-t169 * t225 + (t170 + t282) * t227) * t180) * t201) * t171, 0, -0.2e1 * t274 - 0.2e1 * (t169 * t180 * t171 - (-t171 * t275 - t180 * t274) * t243) * t243;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end