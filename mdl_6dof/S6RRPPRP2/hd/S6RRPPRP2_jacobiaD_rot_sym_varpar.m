% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRP2
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
%   Wie in S6RRPPRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:28
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:46
	% DurationCPUTime: 0.77s
	% Computational Cost: add. (1883->70), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->75)
	t95 = sin(qJ(1));
	t145 = 0.2e1 * t95;
	t90 = 0.1e1 / t95;
	t96 = cos(qJ(1));
	t129 = t90 * t96;
	t123 = qJD(2) * t95;
	t125 = qJD(1) * t96;
	t88 = qJ(2) + pkin(9);
	t86 = sin(t88);
	t115 = t86 * t125;
	t81 = t86 ^ 2;
	t87 = cos(t88);
	t84 = 0.1e1 / t87 ^ 2;
	t132 = t81 * t84;
	t89 = t95 ^ 2;
	t76 = t89 * t132 + 0.1e1;
	t74 = 0.1e1 / t76;
	t83 = 0.1e1 / t87;
	t60 = (-(-t87 * t123 - t115) * t83 + t123 * t132) * t74;
	t144 = t60 - t123;
	t94 = t96 ^ 2;
	t143 = qJD(1) * (0.1e1 / t89 * t94 + 0.1e1) * t129;
	t127 = t95 * t86;
	t73 = atan2(-t127, -t87);
	t71 = sin(t73);
	t117 = t71 * t127;
	t72 = cos(t73);
	t67 = -t72 * t87 - t117;
	t64 = 0.1e1 / t67;
	t65 = 0.1e1 / t67 ^ 2;
	t142 = -0.2e1 * t86;
	t141 = t74 - 0.1e1;
	t134 = t72 * t86;
	t56 = (-t60 * t95 + qJD(2)) * t134 + (t144 * t87 - t115) * t71;
	t140 = t56 * t64 * t65;
	t139 = t60 * t86;
	t138 = t65 * t86;
	t137 = t65 * t96;
	t130 = t83 * t86;
	t80 = t86 * t81;
	t82 = t87 ^ 2;
	t103 = qJD(2) * (t80 * t83 / t82 + t130);
	t108 = t81 * t95 * t125;
	t136 = (t89 * t103 + t84 * t108) / t76 ^ 2;
	t135 = t71 * t95;
	t133 = t81 * t83;
	t131 = t81 * t94;
	t91 = 0.1e1 / t95 ^ 2;
	t128 = t91 * t94;
	t126 = qJD(1) * t95;
	t124 = qJD(2) * t87;
	t122 = qJD(2) * t96;
	t107 = t86 * t94 * t124;
	t63 = t65 * t131 + 0.1e1;
	t121 = 0.2e1 * (-t131 * t140 + (t107 - t108) * t65) / t63 ^ 2;
	t120 = 0.2e1 * t140;
	t79 = t82 * t128 + 0.1e1;
	t119 = 0.2e1 * (-t91 * t107 - t82 * t143) / t79 ^ 2;
	t118 = t86 * t137;
	t116 = t74 * t133;
	t114 = 0.1e1 + t132;
	t113 = 0.1e1 + t128;
	t112 = t86 * t121;
	t111 = t136 * t142;
	t110 = t136 * t145;
	t109 = t95 * t116;
	t106 = t114 * t96;
	t105 = t113 * t86;
	t77 = 0.1e1 / t79;
	t69 = t114 * t95 * t74;
	t61 = 0.1e1 / t63;
	t59 = (t141 * t86 * t71 - t72 * t109) * t96;
	t58 = -t87 * t135 + t134 + (-t72 * t127 + t71 * t87) * t69;
	t57 = -t114 * t110 + (qJD(1) * t106 + t103 * t145) * t74;
	t1 = [t96 * t83 * t111 + (qJD(2) * t106 - t126 * t130) * t74, t57, 0, 0, 0, 0; (t64 * t112 + (-t64 * t124 + (qJD(1) * t59 + t56) * t138) * t61) * t95 + (t65 * t112 * t59 + (-((t60 * t109 + t141 * t124 + t111) * t71 + (t110 * t133 - t139 + (t139 + (-t80 * t84 + t142) * t123) * t74) * t72) * t118 + (t86 * t120 - t65 * t124) * t59 + (-t64 + ((-t89 + t94) * t72 * t116 + t141 * t117) * t65) * t86 * qJD(1)) * t61) * t96, (t58 * t138 - t64 * t87) * t96 * t121 + ((-t64 * t126 + (-qJD(2) * t58 - t56) * t137) * t87 + (-t64 * t122 - (-t57 * t72 * t95 - t144 * t71 + (-qJD(2) * t71 - t125 * t72 + t135 * t60) * t69) * t118 + (t96 * t120 + t65 * t126) * t58 - ((t57 - t125) * t71 + ((-t69 * t95 + 0.1e1) * qJD(2) + (t69 - t95) * t60) * t72) * t87 * t137) * t86) * t61, 0, 0, 0, 0; t113 * t87 * t119 + (qJD(2) * t105 + 0.2e1 * t87 * t143) * t77, t86 * t119 * t129 + (-t87 * t90 * t122 + qJD(1) * t105) * t77, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:46
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (2081->92), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->90)
	t138 = sin(qJ(1));
	t140 = cos(qJ(1));
	t134 = qJ(2) + pkin(9);
	t132 = sin(t134);
	t155 = qJD(1) * t132 + qJD(5);
	t133 = cos(t134);
	t175 = qJD(2) * t133;
	t200 = t155 * t138 - t140 * t175;
	t181 = t138 * t133;
	t122 = atan2(-t181, t132);
	t121 = cos(t122);
	t120 = sin(t122);
	t168 = t120 * t181;
	t106 = t121 * t132 - t168;
	t103 = 0.1e1 / t106;
	t139 = cos(qJ(5));
	t180 = t138 * t139;
	t137 = sin(qJ(5));
	t182 = t137 * t140;
	t117 = t132 * t182 + t180;
	t113 = 0.1e1 / t117;
	t127 = 0.1e1 / t132;
	t104 = 0.1e1 / t106 ^ 2;
	t114 = 0.1e1 / t117 ^ 2;
	t128 = 0.1e1 / t132 ^ 2;
	t135 = t138 ^ 2;
	t131 = t133 ^ 2;
	t186 = t128 * t131;
	t125 = t135 * t186 + 0.1e1;
	t123 = 0.1e1 / t125;
	t199 = t123 - 0.1e1;
	t136 = t140 ^ 2;
	t185 = t131 * t136;
	t100 = t104 * t185 + 0.1e1;
	t98 = 0.1e1 / t100;
	t198 = t104 * t98;
	t174 = qJD(2) * t138;
	t165 = t128 * t174;
	t177 = qJD(1) * t140;
	t166 = t133 * t177;
	t97 = ((t132 * t174 - t166) * t127 + t131 * t165) * t123;
	t157 = -t97 + t174;
	t158 = -t138 * t97 + qJD(2);
	t189 = t121 * t133;
	t92 = t158 * t189 + (t157 * t132 - t166) * t120;
	t197 = t103 * t104 * t92;
	t156 = qJD(5) * t132 + qJD(1);
	t151 = t156 * t140;
	t101 = t137 * t151 + t200 * t139;
	t179 = t139 * t140;
	t183 = t137 * t138;
	t116 = -t132 * t179 + t183;
	t112 = t116 ^ 2;
	t111 = t112 * t114 + 0.1e1;
	t192 = t114 * t116;
	t102 = -t200 * t137 + t139 * t151;
	t194 = t102 * t113 * t114;
	t196 = 0.1e1 / t111 ^ 2 * (t101 * t192 - t112 * t194);
	t191 = t116 * t137;
	t190 = t120 * t132;
	t188 = t127 * t131;
	t187 = t127 * t133;
	t178 = qJD(1) * t138;
	t176 = qJD(2) * t132;
	t153 = t131 * t138 * t177;
	t173 = 0.2e1 * (-t185 * t197 + (-t132 * t136 * t175 - t153) * t104) / t100 ^ 2;
	t172 = 0.2e1 * t197;
	t171 = 0.2e1 * t196;
	t130 = t133 * t131;
	t149 = qJD(2) * (-t127 * t128 * t130 - t187);
	t170 = 0.2e1 * (t128 * t153 + t135 * t149) / t125 ^ 2;
	t169 = t98 * t176;
	t167 = t123 * t188;
	t163 = 0.1e1 + t186;
	t162 = t103 * t173;
	t161 = 0.2e1 * t116 * t194;
	t160 = t133 * t170;
	t159 = t138 * t170;
	t154 = t138 * t167;
	t152 = t163 * t140;
	t150 = t113 * t139 + t114 * t191;
	t148 = t150 * t140;
	t119 = -t132 * t183 + t179;
	t118 = t132 * t180 + t182;
	t109 = 0.1e1 / t111;
	t108 = t163 * t138 * t123;
	t96 = (t199 * t133 * t120 + t121 * t154) * t140;
	t94 = t138 * t190 + t189 + (-t121 * t181 - t190) * t108;
	t93 = -t163 * t159 + (qJD(1) * t152 + 0.2e1 * t138 * t149) * t123;
	t1 = [t127 * t140 * t160 + (qJD(2) * t152 + t178 * t187) * t123, t93, 0, 0, 0, 0; (t103 * t169 + (t162 + (qJD(1) * t96 + t92) * t198) * t133) * t138 + ((t96 * t169 + (t96 * t173 + ((t97 * t154 + t199 * t176 + t160) * t120 + (t159 * t188 + t133 * t97 + (t130 * t165 - (t97 - 0.2e1 * t174) * t133) * t123) * t121) * t98 * t140) * t133) * t104 + (t96 * t172 + (-t103 + ((t135 - t136) * t121 * t167 + t199 * t168) * t104) * qJD(1)) * t133 * t98) * t140, (t103 * t98 * t178 + (t162 + (qJD(2) * t94 + t92) * t198) * t140) * t132 + (((-qJD(2) * t103 + t94 * t172) * t140 + (t94 * t178 + (-(-t108 * t177 - t138 * t93) * t121 - ((t108 * t138 - 0.1e1) * t97 + (-t108 + t138) * qJD(2)) * t120) * t133 * t140) * t104) * t98 + (t94 * t173 - ((-t93 + t177) * t120 + (t157 * t108 - t158) * t121) * t98 * t132) * t104 * t140) * t133, 0, 0, 0, 0; (-t113 * t118 + t119 * t192) * t171 + (t119 * t161 + (-t119 * t101 - t118 * t102 + t156 * t116 * t180 - (-t133 * t174 - t155 * t140) * t191) * t114 + (t155 * t179 + (-t156 * t137 + t139 * t175) * t138) * t113) * t109, t133 * t148 * t171 + (t148 * t176 + (t150 * t178 + ((qJD(5) * t113 + t161) * t137 + (-t101 * t137 + (-qJD(5) * t116 + t102) * t139) * t114) * t140) * t133) * t109, 0, 0, -0.2e1 * t196 + 0.2e1 * (t101 * t109 * t114 + (-t109 * t194 - t114 * t196) * t116) * t116, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:28:45
	% EndTime: 2019-10-10 09:28:46
	% DurationCPUTime: 1.03s
	% Computational Cost: add. (2081->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->92)
	t138 = qJ(2) + pkin(9);
	t136 = sin(t138);
	t132 = 0.1e1 / t136 ^ 2;
	t137 = cos(t138);
	t135 = t137 ^ 2;
	t188 = t132 * t135;
	t204 = t137 * t188;
	t142 = sin(qJ(1));
	t163 = 0.1e1 + t188;
	t203 = t142 * t163;
	t144 = cos(qJ(1));
	t159 = qJD(1) * t136 + qJD(5);
	t175 = qJD(2) * t137;
	t202 = t159 * t142 - t144 * t175;
	t182 = t142 * t137;
	t126 = atan2(-t182, t136);
	t125 = cos(t126);
	t124 = sin(t126);
	t168 = t124 * t182;
	t110 = t125 * t136 - t168;
	t107 = 0.1e1 / t110;
	t143 = cos(qJ(5));
	t181 = t142 * t143;
	t141 = sin(qJ(5));
	t183 = t141 * t144;
	t121 = t136 * t183 + t181;
	t117 = 0.1e1 / t121;
	t131 = 0.1e1 / t136;
	t108 = 0.1e1 / t110 ^ 2;
	t118 = 0.1e1 / t121 ^ 2;
	t139 = t142 ^ 2;
	t129 = t139 * t188 + 0.1e1;
	t127 = 0.1e1 / t129;
	t201 = t127 - 0.1e1;
	t177 = qJD(1) * t144;
	t166 = t137 * t177;
	t174 = qJD(2) * t142;
	t101 = ((t136 * t174 - t166) * t131 + t174 * t188) * t127;
	t190 = t125 * t137;
	t96 = (-t101 * t142 + qJD(2)) * t190 + (-t166 + (-t101 + t174) * t136) * t124;
	t200 = t107 * t108 * t96;
	t160 = qJD(5) * t136 + qJD(1);
	t155 = t160 * t144;
	t105 = t141 * t155 + t202 * t143;
	t180 = t143 * t144;
	t184 = t141 * t142;
	t120 = -t136 * t180 + t184;
	t116 = t120 ^ 2;
	t115 = t116 * t118 + 0.1e1;
	t193 = t118 * t120;
	t106 = -t202 * t141 + t143 * t155;
	t197 = t106 * t117 * t118;
	t199 = 0.1e1 / t115 ^ 2 * (t105 * t193 - t116 * t197);
	t198 = t101 * t137;
	t196 = t108 * t137;
	t195 = t108 * t144;
	t189 = t131 * t137;
	t153 = qJD(2) * (-t131 * t204 - t189);
	t186 = t135 * t142;
	t157 = t177 * t186;
	t194 = (t132 * t157 + t139 * t153) / t129 ^ 2;
	t192 = t120 * t141;
	t191 = t124 * t136;
	t140 = t144 ^ 2;
	t187 = t135 * t140;
	t185 = t137 * t144;
	t112 = t127 * t203;
	t179 = -t112 + t142;
	t178 = qJD(1) * t142;
	t176 = qJD(2) * t136;
	t104 = t108 * t187 + 0.1e1;
	t173 = 0.2e1 / t104 ^ 2 * (-t187 * t200 + (-t136 * t140 * t175 - t157) * t108);
	t172 = 0.2e1 * t200;
	t171 = 0.2e1 * t199;
	t170 = -0.2e1 * t194;
	t169 = t137 * t194;
	t167 = t131 * t186;
	t164 = t112 * t142 - 0.1e1;
	t162 = t137 * t173;
	t161 = 0.2e1 * t120 * t197;
	t158 = t127 * t167;
	t156 = t163 * t144;
	t154 = t117 * t143 + t118 * t192;
	t152 = t154 * t144;
	t123 = -t136 * t184 + t180;
	t122 = t136 * t181 + t183;
	t113 = 0.1e1 / t115;
	t102 = 0.1e1 / t104;
	t100 = (t201 * t137 * t124 + t125 * t158) * t144;
	t98 = t142 * t191 + t190 + (-t125 * t182 - t191) * t112;
	t97 = t170 * t203 + (qJD(1) * t156 + 0.2e1 * t142 * t153) * t127;
	t1 = [0.2e1 * t131 * t144 * t169 + (qJD(2) * t156 + t178 * t189) * t127, t97, 0, 0, 0, 0; (t107 * t162 + (t107 * t176 + (qJD(1) * t100 + t96) * t196) * t102) * t142 + (t108 * t162 * t100 + (-((-t101 * t158 - t201 * t176 - 0.2e1 * t169) * t124 + (t167 * t170 - t198 + (t198 + (-0.2e1 * t137 - t204) * t174) * t127) * t125) * t108 * t185 + (t108 * t176 + t137 * t172) * t100 + (-t107 + ((t139 - t140) * t135 * t131 * t127 * t125 + t201 * t168) * t108) * t137 * qJD(1)) * t102) * t144, (t107 * t136 + t98 * t196) * t144 * t173 + ((t107 * t178 + (qJD(2) * t98 + t96) * t195) * t136 + ((-qJD(2) * t107 + t98 * t172) * t144 + (t98 * t178 + (-(-t112 * t177 - t142 * t97) * t125 - (t179 * qJD(2) + t164 * t101) * t124) * t185) * t108 - ((-t97 + t177) * t124 + (t164 * qJD(2) + t179 * t101) * t125) * t136 * t195) * t137) * t102, 0, 0, 0, 0; (-t117 * t122 + t123 * t193) * t171 + (t123 * t161 + (-t123 * t105 - t122 * t106 + t160 * t120 * t181 - (-t137 * t174 - t159 * t144) * t192) * t118 + (t159 * t180 + (-t160 * t141 + t143 * t175) * t142) * t117) * t113, t137 * t152 * t171 + (t152 * t176 + (t154 * t178 + ((qJD(5) * t117 + t161) * t141 + (-t105 * t141 + (-qJD(5) * t120 + t106) * t143) * t118) * t144) * t137) * t113, 0, 0, -0.2e1 * t199 + 0.2e1 * (t105 * t113 * t118 + (-t113 * t197 - t118 * t199) * t120) * t120, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end