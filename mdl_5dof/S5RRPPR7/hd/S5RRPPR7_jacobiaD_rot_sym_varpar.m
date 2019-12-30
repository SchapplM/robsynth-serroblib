% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRPPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:23
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPPR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:23:04
	% EndTime: 2019-12-29 18:23:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:23:09
	% EndTime: 2019-12-29 18:23:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:23:04
	% EndTime: 2019-12-29 18:23:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:23:04
	% EndTime: 2019-12-29 18:23:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:23:10
	% EndTime: 2019-12-29 18:23:11
	% DurationCPUTime: 1.14s
	% Computational Cost: add. (1883->70), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->75)
	t95 = sin(qJ(1));
	t145 = 0.2e1 * t95;
	t90 = 0.1e1 / t95;
	t96 = cos(qJ(1));
	t129 = t90 * t96;
	t123 = qJD(2) * t95;
	t125 = qJD(1) * t96;
	t88 = qJ(2) + pkin(8);
	t86 = sin(t88);
	t115 = t86 * t125;
	t81 = t86 ^ 2;
	t87 = cos(t88);
	t84 = 0.1e1 / t87 ^ 2;
	t132 = t81 * t84;
	t89 = t95 ^ 2;
	t76 = t132 * t89 + 0.1e1;
	t74 = 0.1e1 / t76;
	t83 = 0.1e1 / t87;
	t60 = (-(-t123 * t87 - t115) * t83 + t123 * t132) * t74;
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
	t136 = (t103 * t89 + t108 * t84) / t76 ^ 2;
	t135 = t71 * t95;
	t133 = t81 * t83;
	t131 = t81 * t94;
	t91 = 0.1e1 / t95 ^ 2;
	t128 = t91 * t94;
	t126 = qJD(1) * t95;
	t124 = qJD(2) * t87;
	t122 = qJD(2) * t96;
	t107 = t86 * t94 * t124;
	t63 = t131 * t65 + 0.1e1;
	t121 = 0.2e1 * (-t131 * t140 + (t107 - t108) * t65) / t63 ^ 2;
	t120 = 0.2e1 * t140;
	t79 = t128 * t82 + 0.1e1;
	t119 = 0.2e1 * (-t107 * t91 - t143 * t82) / t79 ^ 2;
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
	t59 = (t141 * t71 * t86 - t109 * t72) * t96;
	t58 = -t87 * t135 + t134 + (-t127 * t72 + t71 * t87) * t69;
	t57 = -t114 * t110 + (qJD(1) * t106 + t103 * t145) * t74;
	t1 = [t96 * t83 * t111 + (qJD(2) * t106 - t126 * t130) * t74, t57, 0, 0, 0; (t64 * t112 + (-t64 * t124 + (qJD(1) * t59 + t56) * t138) * t61) * t95 + (t65 * t112 * t59 + (-((t109 * t60 + t124 * t141 + t111) * t71 + (t110 * t133 - t139 + (t139 + (-t80 * t84 + t142) * t123) * t74) * t72) * t118 + (t120 * t86 - t124 * t65) * t59 + (-t64 + ((-t89 + t94) * t72 * t116 + t141 * t117) * t65) * t86 * qJD(1)) * t61) * t96, (t138 * t58 - t64 * t87) * t96 * t121 + ((-t64 * t126 + (-qJD(2) * t58 - t56) * t137) * t87 + (-t64 * t122 - (-t57 * t72 * t95 - t144 * t71 + (-qJD(2) * t71 - t125 * t72 + t135 * t60) * t69) * t118 + (t120 * t96 + t126 * t65) * t58 - ((t57 - t125) * t71 + ((-t69 * t95 + 0.1e1) * qJD(2) + (t69 - t95) * t60) * t72) * t87 * t137) * t86) * t61, 0, 0, 0; t113 * t87 * t119 + (qJD(2) * t105 + 0.2e1 * t143 * t87) * t77, t86 * t119 * t129 + (-t122 * t87 * t90 + qJD(1) * t105) * t77, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:23:10
	% EndTime: 2019-12-29 18:23:11
	% DurationCPUTime: 1.54s
	% Computational Cost: add. (2081->92), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->90)
	t138 = sin(qJ(1));
	t140 = cos(qJ(1));
	t134 = qJ(2) + pkin(8);
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
	t1 = [t127 * t140 * t160 + (qJD(2) * t152 + t178 * t187) * t123, t93, 0, 0, 0; (t103 * t169 + (t162 + (qJD(1) * t96 + t92) * t198) * t133) * t138 + ((t96 * t169 + (t96 * t173 + ((t97 * t154 + t199 * t176 + t160) * t120 + (t159 * t188 + t133 * t97 + (t130 * t165 - (t97 - 0.2e1 * t174) * t133) * t123) * t121) * t98 * t140) * t133) * t104 + (t96 * t172 + (-t103 + ((t135 - t136) * t121 * t167 + t199 * t168) * t104) * qJD(1)) * t133 * t98) * t140, (t103 * t98 * t178 + (t162 + (qJD(2) * t94 + t92) * t198) * t140) * t132 + (((-qJD(2) * t103 + t94 * t172) * t140 + (t94 * t178 + (-(-t108 * t177 - t138 * t93) * t121 - ((t108 * t138 - 0.1e1) * t97 + (-t108 + t138) * qJD(2)) * t120) * t133 * t140) * t104) * t98 + (t94 * t173 - ((-t93 + t177) * t120 + (t157 * t108 - t158) * t121) * t98 * t132) * t104 * t140) * t133, 0, 0, 0; (-t113 * t118 + t119 * t192) * t171 + (t119 * t161 + (-t119 * t101 - t118 * t102 + t156 * t116 * t180 - (-t133 * t174 - t155 * t140) * t191) * t114 + (t155 * t179 + (-t156 * t137 + t139 * t175) * t138) * t113) * t109, t133 * t148 * t171 + (t148 * t176 + (t150 * t178 + ((qJD(5) * t113 + t161) * t137 + (-t101 * t137 + (-qJD(5) * t116 + t102) * t139) * t114) * t140) * t133) * t109, 0, 0, -0.2e1 * t196 + 0.2e1 * (t101 * t109 * t114 + (-t109 * t194 - t114 * t196) * t116) * t116;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end