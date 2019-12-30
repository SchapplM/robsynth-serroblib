% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP10
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
%   Wie in S5RRPRP10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPRP10_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP10_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:18
	% EndTime: 2019-12-29 18:59:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:18
	% EndTime: 2019-12-29 18:59:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:13
	% EndTime: 2019-12-29 18:59:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:18
	% EndTime: 2019-12-29 18:59:19
	% DurationCPUTime: 1.13s
	% Computational Cost: add. (774->70), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->75)
	t87 = cos(qJ(1));
	t116 = qJD(1) * t87;
	t85 = sin(qJ(1));
	t138 = 0.2e1 * t85;
	t74 = t85 ^ 2;
	t75 = 0.1e1 / t85;
	t83 = t87 ^ 2;
	t136 = (0.1e1 + 0.1e1 / t74 * t83) * t75 * t116;
	t84 = sin(qJ(2));
	t118 = t85 * t84;
	t86 = cos(qJ(2));
	t65 = atan2(-t118, -t86);
	t63 = sin(t65);
	t108 = t63 * t118;
	t64 = cos(t65);
	t59 = -t64 * t86 - t108;
	t56 = 0.1e1 / t59;
	t79 = 0.1e1 / t86;
	t57 = 0.1e1 / t59 ^ 2;
	t135 = -0.2e1 * t84;
	t73 = t84 ^ 2;
	t80 = 0.1e1 / t86 ^ 2;
	t123 = t73 * t80;
	t70 = t74 * t123 + 0.1e1;
	t66 = 0.1e1 / t70;
	t134 = t66 - 0.1e1;
	t106 = t84 * t116;
	t115 = qJD(2) * t85;
	t125 = t64 * t84;
	t114 = qJD(2) * t86;
	t52 = (-(-t85 * t114 - t106) * t79 + t115 * t123) * t66;
	t48 = (-t52 * t85 + qJD(2)) * t125 + (-t106 + (t52 - t115) * t86) * t63;
	t133 = t48 * t56 * t57;
	t132 = t52 * t63;
	t131 = t52 * t84;
	t130 = t57 * t84;
	t129 = t57 * t87;
	t120 = t79 * t84;
	t72 = t84 * t73;
	t78 = t86 ^ 2;
	t94 = qJD(2) * (t72 * t79 / t78 + t120);
	t99 = t73 * t85 * t116;
	t128 = (t74 * t94 + t80 * t99) / t70 ^ 2;
	t105 = 0.1e1 + t123;
	t62 = t105 * t85 * t66;
	t127 = t62 * t85;
	t126 = t63 * t86;
	t124 = t73 * t79;
	t122 = t73 * t83;
	t76 = 0.1e1 / t85 ^ 2;
	t121 = t76 * t83;
	t119 = t84 * t87;
	t117 = qJD(1) * t85;
	t113 = qJD(2) * t87;
	t55 = t57 * t122 + 0.1e1;
	t98 = t83 * t84 * t114;
	t112 = 0.2e1 * (-t122 * t133 + (t98 - t99) * t57) / t55 ^ 2;
	t111 = 0.2e1 * t133;
	t71 = t78 * t121 + 0.1e1;
	t110 = 0.2e1 * (-t78 * t136 - t76 * t98) / t71 ^ 2;
	t109 = t57 * t119;
	t107 = t66 * t124;
	t104 = 0.1e1 + t121;
	t103 = t84 * t112;
	t102 = t128 * t135;
	t101 = t128 * t138;
	t100 = t85 * t107;
	t97 = t105 * t87;
	t96 = t104 * t84;
	t68 = 0.1e1 / t71;
	t53 = 0.1e1 / t55;
	t51 = (t134 * t84 * t63 - t64 * t100) * t87;
	t50 = -t85 * t126 + t125 + (-t64 * t118 + t126) * t62;
	t49 = -t105 * t101 + (qJD(1) * t97 + t94 * t138) * t66;
	t1 = [t87 * t79 * t102 + (qJD(2) * t97 - t117 * t120) * t66, t49, 0, 0, 0; (t56 * t103 + (-t56 * t114 + (qJD(1) * t51 + t48) * t130) * t53) * t85 + (t57 * t103 * t51 + (-((t52 * t100 + t134 * t114 + t102) * t63 + (t101 * t124 - t131 + (t131 + (-t72 * t80 + t135) * t115) * t66) * t64) * t109 + (t84 * t111 - t57 * t114) * t51 + (-t56 + ((-t74 + t83) * t64 * t107 + t134 * t108) * t57) * t84 * qJD(1)) * t53) * t87, (t50 * t130 - t56 * t86) * t87 * t112 + ((-t56 * t117 + (-qJD(2) * t50 - t48) * t129) * t86 + (-t56 * t113 - (-t49 * t64 * t85 + t63 * t115 + t127 * t132 - t132 + (-qJD(2) * t63 - t116 * t64) * t62) * t109 + (t87 * t111 + t57 * t117) * t50 - ((t49 - t116) * t63 + ((0.1e1 - t127) * qJD(2) + (t62 - t85) * t52) * t64) * t86 * t129) * t84) * t53, 0, 0, 0; t104 * t86 * t110 + (qJD(2) * t96 + 0.2e1 * t86 * t136) * t68, t75 * t110 * t119 + (-t75 * t86 * t113 + qJD(1) * t96) * t68, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:13
	% EndTime: 2019-12-29 18:59:15
	% DurationCPUTime: 1.46s
	% Computational Cost: add. (813->90), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->93)
	t131 = cos(qJ(2));
	t129 = sin(qJ(1));
	t123 = t129 ^ 2;
	t128 = sin(qJ(2));
	t121 = 0.1e1 / t128 ^ 2;
	t125 = t131 ^ 2;
	t179 = t121 * t125;
	t118 = t123 * t179 + 0.1e1;
	t116 = 0.1e1 / t118;
	t120 = 0.1e1 / t128;
	t169 = qJD(2) * t129;
	t157 = t121 * t169;
	t132 = cos(qJ(1));
	t171 = qJD(1) * t132;
	t158 = t131 * t171;
	t90 = ((t128 * t169 - t158) * t120 + t125 * t157) * t116;
	t193 = t131 * t90;
	t148 = -t90 + t169;
	t146 = qJD(1) * t128 + qJD(4);
	t167 = qJD(2) * t132;
	t192 = t146 * t129 - t131 * t167;
	t176 = t129 * t131;
	t115 = atan2(-t176, t128);
	t114 = cos(t115);
	t113 = sin(t115);
	t161 = t113 * t176;
	t99 = t114 * t128 - t161;
	t96 = 0.1e1 / t99;
	t127 = sin(qJ(4));
	t175 = t132 * t127;
	t130 = cos(qJ(4));
	t177 = t129 * t130;
	t110 = t128 * t175 + t177;
	t106 = 0.1e1 / t110;
	t107 = 0.1e1 / t110 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t191 = t116 - 0.1e1;
	t149 = -t129 * t90 + qJD(2);
	t181 = t114 * t131;
	t85 = t149 * t181 + (t148 * t128 - t158) * t113;
	t190 = t85 * t96 * t97;
	t126 = t132 ^ 2;
	t95 = t126 * t125 * t97 + 0.1e1;
	t93 = 0.1e1 / t95;
	t189 = t93 * t97;
	t188 = t96 * t93;
	t174 = t132 * t130;
	t178 = t129 * t127;
	t109 = -t128 * t174 + t178;
	t105 = t109 ^ 2;
	t104 = t105 * t107 + 0.1e1;
	t184 = t107 * t109;
	t147 = qJD(4) * t128 + qJD(1);
	t143 = t147 * t132;
	t92 = -t192 * t127 + t130 * t143;
	t186 = t106 * t107 * t92;
	t91 = t127 * t143 + t192 * t130;
	t187 = 0.1e1 / t104 ^ 2 * (-t105 * t186 + t91 * t184);
	t185 = t132 * t97;
	t183 = t109 * t127;
	t182 = t113 * t129;
	t180 = t120 * t125;
	t173 = qJD(1) * t129;
	t172 = qJD(1) * t131;
	t170 = qJD(2) * t128;
	t168 = qJD(2) * t131;
	t159 = t129 * t171;
	t162 = t97 * t170;
	t166 = 0.2e1 * (-t126 * t131 * t162 + (-t126 * t190 - t97 * t159) * t125) / t95 ^ 2;
	t165 = 0.2e1 * t190;
	t164 = 0.2e1 * t187;
	t124 = t131 * t125;
	t141 = qJD(2) * (-t121 * t124 - t131) * t120;
	t163 = 0.2e1 * (t123 * t141 + t159 * t179) / t118 ^ 2;
	t160 = t116 * t180;
	t155 = t96 * t166;
	t154 = t97 * t166;
	t153 = 0.1e1 + t179;
	t152 = 0.2e1 * t109 * t186;
	t151 = t129 * t163;
	t150 = t131 * t163;
	t145 = t129 * t160;
	t144 = t153 * t132;
	t142 = t106 * t130 + t107 * t183;
	t140 = t142 * t132;
	t112 = -t128 * t178 + t174;
	t111 = t128 * t177 + t175;
	t103 = t153 * t129 * t116;
	t101 = 0.1e1 / t104;
	t89 = (t191 * t131 * t113 + t114 * t145) * t132;
	t88 = t128 * t182 + t181 + (-t113 * t128 - t114 * t176) * t103;
	t86 = -t153 * t151 + (qJD(1) * t144 + 0.2e1 * t129 * t141) * t116;
	t1 = [t132 * t120 * t150 + (t120 * t129 * t172 + qJD(2) * t144) * t116, t86, 0, 0, 0; (t170 * t188 + (t155 + (qJD(1) * t89 + t85) * t189) * t131) * t129 + (t89 * t154 * t131 + (t89 * t162 + (t89 * t165 + ((t90 * t145 + t191 * t170 + t150) * t113 + (t151 * t180 + t193 + (t124 * t157 - (t90 - 0.2e1 * t169) * t131) * t116) * t114) * t185) * t131 + (-t96 + (-(-t123 + t126) * t114 * t160 + t191 * t161) * t97) * t172) * t93) * t132, (t173 * t188 + (t155 + (qJD(2) * t88 + t85) * t189) * t132) * t128 + (t88 * t132 * t154 + (-t96 * t167 + (t132 * t165 + t97 * t173) * t88 + (-t103 * t182 * t193 + (-(-t103 * t171 - t129 * t86) * t131 - (t148 * t103 - t149) * t128) * t114 + (-(-qJD(2) * t103 + t148) * t131 - (-t86 + t171) * t128) * t113) * t185) * t93) * t131, 0, 0, 0; (-t106 * t111 + t112 * t184) * t164 + (t112 * t152 + (-t111 * t92 - t112 * t91 + t147 * t109 * t177 - (-t129 * t168 - t146 * t132) * t183) * t107 + (t146 * t174 + (-t147 * t127 + t130 * t168) * t129) * t106) * t101, t131 * t140 * t164 + (t140 * t170 + (t142 * t173 + ((qJD(4) * t106 + t152) * t127 + (-t127 * t91 + (-qJD(4) * t109 + t92) * t130) * t107) * t132) * t131) * t101, 0, -0.2e1 * t187 + 0.2e1 * (t91 * t107 * t101 + (-t101 * t186 - t107 * t187) * t109) * t109, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:59:13
	% EndTime: 2019-12-29 18:59:14
	% DurationCPUTime: 1.40s
	% Computational Cost: add. (813->91), mult. (2519->208), div. (480->12), fcn. (2968->9), ass. (0->88)
	t133 = sin(qJ(1));
	t136 = cos(qJ(1));
	t132 = sin(qJ(2));
	t151 = qJD(1) * t132 + qJD(4);
	t135 = cos(qJ(2));
	t170 = qJD(2) * t135;
	t193 = t133 * t151 - t136 * t170;
	t178 = t133 * t135;
	t119 = atan2(-t178, t132);
	t118 = cos(t119);
	t117 = sin(t119);
	t164 = t117 * t178;
	t103 = t118 * t132 - t164;
	t100 = 0.1e1 / t103;
	t131 = sin(qJ(4));
	t134 = cos(qJ(4));
	t179 = t133 * t134;
	t180 = t132 * t136;
	t114 = t131 * t180 + t179;
	t110 = 0.1e1 / t114;
	t124 = 0.1e1 / t132;
	t101 = 0.1e1 / t103 ^ 2;
	t111 = 0.1e1 / t114 ^ 2;
	t125 = 0.1e1 / t132 ^ 2;
	t127 = t133 ^ 2;
	t129 = t135 ^ 2;
	t183 = t125 * t129;
	t122 = t127 * t183 + 0.1e1;
	t120 = 0.1e1 / t122;
	t192 = t120 - 0.1e1;
	t130 = t136 ^ 2;
	t182 = t129 * t130;
	t99 = t101 * t182 + 0.1e1;
	t97 = 0.1e1 / t99;
	t191 = t101 * t97;
	t171 = qJD(2) * t133;
	t161 = t125 * t171;
	t173 = qJD(1) * t136;
	t162 = t135 * t173;
	t94 = ((t132 * t171 - t162) * t124 + t129 * t161) * t120;
	t153 = -t94 + t171;
	t154 = -t133 * t94 + qJD(2);
	t185 = t118 * t135;
	t89 = t154 * t185 + (t132 * t153 - t162) * t117;
	t190 = t100 * t101 * t89;
	t177 = t134 * t136;
	t113 = t131 * t133 - t132 * t177;
	t109 = t113 ^ 2;
	t108 = t109 * t111 + 0.1e1;
	t187 = t111 * t113;
	t152 = qJD(4) * t132 + qJD(1);
	t147 = t152 * t136;
	t96 = -t193 * t131 + t134 * t147;
	t188 = t110 * t111 * t96;
	t95 = t131 * t147 + t193 * t134;
	t189 = 0.1e1 / t108 ^ 2 * (-t109 * t188 + t95 * t187);
	t186 = t113 * t131;
	t184 = t124 * t129;
	t181 = t132 * t133;
	t175 = qJD(1) * t133;
	t174 = qJD(1) * t135;
	t172 = qJD(2) * t132;
	t149 = t129 * t133 * t173;
	t169 = 0.2e1 * (-t182 * t190 + (-t130 * t132 * t170 - t149) * t101) / t99 ^ 2;
	t168 = 0.2e1 * t190;
	t167 = 0.2e1 * t189;
	t128 = t135 * t129;
	t145 = qJD(2) * (-t125 * t128 - t135) * t124;
	t166 = 0.2e1 * (t125 * t149 + t127 * t145) / t122 ^ 2;
	t165 = t97 * t172;
	t163 = t120 * t184;
	t159 = 0.1e1 + t183;
	t158 = t100 * t169;
	t157 = 0.2e1 * t113 * t188;
	t156 = t133 * t166;
	t155 = t135 * t166;
	t150 = t133 * t163;
	t148 = t159 * t136;
	t146 = t110 * t134 + t111 * t186;
	t144 = t146 * t136;
	t116 = -t131 * t181 + t177;
	t115 = t131 * t136 + t132 * t179;
	t107 = t159 * t133 * t120;
	t105 = 0.1e1 / t108;
	t93 = (t192 * t135 * t117 + t118 * t150) * t136;
	t92 = t117 * t181 + t185 + (-t117 * t132 - t118 * t178) * t107;
	t90 = -t159 * t156 + (qJD(1) * t148 + 0.2e1 * t133 * t145) * t120;
	t1 = [t124 * t136 * t155 + (t124 * t133 * t174 + qJD(2) * t148) * t120, t90, 0, 0, 0; (t100 * t165 + (t158 + (qJD(1) * t93 + t89) * t191) * t135) * t133 + (t93 * t135 * t97 * t168 + (t93 * t165 + (t93 * t169 + ((t94 * t150 + t192 * t172 + t155) * t117 + (t156 * t184 + t135 * t94 + (t128 * t161 - (t94 - 0.2e1 * t171) * t135) * t120) * t118) * t97 * t136) * t135) * t101 + (-t100 + ((t127 - t130) * t118 * t163 + t192 * t164) * t101) * t97 * t174) * t136, (t100 * t97 * t175 + (t158 + (qJD(2) * t92 + t89) * t191) * t136) * t132 + (t92 * t136 * t101 * t169 + ((-qJD(2) * t100 + t92 * t168) * t136 + (t92 * t175 + (-(-t107 * t173 - t133 * t90) * t118 - ((t107 * t133 - 0.1e1) * t94 + (-t107 + t133) * qJD(2)) * t117) * t135 * t136) * t101) * t97 - ((-t90 + t173) * t117 + (t107 * t153 - t154) * t118) * t180 * t191) * t135, 0, 0, 0; (-t110 * t115 + t116 * t187) * t167 + (t116 * t157 + (-t115 * t96 - t116 * t95 + t152 * t113 * t179 - (-t133 * t170 - t136 * t151) * t186) * t111 + (t151 * t177 + (-t131 * t152 + t134 * t170) * t133) * t110) * t105, t135 * t144 * t167 + (t144 * t172 + (t146 * t175 + ((qJD(4) * t110 + t157) * t131 + (-t131 * t95 + (-qJD(4) * t113 + t96) * t134) * t111) * t136) * t135) * t105, 0, -0.2e1 * t189 + 0.2e1 * (t105 * t111 * t95 + (-t105 * t188 - t111 * t189) * t113) * t113, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end