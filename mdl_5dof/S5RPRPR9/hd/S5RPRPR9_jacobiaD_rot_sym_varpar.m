% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPR9
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
%   Wie in S5RPRPR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRPR9_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR9_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:54:36
	% EndTime: 2019-12-29 16:54:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:54:36
	% EndTime: 2019-12-29 16:54:36
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:54:41
	% EndTime: 2019-12-29 16:54:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:54:41
	% EndTime: 2019-12-29 16:54:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:54:36
	% EndTime: 2019-12-29 16:54:38
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (1452->71), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->76)
	t85 = qJ(1) + pkin(8);
	t83 = sin(t85);
	t143 = 0.2e1 * t83;
	t78 = 0.1e1 / t83;
	t84 = cos(t85);
	t130 = t78 * t84;
	t77 = t83 ^ 2;
	t82 = t84 ^ 2;
	t142 = qJD(1) * (0.1e1 / t77 * t82 + 0.1e1) * t130;
	t92 = sin(qJ(3));
	t127 = t83 * t92;
	t93 = cos(qJ(3));
	t73 = atan2(-t127, -t93);
	t68 = sin(t73);
	t114 = t68 * t127;
	t69 = cos(t73);
	t65 = -t69 * t93 - t114;
	t62 = 0.1e1 / t65;
	t89 = 0.1e1 / t93;
	t63 = 0.1e1 / t65 ^ 2;
	t141 = -0.2e1 * t92;
	t87 = t92 ^ 2;
	t90 = 0.1e1 / t93 ^ 2;
	t124 = t87 * t90;
	t76 = t77 * t124 + 0.1e1;
	t74 = 0.1e1 / t76;
	t140 = t74 - 0.1e1;
	t121 = qJD(1) * t92;
	t112 = t84 * t121;
	t120 = qJD(3) * t83;
	t131 = t69 * t92;
	t119 = qJD(3) * t93;
	t57 = (-(-t83 * t119 - t112) * t89 + t120 * t124) * t74;
	t53 = (-t57 * t83 + qJD(3)) * t131 + (-t112 + (t57 - t120) * t93) * t68;
	t139 = t53 * t62 * t63;
	t138 = t57 * t68;
	t137 = t57 * t92;
	t136 = t63 * t84;
	t135 = t63 * t92;
	t86 = t92 * t87;
	t88 = t93 ^ 2;
	t100 = qJD(3) * (t86 / t88 + t92) * t89;
	t122 = qJD(1) * t84;
	t105 = t83 * t87 * t122;
	t134 = (t77 * t100 + t90 * t105) / t76 ^ 2;
	t110 = 0.1e1 + t124;
	t67 = t110 * t83 * t74;
	t133 = t67 * t83;
	t132 = t68 * t93;
	t79 = 0.1e1 / t83 ^ 2;
	t129 = t79 * t82;
	t128 = t82 * t87;
	t126 = t84 * t92;
	t125 = t87 * t89;
	t123 = qJD(1) * t83;
	t104 = t82 * t92 * t119;
	t60 = t63 * t128 + 0.1e1;
	t118 = 0.2e1 * (-t128 * t139 + (t104 - t105) * t63) / t60 ^ 2;
	t117 = 0.2e1 * t139;
	t72 = t88 * t129 + 0.1e1;
	t116 = 0.2e1 * (-t79 * t104 - t88 * t142) / t72 ^ 2;
	t115 = t63 * t126;
	t113 = t74 * t125;
	t111 = 0.1e1 + t129;
	t109 = t92 * t118;
	t108 = t134 * t143;
	t107 = t134 * t141;
	t106 = t83 * t113;
	t103 = t111 * t92;
	t102 = t110 * t84;
	t70 = 0.1e1 / t72;
	t58 = 0.1e1 / t60;
	t56 = (t140 * t92 * t68 - t69 * t106) * t84;
	t55 = -t83 * t132 + t131 + (-t69 * t127 + t132) * t67;
	t54 = -t110 * t108 + (qJD(1) * t102 + t100 * t143) * t74;
	t1 = [t84 * t89 * t107 + (-t83 * t89 * t121 + qJD(3) * t102) * t74, 0, t54, 0, 0; (t62 * t109 + (-t62 * t119 + (qJD(1) * t56 + t53) * t135) * t58) * t83 + (t63 * t109 * t56 + (-((t57 * t106 + t140 * t119 + t107) * t68 + (t108 * t125 - t137 + (t137 + (-t86 * t90 + t141) * t120) * t74) * t69) * t115 + (t92 * t117 - t63 * t119) * t56 + (-t62 + ((-t77 + t82) * t69 * t113 + t140 * t114) * t63) * t121) * t58) * t84, 0, (t55 * t135 - t62 * t93) * t84 * t118 + ((-t62 * t123 + (-qJD(3) * t55 - t53) * t136) * t93 + (-t84 * qJD(3) * t62 - (-t54 * t69 * t83 + t68 * t120 + t133 * t138 - t138 + (-qJD(3) * t68 - t122 * t69) * t67) * t115 + (t84 * t117 + t63 * t123) * t55 - ((t54 - t122) * t68 + ((0.1e1 - t133) * qJD(3) + (t67 - t83) * t57) * t69) * t93 * t136) * t92) * t58, 0, 0; t111 * t93 * t116 + (qJD(3) * t103 + 0.2e1 * t93 * t142) * t70, 0, t78 * t116 * t126 + (qJD(1) * t103 - t119 * t130) * t70, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:54:37
	% EndTime: 2019-12-29 16:54:38
	% DurationCPUTime: 1.54s
	% Computational Cost: add. (1787->89), mult. (2519->200), div. (480->12), fcn. (2968->9), ass. (0->92)
	t122 = qJ(1) + pkin(8);
	t120 = sin(t122);
	t171 = qJD(3) * t120;
	t118 = t120 ^ 2;
	t129 = sin(qJ(3));
	t124 = 0.1e1 / t129 ^ 2;
	t131 = cos(qJ(3));
	t127 = t131 ^ 2;
	t177 = t124 * t127;
	t116 = t118 * t177 + 0.1e1;
	t114 = 0.1e1 / t116;
	t123 = 0.1e1 / t129;
	t159 = t124 * t171;
	t121 = cos(t122);
	t172 = qJD(1) * t131;
	t160 = t121 * t172;
	t170 = qJD(3) * t129;
	t88 = ((t120 * t170 - t160) * t123 + t127 * t159) * t114;
	t149 = -t88 + t171;
	t150 = -t120 * t88 + qJD(3);
	t128 = sin(qJ(5));
	t130 = cos(qJ(5));
	t148 = qJD(5) * t129 + qJD(1);
	t169 = qJD(3) * t131;
	t193 = t128 * t148 - t130 * t169;
	t192 = t128 * t169 + t130 * t148;
	t179 = t120 * t131;
	t113 = atan2(-t179, t129);
	t112 = cos(t113);
	t111 = sin(t113);
	t163 = t111 * t179;
	t101 = t112 * t129 - t163;
	t95 = 0.1e1 / t101;
	t176 = t128 * t129;
	t108 = t120 * t130 + t121 * t176;
	t104 = 0.1e1 / t108;
	t96 = 0.1e1 / t101 ^ 2;
	t105 = 0.1e1 / t108 ^ 2;
	t191 = t114 - 0.1e1;
	t180 = t112 * t131;
	t83 = t150 * t180 + (t149 * t129 - t160) * t111;
	t190 = t83 * t95 * t96;
	t175 = t129 * t130;
	t107 = t120 * t128 - t121 * t175;
	t103 = t107 ^ 2;
	t102 = t103 * t105 + 0.1e1;
	t182 = t105 * t107;
	t147 = qJD(1) * t129 + qJD(5);
	t143 = t147 * t128;
	t90 = -t120 * t143 + t192 * t121;
	t186 = t104 * t105 * t90;
	t142 = t147 * t130;
	t89 = t120 * t142 + t193 * t121;
	t189 = (-t103 * t186 + t89 * t182) / t102 ^ 2;
	t119 = t121 ^ 2;
	t93 = t119 * t127 * t96 + 0.1e1;
	t91 = 0.1e1 / t93;
	t188 = t91 * t96;
	t187 = t95 * t91;
	t184 = t121 * t96;
	t181 = t111 * t129;
	t178 = t123 * t127;
	t174 = qJD(1) * t120;
	t173 = qJD(1) * t121;
	t145 = t120 * t127 * t173;
	t164 = t96 * t170;
	t168 = 0.2e1 * (-t96 * t145 + (-t127 * t190 - t131 * t164) * t119) / t93 ^ 2;
	t167 = 0.2e1 * t190;
	t166 = 0.2e1 * t189;
	t126 = t131 * t127;
	t140 = qJD(3) * (-t124 * t126 - t131) * t123;
	t165 = 0.2e1 / t116 ^ 2 * (t118 * t140 + t124 * t145);
	t162 = t114 * t178;
	t161 = t120 * t172;
	t156 = t95 * t168;
	t155 = t96 * t168;
	t154 = 0.1e1 + t177;
	t153 = 0.2e1 * t107 * t186;
	t152 = t120 * t165;
	t151 = t131 * t165;
	t146 = t120 * t162;
	t144 = t154 * t121;
	t141 = t104 * t130 + t128 * t182;
	t98 = 0.1e1 / t102;
	t139 = t141 * t98;
	t110 = -t120 * t176 + t121 * t130;
	t109 = t120 * t175 + t121 * t128;
	t100 = t154 * t120 * t114;
	t87 = (t191 * t131 * t111 + t112 * t146) * t121;
	t86 = t120 * t181 + t180 + (-t112 * t179 - t181) * t100;
	t84 = -t154 * t152 + (qJD(1) * t144 + 0.2e1 * t120 * t140) * t114;
	t1 = [t121 * t123 * t151 + (qJD(3) * t144 + t123 * t161) * t114, 0, t84, 0, 0; (t170 * t187 + (t156 + (qJD(1) * t87 + t83) * t188) * t131) * t120 + (t87 * t155 * t131 + (t87 * t164 + (t87 * t167 + ((t88 * t146 + t191 * t170 + t151) * t111 + (t152 * t178 + t131 * t88 + (t126 * t159 - (t88 - 0.2e1 * t171) * t131) * t114) * t112) * t184) * t131 + (-t95 + (-(-t118 + t119) * t112 * t162 + t191 * t163) * t96) * t172) * t91) * t121, 0, (t174 * t187 + (t156 + (qJD(3) * t86 + t83) * t188) * t121) * t129 + (t86 * t121 * t155 + (-t121 * qJD(3) * t95 + (t121 * t167 + t96 * t174) * t86 + (-((-t100 * t173 - t120 * t84) * t112 + (-t150 * t100 + t149) * t111) * t131 - ((-t84 + t173) * t111 + (t149 * t100 - t150) * t112) * t129) * t184) * t91) * t131, 0, 0; (-t104 * t109 + t110 * t182) * t166 + (t110 * t153 + (-t109 * t90 - t110 * t89 + (t192 * t120 + t121 * t143) * t107) * t105 + (-t193 * t120 + t121 * t142) * t104) * t98, 0, t139 * t161 + (t139 * t170 + (t141 * t166 + ((qJD(5) * t104 + t153) * t128 + (-t128 * t89 + (-qJD(5) * t107 + t90) * t130) * t105) * t98) * t131) * t121, 0, -0.2e1 * t189 + 0.2e1 * (t105 * t89 * t98 + (-t105 * t189 - t98 * t186) * t107) * t107;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end