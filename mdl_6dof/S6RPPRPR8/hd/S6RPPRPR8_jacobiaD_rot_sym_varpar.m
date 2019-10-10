% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRPR8
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
%   Wie in S6RPPRPR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRPR8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:04
	% EndTime: 2019-10-09 23:46:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:05
	% EndTime: 2019-10-09 23:46:05
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (1694->69), mult. (1839->158), div. (436->14), fcn. (2165->7), ass. (0->74)
	t93 = cos(qJ(1));
	t141 = 0.2e1 * t93;
	t89 = 0.1e1 / t93;
	t92 = sin(qJ(1));
	t123 = t89 * t92;
	t87 = t92 ^ 2;
	t88 = t93 ^ 2;
	t101 = qJD(1) * (t87 / t88 + 0.1e1) * t123;
	t85 = pkin(9) + qJ(4);
	t83 = sin(t85);
	t119 = qJD(4) * t83;
	t84 = cos(t85);
	t104 = t84 * t87 * t119;
	t90 = 0.1e1 / t93 ^ 2;
	t124 = t87 * t90;
	t77 = t83 ^ 2;
	t76 = t77 * t124 + 0.1e1;
	t140 = -0.2e1 * (t77 * t101 + t90 * t104) / t76 ^ 2;
	t117 = qJD(4) * t93;
	t121 = qJD(1) * t92;
	t111 = t84 * t121;
	t79 = 0.1e1 / t83 ^ 2;
	t82 = t84 ^ 2;
	t126 = t79 * t82;
	t73 = t88 * t126 + 0.1e1;
	t71 = 0.1e1 / t73;
	t78 = 0.1e1 / t83;
	t57 = ((t83 * t117 + t111) * t78 + t117 * t126) * t71;
	t139 = -t57 + t117;
	t122 = t93 * t84;
	t70 = atan2(-t122, t83);
	t68 = sin(t70);
	t69 = cos(t70);
	t64 = -t68 * t122 + t69 * t83;
	t61 = 0.1e1 / t64;
	t62 = 0.1e1 / t64 ^ 2;
	t138 = t71 - 0.1e1;
	t120 = qJD(1) * t93;
	t105 = t82 * t92 * t120;
	t125 = t82 * t87;
	t129 = t69 * t84;
	t53 = (-t57 * t93 + qJD(4)) * t129 + (t139 * t83 + t111) * t68;
	t136 = t53 * t61 * t62;
	t60 = t62 * t125 + 0.1e1;
	t137 = (-t125 * t136 + (-t104 + t105) * t62) / t60 ^ 2;
	t135 = t57 * t84;
	t134 = t62 * t84;
	t133 = t62 * t92;
	t127 = t78 * t84;
	t81 = t84 * t82;
	t100 = qJD(4) * (-t78 / t77 * t81 - t127);
	t132 = (t88 * t100 - t79 * t105) / t73 ^ 2;
	t130 = t68 * t93;
	t128 = t78 * t82;
	t118 = qJD(4) * t92;
	t116 = -0.2e1 * t136;
	t115 = t84 * t137;
	t114 = t84 * t133;
	t113 = t84 * t132;
	t112 = t71 * t128;
	t110 = 0.1e1 + t126;
	t109 = 0.1e1 + t124;
	t108 = t132 * t141;
	t107 = t93 * t112;
	t106 = t138 * t84 * t68;
	t103 = t110 * t92;
	t102 = t109 * t84;
	t74 = 0.1e1 / t76;
	t66 = t110 * t93 * t71;
	t58 = 0.1e1 / t60;
	t56 = (-t69 * t107 - t106) * t92;
	t55 = t83 * t130 + t129 + (-t69 * t122 - t68 * t83) * t66;
	t54 = -t110 * t108 + (-qJD(1) * t103 + t100 * t141) * t71;
	t1 = [-0.2e1 * t92 * t78 * t113 + (-qJD(4) * t103 + t120 * t127) * t71, 0, 0, t54, 0, 0; (0.2e1 * t61 * t115 + (t61 * t119 + (qJD(1) * t56 + t53) * t134) * t58) * t93 + (-0.2e1 * t62 * t115 * t56 + (((t57 * t107 + t138 * t119 + 0.2e1 * t113) * t68 + (t108 * t128 + t135 + (-t135 + (t79 * t81 + 0.2e1 * t84) * t117) * t71) * t69) * t114 + (t84 * t116 - t62 * t119) * t56 + (t61 + ((t87 - t88) * t69 * t112 - t93 * t106) * t62) * t84 * qJD(1)) * t58) * t92, 0, 0, 0.2e1 * (-t55 * t134 - t61 * t83) * t92 * t137 + ((t61 * t120 + (-qJD(4) * t55 - t53) * t133) * t83 + (t61 * t118 + (-t54 * t69 * t93 + t139 * t68 + (-qJD(4) * t68 + t121 * t69 + t130 * t57) * t66) * t114 + (t92 * t116 + t62 * t120) * t55 + ((-t54 - t121) * t68 + ((t66 * t93 - 0.1e1) * qJD(4) + (-t66 + t93) * t57) * t69) * t83 * t133) * t84) * t58, 0, 0; t109 * t83 * t140 + (qJD(4) * t102 + 0.2e1 * t101 * t83) * t74, 0, 0, t84 * t123 * t140 + (-t83 * t89 * t118 + qJD(1) * t102) * t74, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:46:05
	% EndTime: 2019-10-09 23:46:06
	% DurationCPUTime: 1.00s
	% Computational Cost: add. (1892->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->94)
	t137 = cos(qJ(1));
	t133 = t137 ^ 2;
	t131 = pkin(9) + qJ(4);
	t129 = sin(t131);
	t125 = t129 ^ 2;
	t130 = cos(t131);
	t127 = 0.1e1 / t130 ^ 2;
	t182 = t125 * t127;
	t121 = t133 * t182 + 0.1e1;
	t124 = t129 * t125;
	t126 = 0.1e1 / t130;
	t180 = t126 * t129;
	t145 = qJD(4) * (t124 * t126 * t127 + t180);
	t135 = sin(qJ(1));
	t172 = qJD(1) * t137;
	t149 = t125 * t135 * t172;
	t189 = (-t127 * t149 + t133 * t145) / t121 ^ 2;
	t200 = -0.2e1 * t189;
	t158 = 0.1e1 + t182;
	t199 = t137 * t158;
	t152 = qJD(1) * t130 + qJD(6);
	t169 = qJD(4) * t137;
	t198 = t129 * t169 + t152 * t135;
	t174 = t137 * t129;
	t120 = atan2(t174, t130);
	t116 = sin(t120);
	t117 = cos(t120);
	t102 = t116 * t174 + t117 * t130;
	t99 = 0.1e1 / t102;
	t136 = cos(qJ(6));
	t175 = t136 * t137;
	t134 = sin(qJ(6));
	t177 = t135 * t134;
	t115 = -t130 * t177 + t175;
	t109 = 0.1e1 / t115;
	t100 = 0.1e1 / t102 ^ 2;
	t110 = 0.1e1 / t115 ^ 2;
	t118 = 0.1e1 / t121;
	t197 = t118 - 0.1e1;
	t132 = t135 ^ 2;
	t171 = qJD(4) * t130;
	t181 = t125 * t132;
	t161 = t127 * t169;
	t173 = qJD(1) * t135;
	t162 = t129 * t173;
	t93 = ((t130 * t169 - t162) * t126 + t125 * t161) * t118;
	t154 = -t93 + t169;
	t155 = t137 * t93 - qJD(4);
	t184 = t117 * t129;
	t88 = t155 * t184 + (t154 * t130 - t162) * t116;
	t194 = t99 * t100 * t88;
	t96 = t100 * t181 + 0.1e1;
	t196 = (-t181 * t194 + (t129 * t132 * t171 + t149) * t100) / t96 ^ 2;
	t94 = 0.1e1 / t96;
	t195 = t100 * t94;
	t176 = t135 * t136;
	t178 = t134 * t137;
	t114 = t130 * t176 + t178;
	t108 = t114 ^ 2;
	t107 = t108 * t110 + 0.1e1;
	t187 = t110 * t114;
	t153 = qJD(6) * t130 + qJD(1);
	t170 = qJD(4) * t135;
	t98 = -t153 * t176 + (t129 * t170 - t152 * t137) * t134;
	t192 = t109 * t110 * t98;
	t163 = t130 * t175;
	t97 = -qJD(1) * t163 - qJD(6) * t175 + (qJD(4) * t129 * t136 + t153 * t134) * t135;
	t193 = 0.1e1 / t107 ^ 2 * (-t108 * t192 - t97 * t187);
	t188 = t109 * t136;
	t186 = t114 * t134;
	t185 = t116 * t130;
	t183 = t125 * t126;
	t179 = t129 * t135;
	t168 = 0.2e1 * t196;
	t167 = 0.2e1 * t194;
	t166 = 0.2e1 * t193;
	t165 = t94 * t171;
	t164 = t137 * t183;
	t159 = -0.2e1 * t99 * t196;
	t157 = 0.2e1 * t114 * t192;
	t156 = 0.2e1 * t129 * t189;
	t151 = t118 * t164;
	t150 = t197 * t129 * t116;
	t148 = t158 * t135;
	t147 = t153 * t137;
	t146 = t110 * t186 + t188;
	t113 = -t130 * t178 - t176;
	t112 = t163 - t177;
	t105 = 0.1e1 / t107;
	t104 = t118 * t199;
	t92 = (-t117 * t151 + t150) * t135;
	t90 = t137 * t185 - t184 + (t117 * t174 - t185) * t104;
	t89 = t199 * t200 + (-qJD(1) * t148 + 0.2e1 * t137 * t145) * t118;
	t1 = [t135 * t126 * t156 + (-qJD(4) * t148 - t172 * t180) * t118, 0, 0, t89, 0, 0; (t99 * t165 + (t159 + (-qJD(1) * t92 - t88) * t195) * t129) * t137 + ((-t92 * t165 + (t92 * t168 + ((-t93 * t151 - t197 * t171 + t156) * t116 + (t164 * t200 + t129 * t93 + (t124 * t161 - (t93 - 0.2e1 * t169) * t129) * t118) * t117) * t94 * t135) * t129) * t100 + (t92 * t167 + (-t99 + ((-t132 + t133) * t118 * t117 * t183 - t137 * t150) * t100) * qJD(1)) * t129 * t94) * t135, 0, 0, (t99 * t94 * t172 + (t159 + (-qJD(4) * t90 - t88) * t195) * t135) * t130 + (((-qJD(4) * t99 + t90 * t167) * t135 + (-t90 * t172 + (-(-t104 * t173 + t137 * t89) * t117 - ((-t104 * t137 + 0.1e1) * t93 + (t104 - t137) * qJD(4)) * t116) * t179) * t100) * t94 + (t90 * t168 - ((-t89 - t173) * t116 + (t154 * t104 + t155) * t117) * t94 * t130) * t100 * t135) * t129, 0, 0; (-t109 * t112 + t113 * t187) * t166 + (t113 * t157 - t109 * t134 * t147 - t198 * t188 + (t114 * t136 * t147 - t112 * t98 + t113 * t97 - t198 * t186) * t110) * t105, 0, 0, t146 * t166 * t179 + (-t146 * t130 * t170 + (-t146 * t172 + ((qJD(6) * t109 + t157) * t134 + (t134 * t97 + (-qJD(6) * t114 + t98) * t136) * t110) * t135) * t129) * t105, 0, -0.2e1 * t193 + 0.2e1 * (-t97 * t110 * t105 + (-t105 * t192 - t110 * t193) * t114) * t114;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end