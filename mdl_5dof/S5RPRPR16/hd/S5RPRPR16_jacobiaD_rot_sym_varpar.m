% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPR16
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
%   Wie in S5RPRPR16_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRPR16_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR16_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:13:28
	% EndTime: 2019-12-29 17:13:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:13:28
	% EndTime: 2019-12-29 17:13:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:13:33
	% EndTime: 2019-12-29 17:13:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:13:28
	% EndTime: 2019-12-29 17:13:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:13:28
	% EndTime: 2019-12-29 17:13:29
	% DurationCPUTime: 1.05s
	% Computational Cost: add. (585->68), mult. (1839->158), div. (436->14), fcn. (2165->7), ass. (0->74)
	t85 = sin(qJ(1));
	t116 = qJD(1) * t85;
	t87 = cos(qJ(1));
	t137 = 0.2e1 * t87;
	t77 = t85 ^ 2;
	t82 = 0.1e1 / t87 ^ 2;
	t119 = t77 * t82;
	t84 = sin(qJ(3));
	t72 = t84 ^ 2;
	t70 = t72 * t119 + 0.1e1;
	t80 = t87 ^ 2;
	t81 = 0.1e1 / t87;
	t95 = (t77 / t80 + 0.1e1) * t81 * t116;
	t113 = qJD(3) * t84;
	t86 = cos(qJ(3));
	t98 = t77 * t86 * t113;
	t136 = -0.2e1 * (t72 * t95 + t82 * t98) / t70 ^ 2;
	t111 = qJD(3) * t87;
	t115 = qJD(1) * t86;
	t105 = t85 * t115;
	t74 = 0.1e1 / t84 ^ 2;
	t79 = t86 ^ 2;
	t121 = t74 * t79;
	t71 = t80 * t121 + 0.1e1;
	t68 = 0.1e1 / t71;
	t73 = 0.1e1 / t84;
	t52 = ((t84 * t111 + t105) * t73 + t111 * t121) * t68;
	t134 = -t52 + t111;
	t117 = t87 * t86;
	t65 = atan2(-t117, t84);
	t63 = sin(t65);
	t64 = cos(t65);
	t59 = -t63 * t117 + t64 * t84;
	t56 = 0.1e1 / t59;
	t57 = 0.1e1 / t59 ^ 2;
	t133 = t68 - 0.1e1;
	t120 = t77 * t79;
	t124 = t64 * t86;
	t48 = (-t52 * t87 + qJD(3)) * t124 + (t134 * t84 + t105) * t63;
	t131 = t48 * t56 * t57;
	t55 = t57 * t120 + 0.1e1;
	t114 = qJD(1) * t87;
	t99 = t79 * t85 * t114;
	t132 = (-t120 * t131 + (-t98 + t99) * t57) / t55 ^ 2;
	t130 = t52 * t86;
	t129 = t57 * t85;
	t128 = t57 * t86;
	t122 = t73 * t86;
	t78 = t86 * t79;
	t94 = qJD(3) * (-t73 / t72 * t78 - t122);
	t127 = (-t74 * t99 + t80 * t94) / t71 ^ 2;
	t125 = t63 * t87;
	t123 = t73 * t79;
	t118 = t85 * t86;
	t112 = qJD(3) * t85;
	t110 = -0.2e1 * t131;
	t109 = t86 * t132;
	t108 = t57 * t118;
	t107 = t86 * t127;
	t106 = t68 * t123;
	t104 = 0.1e1 + t121;
	t103 = 0.1e1 + t119;
	t102 = t127 * t137;
	t101 = t87 * t106;
	t100 = t133 * t86 * t63;
	t97 = t104 * t85;
	t96 = t103 * t86;
	t66 = 0.1e1 / t70;
	t62 = t104 * t87 * t68;
	t53 = 0.1e1 / t55;
	t51 = (-t64 * t101 - t100) * t85;
	t50 = t84 * t125 + t124 + (-t64 * t117 - t63 * t84) * t62;
	t49 = -t104 * t102 + (-qJD(1) * t97 + t94 * t137) * t68;
	t1 = [-0.2e1 * t85 * t73 * t107 + (-qJD(3) * t97 + t114 * t122) * t68, 0, t49, 0, 0; (0.2e1 * t56 * t109 + (t56 * t113 + (qJD(1) * t51 + t48) * t128) * t53) * t87 + (-0.2e1 * t57 * t109 * t51 + (((t52 * t101 + t133 * t113 + 0.2e1 * t107) * t63 + (t102 * t123 + t130 + (-t130 + (t74 * t78 + 0.2e1 * t86) * t111) * t68) * t64) * t108 + (t86 * t110 - t57 * t113) * t51 + (t56 + ((t77 - t80) * t64 * t106 - t87 * t100) * t57) * t115) * t53) * t85, 0, 0.2e1 * (-t50 * t128 - t56 * t84) * t85 * t132 + ((t56 * t114 + (-qJD(3) * t50 - t48) * t129) * t84 + (t56 * t112 + (-t49 * t64 * t87 + t134 * t63 + (-qJD(3) * t63 + t116 * t64 + t125 * t52) * t62) * t108 + (t110 * t85 + t57 * t114) * t50 + ((-t49 - t116) * t63 + ((t62 * t87 - 0.1e1) * qJD(3) + (-t62 + t87) * t52) * t64) * t84 * t129) * t86) * t53, 0, 0; t103 * t84 * t136 + (qJD(3) * t96 + 0.2e1 * t84 * t95) * t66, 0, t81 * t118 * t136 + (-t112 * t81 * t84 + qJD(1) * t96) * t66, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:13:28
	% EndTime: 2019-12-29 17:13:30
	% DurationCPUTime: 1.46s
	% Computational Cost: add. (624->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->93)
	t129 = cos(qJ(1));
	t123 = t129 ^ 2;
	t125 = sin(qJ(3));
	t118 = t125 ^ 2;
	t128 = cos(qJ(3));
	t121 = 0.1e1 / t128 ^ 2;
	t171 = t118 * t121;
	t114 = t123 * t171 + 0.1e1;
	t117 = t125 * t118;
	t120 = 0.1e1 / t128;
	t170 = t120 * t125;
	t138 = qJD(3) * (t117 * t120 * t121 + t170);
	t126 = sin(qJ(1));
	t163 = qJD(1) * t129;
	t153 = t126 * t163;
	t180 = 0.1e1 / t114 ^ 2 * (t123 * t138 - t153 * t171);
	t192 = -0.2e1 * t180;
	t111 = 0.1e1 / t114;
	t148 = 0.1e1 + t171;
	t190 = t129 * t148;
	t99 = t111 * t190;
	t191 = -t129 * t99 + 0.1e1;
	t144 = qJD(1) * t128 + qJD(5);
	t160 = qJD(3) * t129;
	t189 = t125 * t160 + t144 * t126;
	t166 = t129 * t125;
	t113 = atan2(t166, t128);
	t109 = sin(t113);
	t110 = cos(t113);
	t95 = t109 * t166 + t110 * t128;
	t92 = 0.1e1 / t95;
	t124 = sin(qJ(5));
	t127 = cos(qJ(5));
	t165 = t129 * t127;
	t168 = t126 * t128;
	t108 = -t124 * t168 + t165;
	t102 = 0.1e1 / t108;
	t103 = 0.1e1 / t108 ^ 2;
	t93 = 0.1e1 / t95 ^ 2;
	t188 = t111 - 0.1e1;
	t119 = t126 ^ 2;
	t161 = qJD(3) * t128;
	t157 = t93 * t161;
	t164 = qJD(1) * t126;
	t154 = t125 * t164;
	t173 = t110 * t125;
	t152 = t121 * t160;
	t86 = ((t128 * t160 - t154) * t120 + t118 * t152) * t111;
	t81 = (t129 * t86 - qJD(3)) * t173 + (-t154 + (-t86 + t160) * t128) * t109;
	t186 = t81 * t92 * t93;
	t91 = t119 * t118 * t93 + 0.1e1;
	t187 = (t119 * t125 * t157 + (-t119 * t186 + t93 * t153) * t118) / t91 ^ 2;
	t167 = t129 * t124;
	t107 = t127 * t168 + t167;
	t101 = t107 ^ 2;
	t100 = t101 * t103 + 0.1e1;
	t175 = t103 * t107;
	t145 = qJD(5) * t128 + qJD(1);
	t162 = qJD(3) * t126;
	t169 = t126 * t127;
	t88 = -t145 * t169 + (t125 * t162 - t144 * t129) * t124;
	t182 = t102 * t103 * t88;
	t155 = t128 * t165;
	t87 = -qJD(1) * t155 - qJD(5) * t165 + (qJD(3) * t125 * t127 + t145 * t124) * t126;
	t185 = (-t101 * t182 - t87 * t175) / t100 ^ 2;
	t89 = 0.1e1 / t91;
	t184 = t89 * t93;
	t183 = t92 * t89;
	t179 = t126 * t93;
	t177 = qJD(3) * t99;
	t176 = t102 * t127;
	t174 = t107 * t124;
	t172 = t118 * t120;
	t159 = 0.2e1 * t186;
	t158 = 0.2e1 * t185;
	t156 = t129 * t172;
	t150 = -0.2e1 * t92 * t187;
	t149 = 0.2e1 * t93 * t187;
	t147 = 0.2e1 * t107 * t182;
	t146 = 0.2e1 * t125 * t180;
	t143 = t111 * t156;
	t142 = t188 * t125 * t109;
	t141 = t148 * t126;
	t140 = t145 * t129;
	t139 = t103 * t174 + t176;
	t97 = 0.1e1 / t100;
	t137 = t139 * t97;
	t106 = -t128 * t167 - t169;
	t105 = -t126 * t124 + t155;
	t85 = (-t110 * t143 + t142) * t126;
	t84 = -t191 * t173 + (t129 - t99) * t128 * t109;
	t82 = t190 * t192 + (-qJD(1) * t141 + 0.2e1 * t129 * t138) * t111;
	t1 = [t126 * t120 * t146 + (-qJD(3) * t141 - t163 * t170) * t111, 0, t82, 0, 0; (t161 * t183 + (t150 + (-qJD(1) * t85 - t81) * t184) * t125) * t129 + (t85 * t149 * t125 + (-t85 * t157 + (t85 * t159 + ((-t86 * t143 - t188 * t161 + t146) * t109 + (t156 * t192 + t125 * t86 + (t117 * t152 - (t86 - 0.2e1 * t160) * t125) * t111) * t110) * t179) * t125 + (-t92 + (-(t119 - t123) * t111 * t110 * t172 - t129 * t142) * t93) * t125 * qJD(1)) * t89) * t126, 0, (t163 * t183 + (t150 + (-qJD(3) * t84 - t81) * t184) * t126) * t128 + (t84 * t126 * t149 + (-t92 * t162 - ((t129 * t82 - t164 * t99) * t110 + (t191 * t86 - t160 + t177) * t109) * t125 * t179 + (t126 * t159 - t93 * t163) * t84) * t89 - ((-t82 - t164) * t109 + (-t86 * t99 - qJD(3) + (t86 + t177) * t129) * t110) * t168 * t184) * t125, 0, 0; (-t102 * t105 + t106 * t175) * t158 + (t106 * t147 - t102 * t124 * t140 - t189 * t176 + (t107 * t127 * t140 - t105 * t88 + t106 * t87 - t189 * t174) * t103) * t97, 0, -t126 * t137 * t161 + (-t137 * t163 + (t139 * t158 + ((qJD(5) * t102 + t147) * t124 + (t124 * t87 + (-qJD(5) * t107 + t88) * t127) * t103) * t97) * t126) * t125, 0, -0.2e1 * t185 + 0.2e1 * (-t87 * t103 * t97 + (-t103 * t185 - t97 * t182) * t107) * t107;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end