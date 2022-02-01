% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPR5
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
%   Wie in S5RPRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRPR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:47
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (343->39), mult. (853->108), div. (126->12), fcn. (1047->9), ass. (0->56)
	t102 = cos(pkin(8));
	t105 = cos(qJ(3));
	t106 = cos(qJ(1));
	t121 = t106 * t105;
	t103 = sin(qJ(3));
	t104 = sin(qJ(1));
	t124 = t104 * t103;
	t84 = t102 * t121 + t124;
	t101 = sin(pkin(8));
	t125 = t104 * t101;
	t87 = atan2(-t125, -t102);
	t85 = sin(t87);
	t118 = t85 * t125;
	t86 = cos(t87);
	t71 = -t86 * t102 - t118;
	t68 = 0.1e1 / t71;
	t78 = 0.1e1 / t84;
	t96 = 0.1e1 / t102;
	t97 = 0.1e1 / t102 ^ 2;
	t69 = 0.1e1 / t71 ^ 2;
	t79 = 0.1e1 / t84 ^ 2;
	t95 = t101 ^ 2;
	t99 = t104 ^ 2;
	t90 = t99 * t95 * t97 + 0.1e1;
	t88 = 0.1e1 / t90;
	t135 = t88 - 0.1e1;
	t122 = t106 * t103;
	t123 = t104 * t105;
	t82 = -t102 * t123 + t122;
	t83 = t102 * t122 - t123;
	t73 = t82 * qJD(1) - t83 * qJD(3);
	t132 = t73 * t78 * t79;
	t115 = t102 * t124 + t121;
	t72 = t115 * qJD(1) - t84 * qJD(3);
	t133 = t72 * t79;
	t77 = t83 ^ 2;
	t76 = t77 * t79 + 0.1e1;
	t134 = (-t77 * t132 - t83 * t133) / t76 ^ 2;
	t131 = t82 * t83;
	t89 = 0.1e1 / t90 ^ 2;
	t130 = t89 * t101 * t95;
	t129 = t95 * t96;
	t100 = t106 ^ 2;
	t128 = t100 * t69;
	t127 = t104 * t69;
	t126 = t106 * t69;
	t120 = qJD(1) * t104;
	t119 = t88 * t129;
	t116 = t86 * t119;
	t64 = (t135 * t85 * t101 - t104 * t116) * t106;
	t98 = t96 * t97;
	t74 = 0.1e1 / t76;
	t70 = t68 * t69;
	t67 = t95 * t128 + 0.1e1;
	t63 = qJD(1) * t64;
	t1 = [(-0.2e1 * t100 * t98 * t130 - t101 * t88 * t96) * t120, 0, 0, 0, 0; (0.2e1 * (t104 * t68 + t64 * t126) / t67 ^ 2 * (-t100 * t63 * t70 - t120 * t126) * t95 + ((0.2e1 * t106 * t64 * t70 + t127) * t63 + (t64 * t127 + (-t68 - (-t104 * t85 * t97 * t130 + (-0.2e1 * t119 + (0.2e1 * t95 ^ 2 * t98 * t99 + t129) * t89) * t86) * t128 + (-t99 * t116 + t135 * t118) * t69) * t106) * qJD(1)) / t67) * t101, 0, 0, 0, 0; 0.2e1 * (t115 * t78 + t79 * t131) * t134 + ((-t83 * qJD(1) + t82 * qJD(3)) * t78 + 0.2e1 * t131 * t132 + (t115 * t73 - (-t84 * qJD(1) + t115 * qJD(3)) * t83 + t82 * t72) * t79) * t74, 0, -0.2e1 * t134 + 0.2e1 * (-t74 * t133 + (-t74 * t132 - t79 * t134) * t83) * t83, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:47
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (537->40), mult. (853->108), div. (126->12), fcn. (1047->9), ass. (0->57)
	t110 = cos(pkin(8));
	t104 = 0.1e1 / t110 ^ 2;
	t109 = sin(pkin(8));
	t102 = t109 ^ 2;
	t111 = sin(qJ(1));
	t107 = t111 ^ 2;
	t96 = t102 * t107 * t104 + 0.1e1;
	t95 = 0.1e1 / t96 ^ 2;
	t141 = t104 * t95;
	t106 = qJ(3) + pkin(9);
	t100 = cos(t106);
	t112 = cos(qJ(1));
	t126 = t112 * t100;
	t99 = sin(t106);
	t132 = t111 * t99;
	t89 = t110 * t126 + t132;
	t128 = t109 * t111;
	t93 = atan2(-t128, -t110);
	t90 = sin(t93);
	t124 = t90 * t128;
	t91 = cos(t93);
	t81 = -t91 * t110 - t124;
	t78 = 0.1e1 / t81;
	t83 = 0.1e1 / t89;
	t103 = 0.1e1 / t110;
	t79 = 0.1e1 / t81 ^ 2;
	t84 = 0.1e1 / t89 ^ 2;
	t94 = 0.1e1 / t96;
	t140 = t94 - 0.1e1;
	t127 = t111 * t100;
	t130 = t112 * t99;
	t87 = -t110 * t127 + t130;
	t88 = t110 * t130 - t127;
	t74 = t87 * qJD(1) - t88 * qJD(3);
	t137 = t74 * t83 * t84;
	t121 = t110 * t132 + t126;
	t73 = t121 * qJD(1) - qJD(3) * t89;
	t138 = t73 * t84;
	t82 = t88 ^ 2;
	t77 = t82 * t84 + 0.1e1;
	t139 = (-t82 * t137 - t88 * t138) / t77 ^ 2;
	t136 = t87 * t88;
	t135 = t103 * t141;
	t108 = t112 ^ 2;
	t134 = t108 * t79;
	t133 = t111 * t79;
	t131 = t112 * t79;
	t129 = t102 * t103;
	t125 = qJD(1) * t111;
	t122 = t91 * t94 * t129;
	t69 = (t140 * t90 * t109 - t111 * t122) * t112;
	t101 = t109 * t102;
	t80 = t78 * t79;
	t75 = 0.1e1 / t77;
	t72 = t102 * t134 + 0.1e1;
	t68 = qJD(1) * t69;
	t1 = [(-0.2e1 * t101 * t108 * t135 - t103 * t109 * t94) * t125, 0, 0, 0, 0; (0.2e1 * (t111 * t78 + t69 * t131) / t72 ^ 2 * (-t108 * t68 * t80 - t125 * t131) * t102 + ((0.2e1 * t112 * t69 * t80 + t133) * t68 + (t69 * t133 + (-t78 - (-t101 * t111 * t90 * t141 + (0.2e1 * t107 * t102 ^ 2 * t135 + (-0.2e1 * t94 + t95) * t129) * t91) * t134 + (-t107 * t122 + t140 * t124) * t79) * t112) * qJD(1)) / t72) * t109, 0, 0, 0, 0; 0.2e1 * (t121 * t83 + t84 * t136) * t139 + ((-t88 * qJD(1) + t87 * qJD(3)) * t83 + 0.2e1 * t136 * t137 + (t121 * t74 - (-t89 * qJD(1) + t121 * qJD(3)) * t88 + t87 * t73) * t84) * t75, 0, -0.2e1 * t139 + 0.2e1 * (-t75 * t138 + (-t75 * t137 - t84 * t139) * t88) * t88, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:47
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (1005->41), mult. (1047->107), div. (144->12), fcn. (1257->9), ass. (0->55)
	t134 = cos(pkin(8));
	t128 = 0.1e1 / t134 ^ 2;
	t133 = sin(pkin(8));
	t126 = t133 ^ 2;
	t135 = sin(qJ(1));
	t131 = t135 ^ 2;
	t120 = t131 * t126 * t128 + 0.1e1;
	t136 = cos(qJ(1));
	t132 = t136 ^ 2;
	t155 = 0.1e1 / t120 ^ 2 * t132;
	t163 = t128 * t155;
	t118 = 0.1e1 / t120;
	t162 = (t118 - 0.1e1) * t133;
	t124 = qJ(3) + pkin(9) + qJ(5);
	t123 = cos(t124);
	t149 = t136 * t123;
	t122 = sin(t124);
	t153 = t135 * t122;
	t112 = t134 * t149 + t153;
	t151 = t135 * t133;
	t117 = atan2(-t151, -t134);
	t114 = sin(t117);
	t115 = cos(t117);
	t104 = -t114 * t151 - t115 * t134;
	t101 = 0.1e1 / t104;
	t106 = 0.1e1 / t112;
	t127 = 0.1e1 / t134;
	t102 = 0.1e1 / t104 ^ 2;
	t107 = 0.1e1 / t112 ^ 2;
	t150 = t136 * t122;
	t152 = t135 * t123;
	t111 = t134 * t150 - t152;
	t105 = t111 ^ 2;
	t130 = qJD(3) + qJD(5);
	t145 = t134 * t153 + t149;
	t93 = t145 * qJD(1) - t112 * t130;
	t158 = t93 * t107;
	t110 = -t134 * t152 + t150;
	t94 = qJD(1) * t110 - t111 * t130;
	t159 = t106 * t107 * t94;
	t97 = t105 * t107 + 0.1e1;
	t160 = (-t105 * t159 - t111 * t158) / t97 ^ 2;
	t157 = t102 * t136;
	t156 = t110 * t111;
	t154 = t126 * t127;
	t148 = qJD(1) * t135;
	t147 = t127 * t163;
	t92 = (-t115 * t118 * t135 * t154 + t114 * t162) * t136;
	t125 = t133 * t126;
	t103 = t101 * t102;
	t100 = t132 * t126 * t102 + 0.1e1;
	t95 = 0.1e1 / t97;
	t91 = qJD(1) * t92;
	t88 = -0.2e1 * t160 + 0.2e1 * (-t95 * t158 + (-t107 * t160 - t159 * t95) * t111) * t111;
	t1 = [(-t118 * t127 * t133 - 0.2e1 * t125 * t147) * t148, 0, 0, 0, 0; (0.2e1 * (t101 * t135 + t157 * t92) / t100 ^ 2 * (-t103 * t132 * t91 - t148 * t157) * t126 + ((0.2e1 * t103 * t136 * t92 + t102 * t135) * t91 + (-t136 * t101 + ((t92 + (t125 * t163 + t162) * t136 * t114) * t135 - (0.2e1 * t131 * t126 ^ 2 * t147 + (t155 + (t131 - 0.2e1 * t132) * t118) * t154) * t136 * t115) * t102) * qJD(1)) / t100) * t133, 0, 0, 0, 0; 0.2e1 * (t106 * t145 + t107 * t156) * t160 + ((-qJD(1) * t111 + t110 * t130) * t106 + 0.2e1 * t156 * t159 + (t145 * t94 - (-qJD(1) * t112 + t130 * t145) * t111 + t110 * t93) * t107) * t95, 0, t88, 0, t88;];
	JaD_rot = t1;
end