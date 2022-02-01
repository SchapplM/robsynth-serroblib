% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP4
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
%   Wie in S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:25
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
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (727->41), mult. (1047->106), div. (144->12), fcn. (1257->9), ass. (0->56)
	t129 = cos(pkin(8));
	t122 = 0.1e1 / t129 ^ 2;
	t128 = sin(pkin(8));
	t120 = t128 ^ 2;
	t130 = sin(qJ(1));
	t125 = t130 ^ 2;
	t115 = t125 * t120 * t122 + 0.1e1;
	t131 = cos(qJ(1));
	t126 = t131 ^ 2;
	t150 = 0.1e1 / t115 ^ 2 * t126;
	t159 = t122 * t150;
	t113 = 0.1e1 / t115;
	t158 = (t113 - 0.1e1) * t128;
	t127 = qJ(3) + qJ(4);
	t118 = cos(t127);
	t144 = t131 * t118;
	t117 = sin(t127);
	t148 = t130 * t117;
	t107 = t129 * t144 + t148;
	t146 = t130 * t128;
	t111 = atan2(-t146, -t129);
	t109 = sin(t111);
	t110 = cos(t111);
	t99 = -t109 * t146 - t110 * t129;
	t96 = 0.1e1 / t99;
	t101 = 0.1e1 / t107;
	t121 = 0.1e1 / t129;
	t102 = 0.1e1 / t107 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t145 = t131 * t117;
	t147 = t130 * t118;
	t106 = t129 * t145 - t147;
	t100 = t106 ^ 2;
	t124 = qJD(3) + qJD(4);
	t140 = t129 * t148 + t144;
	t88 = t140 * qJD(1) - t107 * t124;
	t152 = t88 * t102;
	t105 = -t129 * t147 + t145;
	t89 = t105 * qJD(1) - t106 * t124;
	t155 = t101 * t102 * t89;
	t95 = t100 * t102 + 0.1e1;
	t156 = (-t100 * t155 - t106 * t152) / t95 ^ 2;
	t154 = t130 * t97;
	t153 = t131 * t97;
	t151 = t105 * t106;
	t149 = t120 * t121;
	t143 = qJD(1) * t130;
	t142 = t121 * t159;
	t87 = (-t110 * t113 * t130 * t149 + t109 * t158) * t131;
	t119 = t128 * t120;
	t98 = t96 * t97;
	t93 = 0.1e1 / t95;
	t92 = t126 * t120 * t97 + 0.1e1;
	t86 = qJD(1) * t87;
	t83 = -0.2e1 * t156 + 0.2e1 * (-t93 * t152 + (-t102 * t156 - t93 * t155) * t106) * t106;
	t1 = [(-t113 * t121 * t128 - 0.2e1 * t119 * t142) * t143, 0, 0, 0, 0; (0.2e1 * (t130 * t96 + t87 * t153) / t92 ^ 2 * (-t126 * t86 * t98 - t143 * t153) * t120 + ((0.2e1 * t131 * t87 * t98 + t154) * t86 + (t87 * t154 + (-t96 + (t119 * t159 + t158) * t109 * t154 - (0.2e1 * t125 * t120 ^ 2 * t142 + (t150 + (t125 - 0.2e1 * t126) * t113) * t149) * t97 * t110) * t131) * qJD(1)) / t92) * t128, 0, 0, 0, 0; 0.2e1 * (t101 * t140 + t102 * t151) * t156 + ((-t106 * qJD(1) + t105 * t124) * t101 + 0.2e1 * t151 * t155 + (t140 * t89 - (-t107 * qJD(1) + t140 * t124) * t106 + t105 * t88) * t102) * t93, 0, t83, t83, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:25
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (727->41), mult. (1047->107), div. (144->12), fcn. (1257->9), ass. (0->55)
	t133 = cos(pkin(8));
	t126 = 0.1e1 / t133 ^ 2;
	t132 = sin(pkin(8));
	t124 = t132 ^ 2;
	t134 = sin(qJ(1));
	t129 = t134 ^ 2;
	t119 = t129 * t124 * t126 + 0.1e1;
	t135 = cos(qJ(1));
	t130 = t135 ^ 2;
	t154 = 0.1e1 / t119 ^ 2 * t130;
	t162 = t126 * t154;
	t117 = 0.1e1 / t119;
	t161 = (t117 - 0.1e1) * t132;
	t131 = qJ(3) + qJ(4);
	t122 = cos(t131);
	t148 = t135 * t122;
	t121 = sin(t131);
	t152 = t134 * t121;
	t111 = t133 * t148 + t152;
	t150 = t134 * t132;
	t115 = atan2(-t150, -t133);
	t113 = sin(t115);
	t114 = cos(t115);
	t103 = -t113 * t150 - t114 * t133;
	t100 = 0.1e1 / t103;
	t105 = 0.1e1 / t111;
	t125 = 0.1e1 / t133;
	t101 = 0.1e1 / t103 ^ 2;
	t106 = 0.1e1 / t111 ^ 2;
	t149 = t135 * t121;
	t151 = t134 * t122;
	t110 = t133 * t149 - t151;
	t104 = t110 ^ 2;
	t128 = qJD(3) + qJD(4);
	t144 = t133 * t152 + t148;
	t92 = qJD(1) * t144 - t111 * t128;
	t157 = t92 * t106;
	t109 = -t133 * t151 + t149;
	t93 = qJD(1) * t109 - t110 * t128;
	t158 = t105 * t106 * t93;
	t99 = t104 * t106 + 0.1e1;
	t159 = (-t104 * t158 - t110 * t157) / t99 ^ 2;
	t156 = t101 * t135;
	t155 = t109 * t110;
	t153 = t124 * t125;
	t147 = qJD(1) * t134;
	t146 = t125 * t162;
	t91 = (-t114 * t117 * t134 * t153 + t113 * t161) * t135;
	t123 = t132 * t124;
	t102 = t100 * t101;
	t97 = 0.1e1 / t99;
	t96 = t130 * t124 * t101 + 0.1e1;
	t90 = qJD(1) * t91;
	t87 = -0.2e1 * t159 + 0.2e1 * (-t97 * t157 + (-t106 * t159 - t158 * t97) * t110) * t110;
	t1 = [(-t117 * t125 * t132 - 0.2e1 * t123 * t146) * t147, 0, 0, 0, 0; (0.2e1 * (t100 * t134 + t156 * t91) / t96 ^ 2 * (-t102 * t130 * t90 - t147 * t156) * t124 + ((0.2e1 * t102 * t135 * t91 + t101 * t134) * t90 + (-t135 * t100 + ((t91 + (t123 * t162 + t161) * t135 * t113) * t134 - (0.2e1 * t129 * t124 ^ 2 * t146 + (t154 + (t129 - 0.2e1 * t130) * t117) * t153) * t135 * t114) * t101) * qJD(1)) / t96) * t132, 0, 0, 0, 0; 0.2e1 * (t105 * t144 + t106 * t155) * t159 + ((-qJD(1) * t110 + t109 * t128) * t105 + 0.2e1 * t155 * t158 + (t144 * t93 - (-qJD(1) * t111 + t128 * t144) * t110 + t109 * t92) * t106) * t97, 0, t87, t87, 0;];
	JaD_rot = t1;
end