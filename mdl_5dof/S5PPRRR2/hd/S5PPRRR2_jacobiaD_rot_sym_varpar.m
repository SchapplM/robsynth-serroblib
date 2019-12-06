% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPRRR2
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
%   Wie in S5PPRRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PPRRR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:15:17
	% EndTime: 2019-12-05 15:15:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:15:17
	% EndTime: 2019-12-05 15:15:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:15:17
	% EndTime: 2019-12-05 15:15:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:15:17
	% EndTime: 2019-12-05 15:15:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:15:18
	% EndTime: 2019-12-05 15:15:18
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (973->43), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->60)
	t104 = sin(pkin(8));
	t103 = pkin(9) + qJ(3);
	t100 = cos(t103);
	t99 = sin(t103);
	t95 = t99 ^ 2;
	t130 = t95 / t100 ^ 2;
	t120 = 0.1e1 + t130;
	t143 = t104 * t120;
	t142 = qJD(3) * t99;
	t106 = sin(qJ(4));
	t107 = cos(qJ(4));
	t105 = cos(pkin(8));
	t125 = t100 * t105;
	t84 = t104 * t106 + t107 * t125;
	t128 = t104 * t99;
	t87 = atan2(-t128, -t100);
	t85 = sin(t87);
	t86 = cos(t87);
	t70 = -t86 * t100 - t85 * t128;
	t67 = 0.1e1 / t70;
	t80 = 0.1e1 / t84;
	t101 = t104 ^ 2;
	t90 = t101 * t130 + 0.1e1;
	t88 = 0.1e1 / t90;
	t72 = t88 * t143;
	t117 = t104 * t72 - 0.1e1;
	t127 = t104 - t72;
	t129 = t100 * t85;
	t131 = t86 * t99;
	t61 = -t117 * t131 - t127 * t129;
	t140 = 0.2e1 * t61;
	t68 = 0.1e1 / t70 ^ 2;
	t81 = 0.1e1 / t84 ^ 2;
	t102 = t105 ^ 2;
	t123 = qJD(3) * t100;
	t135 = t68 * t99;
	t121 = t86 * t128;
	t71 = qJD(3) * t72;
	t60 = (-t121 + t129) * t71 + (-t104 * t129 + t131) * qJD(3);
	t138 = t60 * t67 * t68;
	t66 = t102 * t95 * t68 + 0.1e1;
	t139 = (t123 * t135 - t95 * t138) * t102 / t66 ^ 2;
	t83 = -t104 * t107 + t106 * t125;
	t132 = t81 * t83;
	t119 = t105 * t142;
	t126 = qJD(4) * t83;
	t77 = -t107 * t119 - t126;
	t133 = t77 * t80 * t81;
	t79 = t83 ^ 2;
	t75 = t79 * t81 + 0.1e1;
	t76 = -t84 * qJD(4) + t106 * t119;
	t137 = (-t76 * t132 - t79 * t133) / t75 ^ 2;
	t64 = 0.1e1 / t66;
	t136 = t64 * t68;
	t134 = t76 * t81;
	t122 = -0.2e1 * t137;
	t116 = -t106 * t80 + t107 * t132;
	t73 = 0.1e1 / t75;
	t62 = 0.2e1 * (t88 - t120 / t90 ^ 2 * t101) / t100 * t142 * t143;
	t1 = [0, 0, t62, 0, 0; 0, 0, ((-0.2e1 * t67 * t139 + (-qJD(3) * t61 - t60) * t136) * t100 + (t68 * t139 * t140 + (t62 * t68 * t121 + t138 * t140 - qJD(3) * t67 - (t127 * qJD(3) + t117 * t71) * t85 * t135) * t64 - (t62 * t85 + (qJD(3) + (-0.2e1 * t104 + t72) * t71) * t86) * t100 * t136) * t99) * t105, 0, 0; 0, 0, (t116 * t99 * t122 + (t116 * t123 + ((t77 - t126) * t81 * t106 + (-qJD(4) * t80 - 0.2e1 * t83 * t133 - t134) * t107) * t99) * t73) * t105, t122 + 0.2e1 * (-t73 * t134 + (-t73 * t133 - t81 * t137) * t83) * t83, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:15:18
	% EndTime: 2019-12-05 15:15:18
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (1374->45), mult. (1346->104), div. (247->12), fcn. (1511->9), ass. (0->63)
	t137 = sin(pkin(8));
	t134 = pkin(9) + qJ(3);
	t128 = sin(t134);
	t124 = t128 ^ 2;
	t129 = cos(t134);
	t162 = t124 / t129 ^ 2;
	t150 = 0.1e1 + t162;
	t174 = t137 * t150;
	t173 = qJD(3) * t128;
	t136 = qJ(4) + qJ(5);
	t130 = sin(t136);
	t131 = cos(t136);
	t138 = cos(pkin(8));
	t161 = t129 * t138;
	t113 = t137 * t130 + t131 * t161;
	t160 = t137 * t128;
	t116 = atan2(-t160, -t129);
	t114 = sin(t116);
	t115 = cos(t116);
	t99 = -t114 * t160 - t115 * t129;
	t96 = 0.1e1 / t99;
	t132 = t137 ^ 2;
	t120 = t132 * t162 + 0.1e1;
	t118 = 0.1e1 / t120;
	t101 = t118 * t174;
	t163 = t114 * t129;
	t147 = t115 * t128 - t137 * t163;
	t155 = t115 * t160;
	t148 = -t155 + t163;
	t90 = t148 * t101 + t147;
	t171 = 0.2e1 * t90;
	t109 = 0.1e1 / t113;
	t110 = 0.1e1 / t113 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t133 = t138 ^ 2;
	t157 = qJD(3) * t129;
	t166 = t128 * t97;
	t100 = qJD(3) * t101;
	t89 = t147 * qJD(3) + t148 * t100;
	t169 = t89 * t96 * t97;
	t95 = t133 * t124 * t97 + 0.1e1;
	t170 = (-t124 * t169 + t157 * t166) * t133 / t95 ^ 2;
	t93 = 0.1e1 / t95;
	t168 = t93 * t97;
	t154 = t130 * t161;
	t112 = -t137 * t131 + t154;
	t108 = t112 ^ 2;
	t104 = t108 * t110 + 0.1e1;
	t135 = qJD(4) + qJD(5);
	t152 = t138 * t173;
	t105 = -t113 * t135 + t130 * t152;
	t164 = t110 * t112;
	t106 = -t135 * t154 + (t135 * t137 - t152) * t131;
	t165 = t106 * t109 * t110;
	t167 = 0.1e1 / t104 ^ 2 * (-t105 * t164 - t108 * t165);
	t158 = t101 - t137;
	t156 = -0.2e1 * t167;
	t151 = t101 * t137 - 0.1e1;
	t149 = -t109 * t130 + t131 * t164;
	t102 = 0.1e1 / t104;
	t92 = 0.2e1 * (t118 - t150 / t120 ^ 2 * t132) / t129 * t173 * t174;
	t87 = t156 + 0.2e1 * (-t102 * t105 * t110 + (-t102 * t165 - t110 * t167) * t112) * t112;
	t1 = [0, 0, t92, 0, 0; 0, 0, ((-0.2e1 * t96 * t170 + (-qJD(3) * t90 - t89) * t168) * t129 + (t97 * t170 * t171 + (t92 * t97 * t155 + t169 * t171 - qJD(3) * t96 - (-t158 * qJD(3) + t151 * t100) * t114 * t166) * t93 - (t114 * t92 + (-t151 * qJD(3) + t158 * t100) * t115) * t129 * t168) * t128) * t138, 0, 0; 0, 0, (t149 * t128 * t156 + (t149 * t157 + ((-t109 * t135 - 0.2e1 * t112 * t165) * t131 + (-t105 * t131 + (-t112 * t135 + t106) * t130) * t110) * t128) * t102) * t138, t87, t87;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end