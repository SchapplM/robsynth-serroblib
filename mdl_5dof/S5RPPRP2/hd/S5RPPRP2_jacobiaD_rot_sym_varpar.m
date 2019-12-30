% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPRP2
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
%   Wie in S5RPPRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:56
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPPRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:56:45
	% EndTime: 2019-12-29 15:56:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:56:45
	% EndTime: 2019-12-29 15:56:45
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:56:45
	% EndTime: 2019-12-29 15:56:45
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:56:50
	% EndTime: 2019-12-29 15:56:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:56:45
	% EndTime: 2019-12-29 15:56:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:56:45
	% EndTime: 2019-12-29 15:56:46
	% DurationCPUTime: 1.19s
	% Computational Cost: add. (2575->73), mult. (1835->159), div. (470->13), fcn. (2177->7), ass. (0->73)
	t96 = qJ(1) + pkin(7);
	t92 = sin(t96);
	t125 = qJD(1) * t92;
	t145 = 0.2e1 * t92;
	t83 = t92 ^ 2;
	t94 = cos(t96);
	t87 = t94 ^ 2;
	t88 = 0.1e1 / t94;
	t143 = (t83 / t87 + 0.1e1) * t88 * t125;
	t95 = pkin(8) + qJ(4);
	t91 = sin(t95);
	t126 = t92 * t91;
	t93 = cos(t95);
	t73 = atan2(-t126, -t93);
	t71 = sin(t73);
	t117 = t71 * t126;
	t72 = cos(t73);
	t67 = -t72 * t93 - t117;
	t64 = 0.1e1 / t67;
	t84 = 0.1e1 / t93;
	t65 = 0.1e1 / t67 ^ 2;
	t142 = -0.2e1 * t91;
	t85 = 0.1e1 / t93 ^ 2;
	t81 = t91 ^ 2;
	t130 = t81 * t85;
	t78 = t83 * t130 + 0.1e1;
	t74 = 0.1e1 / t78;
	t141 = t74 - 0.1e1;
	t124 = qJD(1) * t94;
	t115 = t91 * t124;
	t123 = qJD(4) * t92;
	t132 = t72 * t91;
	t122 = qJD(4) * t93;
	t60 = (-(-t92 * t122 - t115) * t84 + t123 * t130) * t74;
	t56 = (-t60 * t92 + qJD(4)) * t132 + (-t115 + (t60 - t123) * t93) * t71;
	t140 = t56 * t64 * t65;
	t139 = t60 * t71;
	t138 = t60 * t91;
	t137 = t65 * t91;
	t136 = t65 * t94;
	t127 = t84 * t91;
	t80 = t91 * t81;
	t86 = t84 * t85;
	t104 = qJD(4) * (t80 * t86 + t127);
	t108 = t81 * t92 * t124;
	t135 = (t83 * t104 + t85 * t108) / t78 ^ 2;
	t114 = 0.1e1 + t130;
	t70 = t114 * t92 * t74;
	t134 = t70 * t92;
	t133 = t71 * t93;
	t131 = t81 * t84;
	t129 = t81 * t87;
	t128 = t83 / t94 ^ 2;
	t63 = t65 * t129 + 0.1e1;
	t121 = 0.2e1 * (-t129 * t140 + (t87 * t91 * t122 - t108) * t65) / t63 ^ 2;
	t120 = 0.2e1 * t140;
	t79 = t85 * t128 + 0.1e1;
	t119 = 0.2e1 * (t86 * qJD(4) * t91 * t128 + t143 * t85) / t79 ^ 2;
	t118 = t91 * t136;
	t116 = t74 * t131;
	t113 = 0.1e1 + t128;
	t112 = t91 * t121;
	t111 = t135 * t142;
	t110 = t135 * t145;
	t109 = t92 * t116;
	t107 = t114 * t94;
	t105 = t113 * t91 * t85;
	t76 = 0.1e1 / t79;
	t61 = 0.1e1 / t63;
	t59 = (t141 * t91 * t71 - t72 * t109) * t94;
	t58 = -t92 * t133 + t132 + (-t72 * t126 + t133) * t70;
	t57 = -t114 * t110 + (qJD(1) * t107 + t104 * t145) * t74;
	t1 = [t94 * t84 * t111 + (qJD(4) * t107 - t125 * t127) * t74, 0, 0, t57, 0; (t64 * t112 + (-t64 * t122 + (qJD(1) * t59 + t56) * t137) * t61) * t92 + (t65 * t112 * t59 + (-((t60 * t109 + t141 * t122 + t111) * t71 + (t110 * t131 - t138 + (t138 + (-t80 * t85 + t142) * t123) * t74) * t72) * t118 + (t91 * t120 - t65 * t122) * t59 + (-t64 + ((-t83 + t87) * t72 * t116 + t141 * t117) * t65) * t91 * qJD(1)) * t61) * t94, 0, 0, (t58 * t137 - t64 * t93) * t94 * t121 + ((-t64 * t125 + (-qJD(4) * t58 - t56) * t136) * t93 + (-t94 * qJD(4) * t64 - (-t57 * t72 * t92 + t71 * t123 + t134 * t139 - t139 + (-qJD(4) * t71 - t124 * t72) * t70) * t118 + (t94 * t120 + t65 * t125) * t58 - ((t57 - t124) * t71 + ((0.1e1 - t134) * qJD(4) + (t70 - t92) * t60) * t72) * t93 * t136) * t91) * t61, 0; t113 * t84 * t119 + (-qJD(4) * t105 - 0.2e1 * t143 * t84) * t76, 0, 0, t88 * t85 * t119 * t126 + ((-0.2e1 * t81 * t86 - t84) * t88 * t123 - qJD(1) * t105) * t76, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end