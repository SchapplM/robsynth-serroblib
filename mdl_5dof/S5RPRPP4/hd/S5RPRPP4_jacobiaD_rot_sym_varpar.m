% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPP4
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
%   Wie in S5RPRPP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRPP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:40:51
	% EndTime: 2019-12-29 16:40:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:40:51
	% EndTime: 2019-12-29 16:40:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:40:51
	% EndTime: 2019-12-29 16:40:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:40:51
	% EndTime: 2019-12-29 16:40:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:40:51
	% EndTime: 2019-12-29 16:40:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:40:52
	% EndTime: 2019-12-29 16:40:53
	% DurationCPUTime: 1.16s
	% Computational Cost: add. (1704->70), mult. (1835->162), div. (470->13), fcn. (2177->7), ass. (0->72)
	t93 = cos(qJ(1));
	t120 = qJD(1) * t93;
	t85 = qJ(3) + pkin(7);
	t83 = sin(t85);
	t78 = 0.1e1 / t83;
	t141 = 0.2e1 * t78;
	t140 = 0.2e1 * t93;
	t117 = qJD(3) * t93;
	t92 = sin(qJ(1));
	t121 = qJD(1) * t92;
	t84 = cos(t85);
	t111 = t84 * t121;
	t79 = 0.1e1 / t83 ^ 2;
	t82 = t84 ^ 2;
	t125 = t79 * t82;
	t91 = t93 ^ 2;
	t74 = t91 * t125 + 0.1e1;
	t72 = 0.1e1 / t74;
	t58 = ((t83 * t117 + t111) * t78 + t117 * t125) * t72;
	t138 = -t58 + t117;
	t122 = t93 * t84;
	t71 = atan2(-t122, t83);
	t69 = sin(t71);
	t70 = cos(t71);
	t65 = -t69 * t122 + t70 * t83;
	t62 = 0.1e1 / t65;
	t87 = 0.1e1 / t92;
	t63 = 0.1e1 / t65 ^ 2;
	t137 = t72 - 0.1e1;
	t105 = t82 * t92 * t120;
	t118 = qJD(3) * t84;
	t86 = t92 ^ 2;
	t124 = t82 * t86;
	t128 = t70 * t84;
	t54 = (-t58 * t93 + qJD(3)) * t128 + (t138 * t83 + t111) * t69;
	t135 = t54 * t62 * t63;
	t61 = t63 * t124 + 0.1e1;
	t136 = (-t124 * t135 + (-t83 * t86 * t118 + t105) * t63) / t61 ^ 2;
	t134 = t58 * t84;
	t133 = t63 * t84;
	t132 = t63 * t92;
	t126 = t78 * t84;
	t80 = t78 * t79;
	t81 = t84 * t82;
	t101 = qJD(3) * (-t80 * t81 - t126);
	t131 = (t91 * t101 - t79 * t105) / t74 ^ 2;
	t102 = (-0.1e1 - 0.1e1 / t86 * t91) * t87 * t120;
	t123 = 0.1e1 / t92 ^ 2 * t91;
	t77 = t79 * t123 + 0.1e1;
	t130 = (-t80 * t118 * t123 + t79 * t102) / t77 ^ 2;
	t129 = t69 * t93;
	t127 = t78 * t82;
	t119 = qJD(3) * t83;
	t116 = -0.2e1 * t135;
	t115 = t84 * t136;
	t114 = t84 * t132;
	t113 = t84 * t131;
	t112 = t72 * t127;
	t110 = 0.1e1 + t125;
	t109 = -0.1e1 - t123;
	t108 = t131 * t140;
	t107 = t93 * t112;
	t106 = t137 * t84 * t69;
	t104 = t110 * t92;
	t103 = t109 * t84 * t79;
	t75 = 0.1e1 / t77;
	t67 = t110 * t93 * t72;
	t59 = 0.1e1 / t61;
	t57 = (-t70 * t107 - t106) * t92;
	t56 = t83 * t129 + t128 + (-t70 * t122 - t69 * t83) * t67;
	t55 = -t110 * t108 + (-qJD(1) * t104 + t101 * t140) * t72;
	t1 = [-0.2e1 * t92 * t78 * t113 + (-qJD(3) * t104 + t120 * t126) * t72, 0, t55, 0, 0; (0.2e1 * t62 * t115 + (t62 * t119 + (qJD(1) * t57 + t54) * t133) * t59) * t93 + (-0.2e1 * t63 * t115 * t57 + (((t58 * t107 + t137 * t119 + 0.2e1 * t113) * t69 + (t108 * t127 + t134 + (-t134 + (t79 * t81 + 0.2e1 * t84) * t117) * t72) * t70) * t114 + (t84 * t116 - t63 * t119) * t57 + (t62 + ((t86 - t91) * t70 * t112 - t93 * t106) * t63) * t84 * qJD(1)) * t59) * t92, 0, 0.2e1 * (-t56 * t133 - t62 * t83) * t92 * t136 + ((t62 * t120 + (-qJD(3) * t56 - t54) * t132) * t83 + (t92 * qJD(3) * t62 + (-t55 * t70 * t93 + t138 * t69 + (-qJD(3) * t69 + t121 * t70 + t129 * t58) * t67) * t114 + (t92 * t116 + t63 * t120) * t56 + ((-t55 - t121) * t69 + ((t67 * t93 - 0.1e1) * qJD(3) + (-t67 + t93) * t58) * t70) * t83 * t132) * t84) * t59, 0, 0; t109 * t130 * t141 + (qJD(3) * t103 + t102 * t141) * t75, 0, -0.2e1 * t87 * t79 * t122 * t130 + ((-0.2e1 * t80 * t82 - t78) * t87 * t117 + qJD(1) * t103) * t75, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end