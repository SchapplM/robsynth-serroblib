% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPPR5
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
%   Wie in S5RRPPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPPR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:30:44
	% EndTime: 2019-12-31 19:30:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:30:44
	% EndTime: 2019-12-31 19:30:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:30:44
	% EndTime: 2019-12-31 19:30:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:30:44
	% EndTime: 2019-12-31 19:30:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:30:44
	% EndTime: 2019-12-31 19:30:45
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (1893->71), mult. (1835->159), div. (470->13), fcn. (2177->7), ass. (0->71)
	t92 = sin(qJ(1));
	t122 = qJD(1) * t92;
	t141 = 0.2e1 * t92;
	t119 = qJD(2) * t92;
	t93 = cos(qJ(1));
	t121 = qJD(1) * t93;
	t85 = qJ(2) + pkin(8);
	t83 = sin(t85);
	t112 = t83 * t121;
	t79 = t83 ^ 2;
	t84 = cos(t85);
	t81 = 0.1e1 / t84 ^ 2;
	t127 = t79 * t81;
	t87 = t92 ^ 2;
	t74 = t127 * t87 + 0.1e1;
	t72 = 0.1e1 / t74;
	t80 = 0.1e1 / t84;
	t58 = (-(-t119 * t84 - t112) * t80 + t119 * t127) * t72;
	t139 = t58 - t119;
	t88 = t93 ^ 2;
	t89 = 0.1e1 / t93;
	t138 = (t87 / t88 + 0.1e1) * t89 * t122;
	t123 = t92 * t83;
	t71 = atan2(-t123, -t84);
	t69 = sin(t71);
	t114 = t69 * t123;
	t70 = cos(t71);
	t65 = -t70 * t84 - t114;
	t62 = 0.1e1 / t65;
	t63 = 0.1e1 / t65 ^ 2;
	t137 = -0.2e1 * t83;
	t136 = t72 - 0.1e1;
	t129 = t70 * t83;
	t54 = (-t58 * t92 + qJD(2)) * t129 + (t139 * t84 - t112) * t69;
	t135 = t54 * t62 * t63;
	t134 = t58 * t83;
	t133 = t63 * t83;
	t132 = t63 * t93;
	t125 = t80 * t83;
	t78 = t83 * t79;
	t82 = t80 * t81;
	t101 = qJD(2) * (t78 * t82 + t125);
	t105 = t79 * t92 * t121;
	t131 = (t101 * t87 + t105 * t81) / t74 ^ 2;
	t130 = t69 * t92;
	t128 = t79 * t80;
	t126 = t79 * t88;
	t124 = t87 / t93 ^ 2;
	t120 = qJD(2) * t84;
	t61 = t126 * t63 + 0.1e1;
	t118 = 0.2e1 * (-t126 * t135 + (t120 * t83 * t88 - t105) * t63) / t61 ^ 2;
	t117 = 0.2e1 * t135;
	t77 = t124 * t81 + 0.1e1;
	t116 = 0.2e1 * (qJD(2) * t124 * t82 * t83 + t138 * t81) / t77 ^ 2;
	t115 = t83 * t132;
	t113 = t72 * t128;
	t111 = 0.1e1 + t127;
	t110 = 0.1e1 + t124;
	t109 = t83 * t118;
	t108 = t131 * t137;
	t107 = t131 * t141;
	t106 = t92 * t113;
	t104 = t111 * t93;
	t102 = t110 * t83 * t81;
	t75 = 0.1e1 / t77;
	t67 = t111 * t92 * t72;
	t59 = 0.1e1 / t61;
	t57 = (t136 * t69 * t83 - t106 * t70) * t93;
	t56 = -t84 * t130 + t129 + (-t123 * t70 + t69 * t84) * t67;
	t55 = -t111 * t107 + (qJD(1) * t104 + t101 * t141) * t72;
	t1 = [t93 * t80 * t108 + (qJD(2) * t104 - t122 * t125) * t72, t55, 0, 0, 0; (t62 * t109 + (-t62 * t120 + (qJD(1) * t57 + t54) * t133) * t59) * t92 + (t63 * t109 * t57 + (-((t106 * t58 + t120 * t136 + t108) * t69 + (t107 * t128 - t134 + (t134 + (-t78 * t81 + t137) * t119) * t72) * t70) * t115 + (t117 * t83 - t120 * t63) * t57 + (-t62 + ((-t87 + t88) * t70 * t113 + t136 * t114) * t63) * t83 * qJD(1)) * t59) * t93, (t133 * t56 - t62 * t84) * t93 * t118 + ((-t62 * t122 + (-qJD(2) * t56 - t54) * t132) * t84 + (-t93 * qJD(2) * t62 - (-t55 * t70 * t92 - t139 * t69 + (-qJD(2) * t69 - t121 * t70 + t130 * t58) * t67) * t115 + (t117 * t93 + t122 * t63) * t56 - ((t55 - t121) * t69 + ((-t67 * t92 + 0.1e1) * qJD(2) + (t67 - t92) * t58) * t70) * t84 * t132) * t83) * t59, 0, 0, 0; t110 * t80 * t116 + (-qJD(2) * t102 - 0.2e1 * t138 * t80) * t75, t89 * t81 * t116 * t123 + ((-0.2e1 * t79 * t82 - t80) * t89 * t119 - qJD(1) * t102) * t75, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:30:44
	% EndTime: 2019-12-31 19:30:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end