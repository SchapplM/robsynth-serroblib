% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPP4
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
%   Wie in S5PRRPP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:37
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRPP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:37:36
	% EndTime: 2019-12-29 15:37:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:37:36
	% EndTime: 2019-12-29 15:37:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:37:36
	% EndTime: 2019-12-29 15:37:36
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:37:36
	% EndTime: 2019-12-29 15:37:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:37:31
	% EndTime: 2019-12-29 15:37:32
	% DurationCPUTime: 1.13s
	% Computational Cost: add. (1458->72), mult. (1835->161), div. (470->13), fcn. (2177->7), ass. (0->74)
	t83 = pkin(7) + qJ(2);
	t81 = sin(t83);
	t140 = 0.2e1 * t81;
	t82 = cos(t83);
	t78 = 0.1e1 / t82;
	t125 = t78 * t81;
	t76 = t81 ^ 2;
	t77 = t82 ^ 2;
	t139 = qJD(2) * (t76 / t77 + 0.1e1) * t125;
	t89 = sin(qJ(3));
	t124 = t81 * t89;
	t90 = cos(qJ(3));
	t71 = atan2(-t124, -t90);
	t66 = sin(t71);
	t111 = t66 * t124;
	t67 = cos(t71);
	t63 = -t67 * t90 - t111;
	t60 = 0.1e1 / t63;
	t86 = 0.1e1 / t90;
	t61 = 0.1e1 / t63 ^ 2;
	t87 = 0.1e1 / t90 ^ 2;
	t138 = -0.2e1 * t89;
	t85 = t89 ^ 2;
	t122 = t85 * t87;
	t74 = t76 * t122 + 0.1e1;
	t72 = 0.1e1 / t74;
	t137 = t72 - 0.1e1;
	t118 = qJD(2) * t89;
	t109 = t82 * t118;
	t117 = qJD(3) * t81;
	t128 = t67 * t89;
	t116 = qJD(3) * t90;
	t55 = (-(-t81 * t116 - t109) * t86 + t117 * t122) * t72;
	t51 = (-t55 * t81 + qJD(3)) * t128 + (-t109 + (t55 - t117) * t90) * t66;
	t136 = t51 * t60 * t61;
	t135 = t55 * t66;
	t134 = t55 * t89;
	t133 = t61 * t82;
	t132 = t61 * t89;
	t119 = qJD(2) * t82;
	t102 = t81 * t85 * t119;
	t84 = t89 * t85;
	t88 = t86 * t87;
	t98 = qJD(3) * (t84 * t88 + t86 * t89);
	t131 = (t87 * t102 + t76 * t98) / t74 ^ 2;
	t107 = 0.1e1 + t122;
	t65 = t107 * t81 * t72;
	t130 = t65 * t81;
	t129 = t66 * t90;
	t127 = t76 / t82 ^ 2;
	t126 = t77 * t85;
	t123 = t85 * t86;
	t121 = t87 * t89;
	t120 = qJD(2) * t81;
	t58 = t61 * t126 + 0.1e1;
	t115 = 0.2e1 * (-t126 * t136 + (t77 * t89 * t116 - t102) * t61) / t58 ^ 2;
	t114 = 0.2e1 * t136;
	t70 = t87 * t127 + 0.1e1;
	t113 = 0.2e1 * (t88 * qJD(3) * t89 * t127 + t87 * t139) / t70 ^ 2;
	t112 = t82 * t132;
	t110 = t72 * t123;
	t108 = 0.1e1 + t127;
	t106 = t89 * t115;
	t105 = t131 * t140;
	t104 = t131 * t138;
	t103 = t81 * t110;
	t101 = t107 * t82;
	t99 = t108 * t121;
	t68 = 0.1e1 / t70;
	t56 = 0.1e1 / t58;
	t54 = (t137 * t89 * t66 - t67 * t103) * t82;
	t53 = -t81 * t129 + t128 + (-t67 * t124 + t129) * t65;
	t52 = -t107 * t105 + (qJD(2) * t101 + t140 * t98) * t72;
	t1 = [0, t82 * t86 * t104 + (-t118 * t81 * t86 + qJD(3) * t101) * t72, t52, 0, 0; 0, (t60 * t106 + (-t60 * t116 + (qJD(2) * t54 + t51) * t132) * t56) * t81 + (t61 * t106 * t54 + (-((t103 * t55 + t116 * t137 + t104) * t66 + (t105 * t123 - t134 + (t134 + (-t84 * t87 + t138) * t117) * t72) * t67) * t112 + (t114 * t89 - t116 * t61) * t54 + (-t60 + ((-t76 + t77) * t67 * t110 + t137 * t111) * t61) * t118) * t56) * t82, (t132 * t53 - t60 * t90) * t82 * t115 + ((-t60 * t120 + (-qJD(3) * t53 - t51) * t133) * t90 + (-t82 * qJD(3) * t60 - (-t52 * t67 * t81 + t66 * t117 + t130 * t135 - t135 + (-qJD(3) * t66 - t119 * t67) * t65) * t112 + (t114 * t82 + t120 * t61) * t53 - ((t52 - t119) * t66 + ((0.1e1 - t130) * qJD(3) + (t65 - t81) * t55) * t67) * t90 * t133) * t89) * t56, 0, 0; 0, t108 * t86 * t113 + (-qJD(3) * t99 - 0.2e1 * t139 * t86) * t68, t113 * t121 * t125 + ((-0.2e1 * t85 * t88 - t86) * t78 * t117 - qJD(2) * t99) * t68, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:37:36
	% EndTime: 2019-12-29 15:37:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end