% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPRP6
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
%   Wie in S5RPPRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPPRP6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:07:55
	% EndTime: 2019-12-29 16:07:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:07:50
	% EndTime: 2019-12-29 16:07:50
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:07:49
	% EndTime: 2019-12-29 16:07:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:07:50
	% EndTime: 2019-12-29 16:07:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:07:50
	% EndTime: 2019-12-29 16:07:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:07:50
	% EndTime: 2019-12-29 16:07:51
	% DurationCPUTime: 1.17s
	% Computational Cost: add. (1704->70), mult. (1835->162), div. (470->13), fcn. (2177->7), ass. (0->72)
	t90 = cos(qJ(1));
	t117 = qJD(1) * t90;
	t82 = pkin(7) + qJ(4);
	t80 = sin(t82);
	t75 = 0.1e1 / t80;
	t138 = 0.2e1 * t75;
	t137 = 0.2e1 * t90;
	t114 = qJD(4) * t90;
	t89 = sin(qJ(1));
	t118 = qJD(1) * t89;
	t81 = cos(t82);
	t108 = t81 * t118;
	t76 = 0.1e1 / t80 ^ 2;
	t79 = t81 ^ 2;
	t122 = t76 * t79;
	t88 = t90 ^ 2;
	t71 = t88 * t122 + 0.1e1;
	t69 = 0.1e1 / t71;
	t55 = ((t80 * t114 + t108) * t75 + t114 * t122) * t69;
	t135 = -t55 + t114;
	t119 = t90 * t81;
	t68 = atan2(-t119, t80);
	t66 = sin(t68);
	t67 = cos(t68);
	t62 = -t66 * t119 + t67 * t80;
	t59 = 0.1e1 / t62;
	t84 = 0.1e1 / t89;
	t60 = 0.1e1 / t62 ^ 2;
	t134 = t69 - 0.1e1;
	t102 = t79 * t89 * t117;
	t115 = qJD(4) * t81;
	t83 = t89 ^ 2;
	t121 = t79 * t83;
	t125 = t67 * t81;
	t51 = (-t55 * t90 + qJD(4)) * t125 + (t135 * t80 + t108) * t66;
	t132 = t51 * t59 * t60;
	t58 = t60 * t121 + 0.1e1;
	t133 = (-t121 * t132 + (-t80 * t83 * t115 + t102) * t60) / t58 ^ 2;
	t131 = t55 * t81;
	t130 = t60 * t81;
	t129 = t60 * t89;
	t123 = t75 * t81;
	t77 = t75 * t76;
	t78 = t81 * t79;
	t98 = qJD(4) * (-t77 * t78 - t123);
	t128 = (-t76 * t102 + t88 * t98) / t71 ^ 2;
	t120 = 0.1e1 / t89 ^ 2 * t88;
	t74 = t76 * t120 + 0.1e1;
	t99 = (-0.1e1 - 0.1e1 / t83 * t88) * t84 * t117;
	t127 = (-t77 * t115 * t120 + t76 * t99) / t74 ^ 2;
	t126 = t66 * t90;
	t124 = t75 * t79;
	t116 = qJD(4) * t80;
	t113 = -0.2e1 * t132;
	t112 = t81 * t133;
	t111 = t81 * t129;
	t110 = t81 * t128;
	t109 = t69 * t124;
	t107 = 0.1e1 + t122;
	t106 = -0.1e1 - t120;
	t105 = t128 * t137;
	t104 = t90 * t109;
	t103 = t134 * t81 * t66;
	t101 = t107 * t89;
	t100 = t106 * t81 * t76;
	t72 = 0.1e1 / t74;
	t64 = t107 * t90 * t69;
	t56 = 0.1e1 / t58;
	t54 = (-t67 * t104 - t103) * t89;
	t53 = t80 * t126 + t125 + (-t67 * t119 - t66 * t80) * t64;
	t52 = -t107 * t105 + (-qJD(1) * t101 + t98 * t137) * t69;
	t1 = [-0.2e1 * t89 * t75 * t110 + (-qJD(4) * t101 + t117 * t123) * t69, 0, 0, t52, 0; (0.2e1 * t59 * t112 + (t59 * t116 + (qJD(1) * t54 + t51) * t130) * t56) * t90 + (-0.2e1 * t60 * t112 * t54 + (((t55 * t104 + t134 * t116 + 0.2e1 * t110) * t66 + (t105 * t124 + t131 + (-t131 + (t76 * t78 + 0.2e1 * t81) * t114) * t69) * t67) * t111 + (t81 * t113 - t60 * t116) * t54 + (t59 + ((t83 - t88) * t67 * t109 - t90 * t103) * t60) * t81 * qJD(1)) * t56) * t89, 0, 0, 0.2e1 * (-t53 * t130 - t59 * t80) * t89 * t133 + ((t59 * t117 + (-qJD(4) * t53 - t51) * t129) * t80 + (t89 * qJD(4) * t59 + (-t52 * t67 * t90 + t135 * t66 + (-qJD(4) * t66 + t118 * t67 + t126 * t55) * t64) * t111 + (t89 * t113 + t60 * t117) * t53 + ((-t52 - t118) * t66 + ((t64 * t90 - 0.1e1) * qJD(4) + (-t64 + t90) * t55) * t67) * t80 * t129) * t81) * t56, 0; t106 * t127 * t138 + (qJD(4) * t100 + t99 * t138) * t72, 0, 0, -0.2e1 * t84 * t76 * t119 * t127 + ((-0.2e1 * t77 * t79 - t75) * t84 * t114 + qJD(1) * t100) * t72, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end