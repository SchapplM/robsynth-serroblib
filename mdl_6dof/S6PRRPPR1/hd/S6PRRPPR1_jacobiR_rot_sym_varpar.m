% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(10));
	t20 = sin(pkin(6));
	t19 = sin(pkin(10));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (19->13), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t66 = sin(pkin(6));
	t69 = sin(qJ(3));
	t78 = t66 * t69;
	t70 = sin(qJ(2));
	t77 = t66 * t70;
	t71 = cos(qJ(3));
	t76 = t66 * t71;
	t72 = cos(qJ(2));
	t75 = t66 * t72;
	t68 = cos(pkin(6));
	t74 = t68 * t70;
	t73 = t68 * t72;
	t67 = cos(pkin(10));
	t65 = sin(pkin(10));
	t64 = -t65 * t74 + t67 * t72;
	t63 = -t65 * t73 - t67 * t70;
	t62 = t65 * t72 + t67 * t74;
	t61 = -t65 * t70 + t67 * t73;
	t1 = [0, t63 * t71, -t64 * t69 + t65 * t76, 0, 0, 0; 0, t61 * t71, -t62 * t69 - t67 * t76, 0, 0, 0; 0, t71 * t75, t68 * t71 - t69 * t77, 0, 0, 0; 0, -t63 * t69, -t64 * t71 - t65 * t78, 0, 0, 0; 0, -t61 * t69, -t62 * t71 + t67 * t78, 0, 0, 0; 0, -t69 * t75, -t68 * t69 - t70 * t76, 0, 0, 0; 0, t64, 0, 0, 0, 0; 0, t62, 0, 0, 0, 0; 0, t77, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (37->14), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->20)
	t75 = sin(pkin(10));
	t76 = sin(pkin(6));
	t86 = t75 * t76;
	t77 = cos(pkin(10));
	t85 = t76 * t77;
	t79 = sin(qJ(2));
	t84 = t76 * t79;
	t80 = cos(qJ(2));
	t83 = t76 * t80;
	t78 = cos(pkin(6));
	t82 = t78 * t79;
	t81 = t78 * t80;
	t74 = qJ(3) + pkin(11);
	t73 = cos(t74);
	t72 = sin(t74);
	t71 = -t75 * t82 + t77 * t80;
	t70 = -t75 * t81 - t77 * t79;
	t69 = t75 * t80 + t77 * t82;
	t68 = -t75 * t79 + t77 * t81;
	t1 = [0, t70 * t73, -t71 * t72 + t73 * t86, 0, 0, 0; 0, t68 * t73, -t69 * t72 - t73 * t85, 0, 0, 0; 0, t73 * t83, -t72 * t84 + t78 * t73, 0, 0, 0; 0, -t70 * t72, -t71 * t73 - t72 * t86, 0, 0, 0; 0, -t68 * t72, -t69 * t73 + t72 * t85, 0, 0, 0; 0, -t72 * t83, -t78 * t72 - t73 * t84, 0, 0, 0; 0, t71, 0, 0, 0, 0; 0, t69, 0, 0, 0, 0; 0, t84, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (61->20), mult. (107->52), div. (0->0), fcn. (156->10), ass. (0->27)
	t108 = qJ(3) + pkin(11);
	t107 = cos(t108);
	t109 = sin(pkin(12));
	t124 = t107 * t109;
	t112 = cos(pkin(12));
	t123 = t107 * t112;
	t116 = cos(qJ(2));
	t122 = t107 * t116;
	t110 = sin(pkin(10));
	t111 = sin(pkin(6));
	t121 = t110 * t111;
	t113 = cos(pkin(10));
	t120 = t111 * t113;
	t115 = sin(qJ(2));
	t119 = t111 * t115;
	t114 = cos(pkin(6));
	t118 = t114 * t115;
	t117 = t114 * t116;
	t106 = sin(t108);
	t105 = -t110 * t118 + t113 * t116;
	t104 = -t110 * t117 - t113 * t115;
	t103 = t110 * t116 + t113 * t118;
	t102 = -t110 * t115 + t113 * t117;
	t101 = -t106 * t119 + t114 * t107;
	t100 = -t105 * t106 + t107 * t121;
	t99 = -t103 * t106 - t107 * t120;
	t1 = [0, t104 * t123 + t105 * t109, t100 * t112, 0, 0, 0; 0, t102 * t123 + t103 * t109, t99 * t112, 0, 0, 0; 0, (t109 * t115 + t112 * t122) * t111, t101 * t112, 0, 0, 0; 0, -t104 * t124 + t105 * t112, -t100 * t109, 0, 0, 0; 0, -t102 * t124 + t103 * t112, -t99 * t109, 0, 0, 0; 0, (-t109 * t122 + t112 * t115) * t111, -t101 * t109, 0, 0, 0; 0, t104 * t106, t105 * t107 + t106 * t121, 0, 0, 0; 0, t102 * t106, t103 * t107 - t106 * t120, 0, 0, 0; 0, t111 * t116 * t106, t114 * t106 + t107 * t119, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:16
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (123->29), mult. (163->64), div. (0->0), fcn. (238->10), ass. (0->32)
	t126 = pkin(12) + qJ(6);
	t122 = sin(t126);
	t127 = qJ(3) + pkin(11);
	t125 = cos(t127);
	t142 = t122 * t125;
	t124 = cos(t126);
	t141 = t124 * t125;
	t133 = cos(qJ(2));
	t140 = t125 * t133;
	t128 = sin(pkin(10));
	t129 = sin(pkin(6));
	t139 = t128 * t129;
	t130 = cos(pkin(10));
	t138 = t129 * t130;
	t132 = sin(qJ(2));
	t137 = t129 * t132;
	t136 = t129 * t133;
	t131 = cos(pkin(6));
	t135 = t131 * t132;
	t134 = t131 * t133;
	t123 = sin(t127);
	t120 = -t128 * t135 + t130 * t133;
	t119 = t128 * t134 + t130 * t132;
	t118 = t128 * t133 + t130 * t135;
	t117 = t128 * t132 - t130 * t134;
	t116 = t123 * t131 + t125 * t137;
	t115 = -t123 * t137 + t125 * t131;
	t114 = t120 * t125 + t123 * t139;
	t113 = -t120 * t123 + t125 * t139;
	t112 = t118 * t125 - t123 * t138;
	t111 = -t118 * t123 - t125 * t138;
	t1 = [0, -t119 * t141 + t120 * t122, t113 * t124, 0, 0, -t114 * t122 + t119 * t124; 0, -t117 * t141 + t118 * t122, t111 * t124, 0, 0, -t112 * t122 + t117 * t124; 0, (t122 * t132 + t124 * t140) * t129, t115 * t124, 0, 0, -t116 * t122 - t124 * t136; 0, t119 * t142 + t120 * t124, -t113 * t122, 0, 0, -t114 * t124 - t119 * t122; 0, t117 * t142 + t118 * t124, -t111 * t122, 0, 0, -t112 * t124 - t117 * t122; 0, (-t122 * t140 + t124 * t132) * t129, -t115 * t122, 0, 0, -t116 * t124 + t122 * t136; 0, -t119 * t123, t114, 0, 0, 0; 0, -t117 * t123, t112, 0, 0, 0; 0, t123 * t136, t116, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end