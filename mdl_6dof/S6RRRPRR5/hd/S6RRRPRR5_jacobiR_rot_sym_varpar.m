% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:00
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0, 0; t13, -t16, 0, 0, 0, 0; 0, t11, 0, 0, 0, 0; t16, -t13, 0, 0, 0, 0; -t15, -t14, 0, 0, 0, 0; 0, -t9, 0, 0, 0, 0; t12, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t22 = qJ(2) + qJ(3);
	t20 = sin(t22);
	t23 = sin(qJ(1));
	t28 = t23 * t20;
	t21 = cos(t22);
	t27 = t23 * t21;
	t24 = cos(qJ(1));
	t26 = t24 * t20;
	t25 = t24 * t21;
	t1 = [-t27, -t26, -t26, 0, 0, 0; t25, -t28, -t28, 0, 0, 0; 0, t21, t21, 0, 0, 0; t28, -t25, -t25, 0, 0, 0; -t26, -t27, -t27, 0, 0, 0; 0, -t20, -t20, 0, 0, 0; t24, 0, 0, 0, 0, 0; t23, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (20->5), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t76 = cos(qJ(1));
	t75 = sin(qJ(1));
	t74 = qJ(2) + qJ(3);
	t73 = cos(t74);
	t72 = sin(t74);
	t71 = t76 * t73;
	t70 = t76 * t72;
	t69 = t75 * t73;
	t68 = t75 * t72;
	t1 = [t76, 0, 0, 0, 0, 0; t75, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t69, t70, t70, 0, 0, 0; -t71, t68, t68, 0, 0, 0; 0, -t73, -t73, 0, 0, 0; -t68, t71, t71, 0, 0, 0; t70, t69, t69, 0, 0, 0; 0, t72, t72, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (44->13), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t103 = qJ(2) + qJ(3);
	t101 = sin(t103);
	t105 = sin(qJ(1));
	t113 = t105 * t101;
	t104 = sin(qJ(5));
	t112 = t105 * t104;
	t106 = cos(qJ(5));
	t111 = t105 * t106;
	t107 = cos(qJ(1));
	t110 = t107 * t101;
	t109 = t107 * t104;
	t108 = t107 * t106;
	t102 = cos(t103);
	t100 = t101 * t106;
	t99 = t101 * t104;
	t98 = t102 * t108;
	t97 = t102 * t109;
	t96 = t102 * t111;
	t95 = t102 * t112;
	t94 = -t101 * t112 + t108;
	t93 = t101 * t111 + t109;
	t92 = t101 * t109 + t111;
	t91 = t101 * t108 - t112;
	t1 = [t94, t97, t97, 0, t91, 0; t92, t95, t95, 0, t93, 0; 0, t99, t99, 0, -t102 * t106, 0; -t93, t98, t98, 0, -t92, 0; t91, t96, t96, 0, t94, 0; 0, t100, t100, 0, t102 * t104, 0; -t105 * t102, -t110, -t110, 0, 0, 0; t107 * t102, -t113, -t113, 0, 0, 0; 0, t102, t102, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (95->16), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
	t131 = qJ(5) + qJ(6);
	t129 = cos(t131);
	t132 = qJ(2) + qJ(3);
	t130 = cos(t132);
	t141 = t130 * t129;
	t128 = sin(t132);
	t133 = sin(qJ(1));
	t140 = t133 * t128;
	t139 = t133 * t129;
	t138 = t133 * t130;
	t134 = cos(qJ(1));
	t137 = t134 * t128;
	t136 = t134 * t129;
	t135 = t134 * t130;
	t127 = sin(t131);
	t126 = t130 * t127;
	t125 = t128 * t129;
	t124 = t128 * t127;
	t123 = t129 * t135;
	t122 = t127 * t135;
	t121 = t129 * t138;
	t120 = t127 * t138;
	t119 = -t127 * t140 + t136;
	t118 = t134 * t127 + t128 * t139;
	t117 = t127 * t137 + t139;
	t116 = -t133 * t127 + t128 * t136;
	t1 = [t119, t122, t122, 0, t116, t116; t117, t120, t120, 0, t118, t118; 0, t124, t124, 0, -t141, -t141; -t118, t123, t123, 0, -t117, -t117; t116, t121, t121, 0, t119, t119; 0, t125, t125, 0, t126, t126; -t138, -t137, -t137, 0, 0, 0; t135, -t140, -t140, 0, 0, 0; 0, t130, t130, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end