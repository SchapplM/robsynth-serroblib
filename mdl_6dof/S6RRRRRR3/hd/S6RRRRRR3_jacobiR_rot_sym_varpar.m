% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:20:13
	% EndTime: 2019-10-10 13:20:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:20:13
	% EndTime: 2019-10-10 13:20:13
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
	% StartTime: 2019-10-10 13:20:13
	% EndTime: 2019-10-10 13:20:13
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
	% StartTime: 2019-10-10 13:20:13
	% EndTime: 2019-10-10 13:20:13
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
	% StartTime: 2019-10-10 13:20:13
	% EndTime: 2019-10-10 13:20:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (47->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t101 = sin(qJ(4));
	t100 = qJ(2) + qJ(3);
	t99 = cos(t100);
	t111 = t99 * t101;
	t102 = sin(qJ(1));
	t110 = t102 * t101;
	t103 = cos(qJ(4));
	t109 = t102 * t103;
	t104 = cos(qJ(1));
	t108 = t104 * t101;
	t107 = t104 * t103;
	t98 = sin(t100);
	t106 = t98 * t109;
	t105 = t98 * t107;
	t97 = t104 * t99;
	t96 = t99 * t103;
	t95 = t102 * t99;
	t94 = t98 * t108;
	t93 = t98 * t110;
	t92 = t99 * t107 + t110;
	t91 = -t99 * t108 + t109;
	t90 = -t99 * t109 + t108;
	t89 = t99 * t110 + t107;
	t1 = [t90, -t105, -t105, t91, 0, 0; t92, -t106, -t106, -t89, 0, 0; 0, t96, t96, -t98 * t101, 0, 0; t89, t94, t94, -t92, 0, 0; t91, t93, t93, t90, 0, 0; 0, -t111, -t111, -t98 * t103, 0, 0; -t102 * t98, t97, t97, 0, 0, 0; t104 * t98, t95, t95, 0, 0, 0; 0, t98, t98, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:20:13
	% EndTime: 2019-10-10 13:20:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (99->21), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
	t115 = qJ(4) + qJ(5);
	t111 = sin(t115);
	t116 = qJ(2) + qJ(3);
	t112 = sin(t116);
	t127 = t112 * t111;
	t113 = cos(t115);
	t126 = t112 * t113;
	t114 = cos(t116);
	t125 = t114 * t111;
	t117 = sin(qJ(1));
	t124 = t117 * t112;
	t123 = t117 * t113;
	t109 = t117 * t114;
	t118 = cos(qJ(1));
	t122 = t118 * t112;
	t121 = t118 * t113;
	t110 = t118 * t114;
	t120 = t112 * t123;
	t119 = t112 * t121;
	t108 = t114 * t113;
	t107 = t111 * t122;
	t106 = t111 * t124;
	t105 = t113 * t110 + t117 * t111;
	t104 = -t111 * t110 + t123;
	t103 = -t113 * t109 + t118 * t111;
	t102 = t111 * t109 + t121;
	t1 = [t103, -t119, -t119, t104, t104, 0; t105, -t120, -t120, -t102, -t102, 0; 0, t108, t108, -t127, -t127, 0; t102, t107, t107, -t105, -t105, 0; t104, t106, t106, t103, t103, 0; 0, -t125, -t125, -t126, -t126, 0; -t124, t110, t110, 0, 0, 0; t122, t109, t109, 0, 0, 0; 0, t112, t112, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:20:13
	% EndTime: 2019-10-10 13:20:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (171->25), mult. (80->20), div. (0->0), fcn. (138->6), ass. (0->25)
	t120 = qJ(4) + qJ(5) + qJ(6);
	t116 = sin(t120);
	t121 = qJ(2) + qJ(3);
	t118 = sin(t121);
	t130 = t118 * t116;
	t117 = cos(t120);
	t129 = t118 * t117;
	t119 = cos(t121);
	t128 = t119 * t116;
	t122 = sin(qJ(1));
	t127 = t122 * t118;
	t114 = t122 * t119;
	t123 = cos(qJ(1));
	t126 = t123 * t118;
	t115 = t123 * t119;
	t125 = t117 * t127;
	t124 = t117 * t126;
	t113 = t119 * t117;
	t112 = t116 * t126;
	t111 = t116 * t127;
	t110 = t117 * t115 + t122 * t116;
	t109 = -t116 * t115 + t122 * t117;
	t108 = -t117 * t114 + t123 * t116;
	t107 = t116 * t114 + t123 * t117;
	t1 = [t108, -t124, -t124, t109, t109, t109; t110, -t125, -t125, -t107, -t107, -t107; 0, t113, t113, -t130, -t130, -t130; t107, t112, t112, -t110, -t110, -t110; t109, t111, t111, t108, t108, t108; 0, -t128, -t128, -t129, -t129, -t129; -t127, t115, t115, 0, 0, 0; t126, t114, t114, 0, 0, 0; 0, t118, t118, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end