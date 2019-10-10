% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:54
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:01
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
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:02
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
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t17 = qJ(2) + pkin(11);
	t15 = sin(t17);
	t18 = sin(qJ(1));
	t23 = t18 * t15;
	t16 = cos(t17);
	t22 = t18 * t16;
	t19 = cos(qJ(1));
	t21 = t19 * t15;
	t20 = t19 * t16;
	t1 = [-t22, -t21, 0, 0, 0, 0; t20, -t23, 0, 0, 0, 0; 0, t16, 0, 0, 0, 0; t23, -t20, 0, 0, 0, 0; -t21, -t22, 0, 0, 0, 0; 0, -t15, 0, 0, 0, 0; t19, 0, 0, 0, 0, 0; t18, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:01
	% EndTime: 2019-10-10 10:54:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (44->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t25 = qJ(2) + pkin(11) + qJ(4);
	t23 = sin(t25);
	t26 = sin(qJ(1));
	t31 = t26 * t23;
	t24 = cos(t25);
	t30 = t26 * t24;
	t27 = cos(qJ(1));
	t29 = t27 * t23;
	t28 = t27 * t24;
	t1 = [-t30, -t29, 0, -t29, 0, 0; t28, -t31, 0, -t31, 0, 0; 0, t24, 0, t24, 0, 0; t31, -t28, 0, -t28, 0, 0; -t29, -t30, 0, -t30, 0, 0; 0, -t23, 0, -t23, 0, 0; t27, 0, 0, 0, 0, 0; t26, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:02
	% EndTime: 2019-10-10 10:54:02
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (77->16), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
	t103 = qJ(2) + pkin(11) + qJ(4);
	t102 = cos(t103);
	t104 = sin(qJ(5));
	t114 = t102 * t104;
	t105 = sin(qJ(1));
	t113 = t105 * t104;
	t106 = cos(qJ(5));
	t112 = t105 * t106;
	t107 = cos(qJ(1));
	t111 = t107 * t104;
	t110 = t107 * t106;
	t101 = sin(t103);
	t109 = t101 * t112;
	t108 = t101 * t110;
	t100 = t107 * t102;
	t99 = t102 * t106;
	t98 = t105 * t102;
	t97 = t101 * t111;
	t96 = t101 * t113;
	t95 = t102 * t110 + t113;
	t94 = -t102 * t111 + t112;
	t93 = -t102 * t112 + t111;
	t92 = t102 * t113 + t110;
	t1 = [t93, -t108, 0, -t108, t94, 0; t95, -t109, 0, -t109, -t92, 0; 0, t99, 0, t99, -t101 * t104, 0; t92, t97, 0, t97, -t95, 0; t94, t96, 0, t96, t93, 0; 0, -t114, 0, -t114, -t101 * t106, 0; -t105 * t101, t100, 0, t100, 0, 0; t107 * t101, t98, 0, t98, 0, 0; 0, t101, 0, t101, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:54:02
	% EndTime: 2019-10-10 10:54:02
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (135->21), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
	t116 = qJ(2) + pkin(11) + qJ(4);
	t114 = sin(t116);
	t119 = qJ(5) + qJ(6);
	t117 = sin(t119);
	t130 = t114 * t117;
	t118 = cos(t119);
	t129 = t114 * t118;
	t115 = cos(t116);
	t128 = t115 * t117;
	t120 = sin(qJ(1));
	t127 = t120 * t117;
	t126 = t120 * t118;
	t121 = cos(qJ(1));
	t125 = t121 * t117;
	t124 = t121 * t118;
	t123 = t114 * t126;
	t122 = t114 * t124;
	t113 = t121 * t115;
	t112 = t120 * t115;
	t111 = t115 * t118;
	t110 = t114 * t125;
	t109 = t114 * t127;
	t108 = t115 * t124 + t127;
	t107 = -t115 * t125 + t126;
	t106 = -t115 * t126 + t125;
	t105 = t115 * t127 + t124;
	t1 = [t106, -t122, 0, -t122, t107, t107; t108, -t123, 0, -t123, -t105, -t105; 0, t111, 0, t111, -t130, -t130; t105, t110, 0, t110, -t108, -t108; t107, t109, 0, t109, t106, t106; 0, -t128, 0, -t128, -t129, -t129; -t120 * t114, t113, 0, t113, 0, 0; t121 * t114, t112, 0, t112, 0, 0; 0, t114, 0, t114, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end