% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
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
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(10);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [-t10, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t11, 0, 0, 0, 0, 0; -t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t13 = qJ(1) + pkin(10);
	t11 = sin(t13);
	t14 = sin(qJ(3));
	t19 = t11 * t14;
	t15 = cos(qJ(3));
	t18 = t11 * t15;
	t12 = cos(t13);
	t17 = t12 * t14;
	t16 = t12 * t15;
	t1 = [-t18, 0, -t17, 0, 0, 0; t16, 0, -t19, 0, 0, 0; 0, 0, t15, 0, 0, 0; t19, 0, -t16, 0, 0, 0; -t17, 0, -t18, 0, 0, 0; 0, 0, -t14, 0, 0, 0; t12, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (42->14), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t26 = qJ(1) + pkin(10);
	t22 = sin(t26);
	t27 = qJ(3) + qJ(4);
	t24 = sin(t27);
	t31 = t22 * t24;
	t25 = cos(t27);
	t30 = t22 * t25;
	t23 = cos(t26);
	t29 = t23 * t24;
	t28 = t23 * t25;
	t1 = [-t30, 0, -t29, -t29, 0, 0; t28, 0, -t31, -t31, 0, 0; 0, 0, t25, t25, 0, 0; t31, 0, -t28, -t28, 0, 0; -t29, 0, -t30, -t30, 0, 0; 0, 0, -t24, -t24, 0, 0; t23, 0, 0, 0, 0, 0; t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:31
	% EndTime: 2019-10-10 01:23:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (58->14), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t31 = qJ(3) + qJ(4) + pkin(11);
	t27 = sin(t31);
	t32 = qJ(1) + pkin(10);
	t29 = sin(t32);
	t36 = t29 * t27;
	t28 = cos(t31);
	t35 = t29 * t28;
	t30 = cos(t32);
	t34 = t30 * t27;
	t33 = t30 * t28;
	t1 = [-t35, 0, -t34, -t34, 0, 0; t33, 0, -t36, -t36, 0, 0; 0, 0, t28, t28, 0, 0; t36, 0, -t33, -t33, 0, 0; -t34, 0, -t35, -t35, 0, 0; 0, 0, -t27, -t27, 0, 0; t30, 0, 0, 0, 0, 0; t29, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:23:32
	% EndTime: 2019-10-10 01:23:32
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (107->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->25)
	t111 = qJ(3) + qJ(4) + pkin(11);
	t108 = cos(t111);
	t113 = sin(qJ(6));
	t121 = t108 * t113;
	t112 = qJ(1) + pkin(10);
	t109 = sin(t112);
	t120 = t109 * t113;
	t114 = cos(qJ(6));
	t119 = t109 * t114;
	t110 = cos(t112);
	t118 = t110 * t113;
	t117 = t110 * t114;
	t107 = sin(t111);
	t116 = t107 * t119;
	t115 = t107 * t117;
	t106 = t108 * t114;
	t105 = t110 * t108;
	t104 = t109 * t108;
	t103 = t107 * t118;
	t102 = t107 * t120;
	t101 = t108 * t117 + t120;
	t100 = -t108 * t118 + t119;
	t99 = -t108 * t119 + t118;
	t98 = t108 * t120 + t117;
	t1 = [t99, 0, -t115, -t115, 0, t100; t101, 0, -t116, -t116, 0, -t98; 0, 0, t106, t106, 0, -t107 * t113; t98, 0, t103, t103, 0, -t101; t100, 0, t102, t102, 0, t99; 0, 0, -t121, -t121, 0, -t107 * t114; -t109 * t107, 0, t105, t105, 0, 0; t110 * t107, 0, t104, t104, 0, 0; 0, 0, t107, t107, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end