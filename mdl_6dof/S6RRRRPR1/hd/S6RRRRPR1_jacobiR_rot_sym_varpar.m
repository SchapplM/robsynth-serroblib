% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:06
	% EndTime: 2019-10-10 12:35:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:06
	% EndTime: 2019-10-10 12:35:06
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
	% StartTime: 2019-10-10 12:35:06
	% EndTime: 2019-10-10 12:35:07
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
	% StartTime: 2019-10-10 12:35:06
	% EndTime: 2019-10-10 12:35:07
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
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (61->18), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t28 = qJ(2) + qJ(3) + qJ(4);
	t26 = sin(t28);
	t29 = sin(qJ(1));
	t34 = t29 * t26;
	t27 = cos(t28);
	t33 = t29 * t27;
	t30 = cos(qJ(1));
	t32 = t30 * t26;
	t31 = t30 * t27;
	t1 = [-t33, -t32, -t32, -t32, 0, 0; t31, -t34, -t34, -t34, 0, 0; 0, t27, t27, t27, 0, 0; t34, -t31, -t31, -t31, 0, 0; -t32, -t33, -t33, -t33, 0, 0; 0, -t26, -t26, -t26, 0, 0; t30, 0, 0, 0, 0, 0; t29, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (83->18), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t33 = qJ(2) + qJ(3) + qJ(4) + pkin(11);
	t31 = sin(t33);
	t34 = sin(qJ(1));
	t39 = t34 * t31;
	t32 = cos(t33);
	t38 = t34 * t32;
	t35 = cos(qJ(1));
	t37 = t35 * t31;
	t36 = t35 * t32;
	t1 = [-t38, -t37, -t37, -t37, 0, 0; t36, -t39, -t39, -t39, 0, 0; 0, t32, t32, t32, 0, 0; t39, -t36, -t36, -t36, 0, 0; -t37, -t38, -t38, -t38, 0, 0; 0, -t31, -t31, -t31, 0, 0; t35, 0, 0, 0, 0, 0; t34, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:35:07
	% EndTime: 2019-10-10 12:35:07
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (137->19), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
	t112 = qJ(2) + qJ(3) + qJ(4) + pkin(11);
	t111 = cos(t112);
	t113 = sin(qJ(6));
	t123 = t111 * t113;
	t114 = sin(qJ(1));
	t122 = t114 * t113;
	t115 = cos(qJ(6));
	t121 = t114 * t115;
	t116 = cos(qJ(1));
	t120 = t116 * t113;
	t119 = t116 * t115;
	t110 = sin(t112);
	t118 = t110 * t121;
	t117 = t110 * t119;
	t109 = t116 * t111;
	t108 = t111 * t115;
	t107 = t114 * t111;
	t106 = t110 * t120;
	t105 = t110 * t122;
	t104 = t111 * t119 + t122;
	t103 = -t111 * t120 + t121;
	t102 = -t111 * t121 + t120;
	t101 = t111 * t122 + t119;
	t1 = [t102, -t117, -t117, -t117, 0, t103; t104, -t118, -t118, -t118, 0, -t101; 0, t108, t108, t108, 0, -t110 * t113; t101, t106, t106, t106, 0, -t104; t103, t105, t105, t105, 0, t102; 0, -t123, -t123, -t123, 0, -t110 * t115; -t114 * t110, t109, t109, t109, 0, 0; t116 * t110, t107, t107, t107, 0, 0; 0, t110, t110, t110, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end