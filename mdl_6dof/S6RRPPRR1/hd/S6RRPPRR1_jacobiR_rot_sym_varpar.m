% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
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
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
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
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->8), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t17 = qJ(2) + pkin(10);
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
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (14->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t47 = qJ(2) + pkin(10);
	t45 = sin(t47);
	t48 = sin(qJ(1));
	t52 = t48 * t45;
	t46 = cos(t47);
	t51 = t48 * t46;
	t49 = cos(qJ(1));
	t50 = t49 * t45;
	t44 = t49 * t46;
	t1 = [-t51, -t50, 0, 0, 0, 0; t44, -t52, 0, 0, 0, 0; 0, t46, 0, 0, 0, 0; t49, 0, 0, 0, 0, 0; t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t52, t44, 0, 0, 0, 0; t50, t51, 0, 0, 0, 0; 0, t45, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (50->13), mult. (56->8), div. (0->0), fcn. (90->6), ass. (0->14)
	t44 = qJ(2) + pkin(10);
	t42 = sin(t44);
	t43 = cos(t44);
	t45 = sin(qJ(5));
	t47 = cos(qJ(5));
	t52 = -t42 * t47 + t43 * t45;
	t49 = t42 * t45 + t43 * t47;
	t48 = cos(qJ(1));
	t46 = sin(qJ(1));
	t38 = t49 * t48;
	t37 = t52 * t48;
	t36 = t49 * t46;
	t35 = t52 * t46;
	t1 = [-t36, t37, 0, 0, -t37, 0; t38, t35, 0, 0, -t35, 0; 0, t49, 0, 0, -t49, 0; t35, t38, 0, 0, -t38, 0; -t37, t36, 0, 0, -t36, 0; 0, -t52, 0, 0, t52, 0; -t48, 0, 0, 0, 0, 0; -t46, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:35:45
	% EndTime: 2019-10-10 09:35:45
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (106->21), mult. (134->24), div. (0->0), fcn. (202->8), ass. (0->26)
	t123 = qJ(2) + pkin(10);
	t122 = cos(t123);
	t125 = sin(qJ(5));
	t141 = sin(t123);
	t142 = cos(qJ(5));
	t117 = -t122 * t125 + t141 * t142;
	t116 = t122 * t142 + t141 * t125;
	t126 = sin(qJ(1));
	t112 = t117 * t126;
	t124 = sin(qJ(6));
	t140 = t112 * t124;
	t127 = cos(qJ(6));
	t139 = t112 * t127;
	t128 = cos(qJ(1));
	t115 = t117 * t128;
	t138 = t115 * t124;
	t137 = t115 * t127;
	t136 = t116 * t124;
	t135 = t116 * t127;
	t113 = t116 * t126;
	t130 = -t113 * t127 - t128 * t124;
	t129 = t113 * t124 - t128 * t127;
	t114 = t116 * t128;
	t111 = t114 * t127 - t126 * t124;
	t110 = -t114 * t124 - t126 * t127;
	t1 = [t130, -t137, 0, 0, t137, t110; t111, -t139, 0, 0, t139, -t129; 0, t135, 0, 0, -t135, -t117 * t124; t129, t138, 0, 0, -t138, -t111; t110, t140, 0, 0, -t140, t130; 0, -t136, 0, 0, t136, -t117 * t127; t112, -t114, 0, 0, t114, 0; -t115, -t113, 0, 0, t113, 0; 0, -t117, 0, 0, t117, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end