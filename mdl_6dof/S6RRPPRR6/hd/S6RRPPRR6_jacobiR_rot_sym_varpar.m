% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
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
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
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
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t38 = sin(qJ(2));
	t39 = sin(qJ(1));
	t44 = t39 * t38;
	t40 = cos(qJ(2));
	t43 = t39 * t40;
	t41 = cos(qJ(1));
	t42 = t41 * t38;
	t37 = t41 * t40;
	t1 = [-t43, -t42, 0, 0, 0, 0; t37, -t44, 0, 0, 0, 0; 0, t40, 0, 0, 0, 0; t41, 0, 0, 0, 0, 0; t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t44, t37, 0, 0, 0, 0; t42, t43, 0, 0, 0, 0; 0, t38, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (12->7), mult. (36->8), div. (0->0), fcn. (58->6), ass. (0->13)
	t22 = sin(pkin(10));
	t23 = cos(pkin(10));
	t24 = sin(qJ(2));
	t26 = cos(qJ(2));
	t29 = t26 * t22 - t24 * t23;
	t28 = t24 * t22 + t26 * t23;
	t27 = cos(qJ(1));
	t25 = sin(qJ(1));
	t21 = t28 * t27;
	t20 = t29 * t27;
	t19 = t28 * t25;
	t18 = t29 * t25;
	t1 = [-t19, t20, 0, 0, 0, 0; t21, t18, 0, 0, 0, 0; 0, t28, 0, 0, 0, 0; t18, t21, 0, 0, 0, 0; -t20, t19, 0, 0, 0, 0; 0, -t29, 0, 0, 0, 0; -t27, 0, 0, 0, 0, 0; -t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (50->13), mult. (56->8), div. (0->0), fcn. (90->6), ass. (0->14)
	t41 = pkin(10) + qJ(5);
	t39 = sin(t41);
	t40 = cos(t41);
	t42 = sin(qJ(2));
	t44 = cos(qJ(2));
	t47 = t44 * t39 - t42 * t40;
	t46 = t42 * t39 + t44 * t40;
	t45 = cos(qJ(1));
	t43 = sin(qJ(1));
	t35 = t46 * t45;
	t34 = t47 * t45;
	t33 = t46 * t43;
	t32 = t47 * t43;
	t1 = [-t33, t34, 0, 0, -t34, 0; t35, t32, 0, 0, -t32, 0; 0, t46, 0, 0, -t46, 0; t32, t35, 0, 0, -t35, 0; -t34, t33, 0, 0, -t33, 0; 0, -t47, 0, 0, t47, 0; -t45, 0, 0, 0, 0, 0; -t43, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (106->21), mult. (134->24), div. (0->0), fcn. (202->8), ass. (0->26)
	t121 = pkin(10) + qJ(5);
	t120 = sin(t121);
	t125 = cos(qJ(2));
	t139 = sin(qJ(2));
	t140 = cos(t121);
	t115 = -t125 * t120 + t139 * t140;
	t114 = t139 * t120 + t125 * t140;
	t123 = sin(qJ(1));
	t110 = t115 * t123;
	t122 = sin(qJ(6));
	t138 = t110 * t122;
	t124 = cos(qJ(6));
	t137 = t110 * t124;
	t126 = cos(qJ(1));
	t113 = t115 * t126;
	t136 = t113 * t122;
	t135 = t113 * t124;
	t134 = t114 * t122;
	t133 = t114 * t124;
	t111 = t114 * t123;
	t128 = -t111 * t124 - t126 * t122;
	t127 = t111 * t122 - t126 * t124;
	t112 = t114 * t126;
	t109 = t112 * t124 - t123 * t122;
	t108 = -t112 * t122 - t123 * t124;
	t1 = [t128, -t135, 0, 0, t135, t108; t109, -t137, 0, 0, t137, -t127; 0, t133, 0, 0, -t133, -t115 * t122; t127, t136, 0, 0, -t136, -t109; t108, t138, 0, 0, -t138, t128; 0, -t134, 0, 0, t134, -t115 * t124; t110, -t112, 0, 0, t112, 0; -t113, -t111, 0, 0, t111, 0; 0, -t115, 0, 0, t115, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end