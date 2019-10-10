% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:15
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
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
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
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
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (18->12), mult. (56->8), div. (0->0), fcn. (90->6), ass. (0->13)
	t35 = sin(qJ(4));
	t36 = sin(qJ(2));
	t38 = cos(qJ(4));
	t39 = cos(qJ(2));
	t42 = t39 * t35 - t36 * t38;
	t41 = t36 * t35 + t39 * t38;
	t40 = cos(qJ(1));
	t37 = sin(qJ(1));
	t31 = t41 * t40;
	t30 = t42 * t40;
	t29 = t41 * t37;
	t28 = t42 * t37;
	t1 = [-t29, t30, 0, -t30, 0, 0; t31, t28, 0, -t28, 0, 0; 0, t41, 0, -t41, 0, 0; t28, t31, 0, -t31, 0, 0; -t30, t29, 0, -t29, 0, 0; 0, -t42, 0, t42, 0, 0; -t40, 0, 0, 0, 0, 0; -t37, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (50->13), mult. (56->8), div. (0->0), fcn. (90->6), ass. (0->14)
	t43 = qJ(4) + pkin(10);
	t41 = sin(t43);
	t42 = cos(t43);
	t44 = sin(qJ(2));
	t46 = cos(qJ(2));
	t49 = t46 * t41 - t44 * t42;
	t48 = t44 * t41 + t46 * t42;
	t47 = cos(qJ(1));
	t45 = sin(qJ(1));
	t37 = t48 * t47;
	t36 = t49 * t47;
	t35 = t48 * t45;
	t34 = t49 * t45;
	t1 = [-t35, t36, 0, -t36, 0, 0; t37, t34, 0, -t34, 0, 0; 0, t48, 0, -t48, 0, 0; t34, t37, 0, -t37, 0, 0; -t36, t35, 0, -t35, 0, 0; 0, -t49, 0, t49, 0, 0; -t47, 0, 0, 0, 0, 0; -t45, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:24
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (106->21), mult. (134->24), div. (0->0), fcn. (202->8), ass. (0->26)
	t123 = qJ(4) + pkin(10);
	t122 = sin(t123);
	t127 = cos(qJ(2));
	t141 = sin(qJ(2));
	t142 = cos(t123);
	t117 = -t127 * t122 + t141 * t142;
	t116 = t141 * t122 + t127 * t142;
	t125 = sin(qJ(1));
	t112 = t117 * t125;
	t124 = sin(qJ(6));
	t140 = t112 * t124;
	t126 = cos(qJ(6));
	t139 = t112 * t126;
	t128 = cos(qJ(1));
	t115 = t117 * t128;
	t138 = t115 * t124;
	t137 = t115 * t126;
	t136 = t116 * t124;
	t135 = t116 * t126;
	t113 = t116 * t125;
	t130 = -t113 * t126 - t128 * t124;
	t129 = t113 * t124 - t128 * t126;
	t114 = t116 * t128;
	t111 = t114 * t126 - t125 * t124;
	t110 = -t114 * t124 - t125 * t126;
	t1 = [t130, -t137, 0, t137, 0, t110; t111, -t139, 0, t139, 0, -t129; 0, t135, 0, -t135, 0, -t117 * t124; t129, t138, 0, -t138, 0, -t111; t110, t140, 0, -t140, 0, t130; 0, -t136, 0, t136, 0, -t117 * t126; t112, -t114, 0, t114, 0, 0; -t115, -t113, 0, t113, 0, 0; 0, -t117, 0, t117, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end