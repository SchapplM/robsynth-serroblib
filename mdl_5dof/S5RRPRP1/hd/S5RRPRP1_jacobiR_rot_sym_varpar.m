% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPRP1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:59:31
	% EndTime: 2020-01-03 11:59:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:59:31
	% EndTime: 2020-01-03 11:59:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t4, 0, 0, 0, 0; t3, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t3, 0, 0, 0, 0; t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:59:31
	% EndTime: 2020-01-03 11:59:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->3), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t11 = qJ(1) + qJ(2);
	t10 = cos(t11);
	t9 = sin(t11);
	t1 = [0, 0, 0, 0, 0; t10, t10, 0, 0, 0; t9, t9, 0, 0, 0; 0, 0, 0, 0, 0; -t9, -t9, 0, 0, 0; t10, t10, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:59:31
	% EndTime: 2020-01-03 11:59:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->3), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t14 = qJ(1) + qJ(2) + pkin(8);
	t13 = cos(t14);
	t12 = sin(t14);
	t1 = [0, 0, 0, 0, 0; t13, t13, 0, 0, 0; t12, t12, 0, 0, 0; 0, 0, 0, 0, 0; -t12, -t12, 0, 0, 0; t13, t13, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:59:31
	% EndTime: 2020-01-03 11:59:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (41->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t40 = qJ(1) + qJ(2) + pkin(8);
	t38 = sin(t40);
	t41 = sin(qJ(4));
	t44 = t38 * t41;
	t39 = cos(t40);
	t43 = t39 * t41;
	t42 = cos(qJ(4));
	t37 = t39 * t42;
	t36 = t38 * t42;
	t1 = [0, 0, 0, t42, 0; t37, t37, 0, -t44, 0; t36, t36, 0, t43, 0; 0, 0, 0, -t41, 0; -t43, -t43, 0, -t36, 0; -t44, -t44, 0, t37, 0; 0, 0, 0, 0, 0; t38, t38, 0, 0, 0; -t39, -t39, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:59:31
	% EndTime: 2020-01-03 11:59:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (41->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t47 = qJ(1) + qJ(2) + pkin(8);
	t45 = sin(t47);
	t48 = sin(qJ(4));
	t51 = t45 * t48;
	t46 = cos(t47);
	t50 = t46 * t48;
	t49 = cos(qJ(4));
	t44 = t46 * t49;
	t43 = t45 * t49;
	t1 = [0, 0, 0, t49, 0; t44, t44, 0, -t51, 0; t43, t43, 0, t50, 0; 0, 0, 0, -t48, 0; -t50, -t50, 0, -t43, 0; -t51, -t51, 0, t44, 0; 0, 0, 0, 0, 0; t45, t45, 0, 0, 0; -t46, -t46, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end