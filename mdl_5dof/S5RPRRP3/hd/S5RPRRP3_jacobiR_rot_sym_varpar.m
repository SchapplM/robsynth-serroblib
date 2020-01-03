% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRRP3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t4, 0, 0, 0, 0; t3, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t3, 0, 0, 0, 0; t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (5->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t6 = qJ(1) + pkin(8);
	t5 = cos(t6);
	t4 = sin(t6);
	t1 = [0, 0, 0, 0, 0; t5, 0, 0, 0, 0; t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t4, 0, 0, 0, 0; t5, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t25 = qJ(1) + pkin(8);
	t23 = sin(t25);
	t26 = sin(qJ(3));
	t30 = t23 * t26;
	t27 = cos(qJ(3));
	t29 = t23 * t27;
	t24 = cos(t25);
	t28 = t24 * t26;
	t22 = t24 * t27;
	t1 = [0, 0, t27, 0, 0; t22, 0, -t30, 0, 0; t29, 0, t28, 0, 0; 0, 0, -t26, 0, 0; -t28, 0, -t29, 0, 0; -t30, 0, t22, 0, 0; 0, 0, 0, 0, 0; t23, 0, 0, 0, 0; -t24, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (39->11), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t47 = qJ(1) + pkin(8);
	t43 = sin(t47);
	t48 = qJ(3) + qJ(4);
	t45 = sin(t48);
	t50 = t43 * t45;
	t46 = cos(t48);
	t49 = t43 * t46;
	t44 = cos(t47);
	t42 = t44 * t46;
	t41 = t44 * t45;
	t1 = [0, 0, t46, t46, 0; t42, 0, -t50, -t50, 0; t49, 0, t41, t41, 0; 0, 0, -t45, -t45, 0; -t41, 0, -t49, -t49, 0; -t50, 0, t42, t42, 0; 0, 0, 0, 0, 0; t43, 0, 0, 0, 0; -t44, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (39->11), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t46 = qJ(1) + pkin(8);
	t42 = sin(t46);
	t47 = qJ(3) + qJ(4);
	t44 = sin(t47);
	t49 = t42 * t44;
	t45 = cos(t47);
	t48 = t42 * t45;
	t43 = cos(t46);
	t41 = t43 * t45;
	t40 = t43 * t44;
	t1 = [0, 0, t45, t45, 0; t41, 0, -t49, -t49, 0; t48, 0, t40, t40, 0; 0, 0, -t44, -t44, 0; -t40, 0, -t48, -t48, 0; -t49, 0, t41, t41, 0; 0, 0, 0, 0, 0; t42, 0, 0, 0, 0; -t43, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end