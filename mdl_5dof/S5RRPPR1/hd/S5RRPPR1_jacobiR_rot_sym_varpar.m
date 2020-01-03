% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:56:27
	% EndTime: 2020-01-03 11:56:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:56:26
	% EndTime: 2020-01-03 11:56:26
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
	% StartTime: 2020-01-03 11:56:27
	% EndTime: 2020-01-03 11:56:27
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
	% StartTime: 2020-01-03 11:56:26
	% EndTime: 2020-01-03 11:56:27
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
	% StartTime: 2020-01-03 11:56:27
	% EndTime: 2020-01-03 11:56:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t32 = qJ(1) + qJ(2) + pkin(8);
	t30 = sin(t32);
	t33 = sin(pkin(9));
	t36 = t30 * t33;
	t31 = cos(t32);
	t35 = t31 * t33;
	t34 = cos(pkin(9));
	t29 = t31 * t34;
	t28 = t30 * t34;
	t1 = [0, 0, 0, 0, 0; t29, t29, 0, 0, 0; t28, t28, 0, 0, 0; 0, 0, 0, 0, 0; -t35, -t35, 0, 0, 0; -t36, -t36, 0, 0, 0; 0, 0, 0, 0, 0; t30, t30, 0, 0, 0; -t31, -t31, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:56:27
	% EndTime: 2020-01-03 11:56:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (55->11), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t46 = qJ(1) + qJ(2) + pkin(8);
	t42 = sin(t46);
	t47 = pkin(9) + qJ(5);
	t44 = sin(t47);
	t49 = t42 * t44;
	t43 = cos(t46);
	t48 = t43 * t44;
	t45 = cos(t47);
	t41 = t43 * t45;
	t40 = t42 * t45;
	t1 = [0, 0, 0, 0, t45; t41, t41, 0, 0, -t49; t40, t40, 0, 0, t48; 0, 0, 0, 0, -t44; -t48, -t48, 0, 0, -t40; -t49, -t49, 0, 0, t41; 0, 0, 0, 0, 0; t42, t42, 0, 0, 0; -t43, -t43, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end