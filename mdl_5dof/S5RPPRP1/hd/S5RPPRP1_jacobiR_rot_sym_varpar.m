% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPRP1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
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
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t6 = qJ(1) + pkin(7);
	t5 = cos(t6);
	t4 = sin(t6);
	t1 = [0, 0, 0, 0, 0; t5, 0, 0, 0, 0; t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t4, 0, 0, 0, 0; t5, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->4), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t18 = cos(pkin(8));
	t17 = sin(pkin(8));
	t16 = qJ(1) + pkin(7);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [0, 0, 0, 0, 0; t15 * t18, 0, 0, 0, 0; t14 * t18, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t15 * t17, 0, 0, 0, 0; -t14 * t17, 0, 0, 0, 0; 0, 0, 0, 0, 0; t14, 0, 0, 0, 0; -t15, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:16
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (28->9), mult. (28->14), div. (0->0), fcn. (48->6), ass. (0->14)
	t46 = cos(pkin(8));
	t47 = sin(qJ(4));
	t50 = t46 * t47;
	t48 = cos(qJ(4));
	t49 = t46 * t48;
	t45 = sin(pkin(8));
	t44 = qJ(1) + pkin(7);
	t43 = cos(t44);
	t42 = sin(t44);
	t41 = t42 * t47 + t43 * t49;
	t40 = -t42 * t48 + t43 * t50;
	t39 = t42 * t49 - t43 * t47;
	t38 = -t42 * t50 - t43 * t48;
	t1 = [0, 0, 0, -t45 * t47, 0; t41, 0, 0, t38, 0; t39, 0, 0, t40, 0; 0, 0, 0, -t45 * t48, 0; -t40, 0, 0, -t39, 0; t38, 0, 0, t41, 0; 0, 0, 0, 0, 0; t43 * t45, 0, 0, 0, 0; t42 * t45, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:26:17
	% EndTime: 2020-01-03 11:26:17
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (28->9), mult. (28->14), div. (0->0), fcn. (48->6), ass. (0->14)
	t48 = cos(pkin(8));
	t49 = sin(qJ(4));
	t52 = t48 * t49;
	t50 = cos(qJ(4));
	t51 = t48 * t50;
	t47 = sin(pkin(8));
	t46 = qJ(1) + pkin(7);
	t45 = cos(t46);
	t44 = sin(t46);
	t43 = t44 * t49 + t45 * t51;
	t42 = -t44 * t50 + t45 * t52;
	t41 = t44 * t51 - t45 * t49;
	t40 = -t44 * t52 - t45 * t50;
	t1 = [0, 0, 0, -t47 * t49, 0; t43, 0, 0, t40, 0; t41, 0, 0, t42, 0; 0, 0, 0, -t47 * t50, 0; -t42, 0, 0, -t41, 0; t40, 0, 0, t43, 0; 0, 0, 0, 0, 0; t45 * t47, 0, 0, 0, 0; t44 * t47, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end