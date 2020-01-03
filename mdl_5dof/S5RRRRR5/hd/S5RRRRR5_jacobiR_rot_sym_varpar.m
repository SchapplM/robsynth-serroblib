% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
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
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
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
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (27->4), mult. (0->0), div. (0->0), fcn. (12->2), ass. (0->4)
	t16 = qJ(1) + qJ(2) + qJ(3);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [0, 0, 0, 0, 0; t15, t15, t15, 0, 0; t14, t14, t14, 0, 0; 0, 0, 0, 0, 0; -t14, -t14, -t14, 0, 0; t15, t15, t15, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (56->13), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t42 = qJ(1) + qJ(2) + qJ(3);
	t40 = sin(t42);
	t43 = sin(qJ(4));
	t46 = t40 * t43;
	t41 = cos(t42);
	t45 = t41 * t43;
	t44 = cos(qJ(4));
	t39 = t41 * t44;
	t38 = t40 * t44;
	t1 = [0, 0, 0, t44, 0; t39, t39, t39, -t46, 0; t38, t38, t38, t45, 0; 0, 0, 0, -t43, 0; -t45, -t45, -t45, -t38, 0; -t46, -t46, -t46, t39, 0; 0, 0, 0, 0, 0; t40, t40, t40, 0, 0; -t41, -t41, -t41, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (91->17), mult. (20->4), div. (0->0), fcn. (50->4), ass. (0->11)
	t59 = qJ(1) + qJ(2) + qJ(3);
	t55 = sin(t59);
	t60 = qJ(4) + qJ(5);
	t57 = sin(t60);
	t61 = t55 * t57;
	t58 = cos(t60);
	t56 = cos(t59);
	t54 = t56 * t58;
	t53 = t56 * t57;
	t52 = t55 * t58;
	t1 = [0, 0, 0, t58, t58; t54, t54, t54, -t61, -t61; t52, t52, t52, t53, t53; 0, 0, 0, -t57, -t57; -t53, -t53, -t53, -t52, -t52; -t61, -t61, -t61, t54, t54; 0, 0, 0, 0, 0; t55, t55, t55, 0, 0; -t56, -t56, -t56, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end