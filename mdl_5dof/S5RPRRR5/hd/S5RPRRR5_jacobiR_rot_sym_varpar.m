% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRRR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
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
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t6 = qJ(1) + pkin(9);
	t5 = cos(t6);
	t4 = sin(t6);
	t1 = [0, 0, 0, 0, 0; t5, 0, 0, 0, 0; t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t4, 0, 0, 0, 0; t5, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->3), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t12 = qJ(1) + pkin(9) + qJ(3);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [0, 0, 0, 0, 0; t11, 0, t11, 0, 0; t10, 0, t10, 0, 0; 0, 0, 0, 0, 0; -t10, 0, -t10, 0, 0; t11, 0, t11, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (41->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t38 = qJ(1) + pkin(9) + qJ(3);
	t36 = sin(t38);
	t39 = sin(qJ(4));
	t42 = t36 * t39;
	t37 = cos(t38);
	t41 = t37 * t39;
	t40 = cos(qJ(4));
	t35 = t37 * t40;
	t34 = t36 * t40;
	t1 = [0, 0, 0, t40, 0; t35, 0, t35, -t42, 0; t34, 0, t34, t41, 0; 0, 0, 0, -t39, 0; -t41, 0, -t41, -t34, 0; -t42, 0, -t42, t35, 0; 0, 0, 0, 0, 0; t36, 0, t36, 0, 0; -t37, 0, -t37, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (72->14), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->11)
	t53 = qJ(1) + pkin(9) + qJ(3);
	t51 = sin(t53);
	t56 = qJ(4) + qJ(5);
	t54 = sin(t56);
	t57 = t51 * t54;
	t55 = cos(t56);
	t52 = cos(t53);
	t50 = t52 * t55;
	t49 = t52 * t54;
	t48 = t51 * t55;
	t1 = [0, 0, 0, t55, t55; t50, 0, t50, -t57, -t57; t48, 0, t48, t49, t49; 0, 0, 0, -t54, -t54; -t49, 0, -t49, -t48, -t48; -t57, 0, -t57, t50, t50; 0, 0, 0, 0, 0; t51, 0, t51, 0, 0; -t52, 0, -t52, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end