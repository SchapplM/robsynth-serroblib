% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR6
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
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:16:04
	% EndTime: 2020-01-03 12:16:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:16:04
	% EndTime: 2020-01-03 12:16:04
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
	% StartTime: 2020-01-03 12:16:04
	% EndTime: 2020-01-03 12:16:04
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
	% StartTime: 2020-01-03 12:16:04
	% EndTime: 2020-01-03 12:16:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (25->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t37 = qJ(1) + qJ(2);
	t35 = sin(t37);
	t38 = sin(qJ(3));
	t41 = t35 * t38;
	t36 = cos(t37);
	t40 = t36 * t38;
	t39 = cos(qJ(3));
	t34 = t36 * t39;
	t33 = t35 * t39;
	t1 = [0, 0, t39, 0, 0; t34, t34, -t41, 0, 0; t33, t33, t40, 0, 0; 0, 0, -t38, 0, 0; -t40, -t40, -t33, 0, 0; -t41, -t41, t34, 0, 0; 0, 0, 0, 0, 0; t35, t35, 0, 0, 0; -t36, -t36, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:16:04
	% EndTime: 2020-01-03 12:16:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (52->14), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->11)
	t54 = qJ(3) + qJ(4);
	t50 = sin(t54);
	t55 = qJ(1) + qJ(2);
	t51 = sin(t55);
	t56 = t51 * t50;
	t53 = cos(t55);
	t52 = cos(t54);
	t49 = t53 * t52;
	t48 = t53 * t50;
	t47 = t51 * t52;
	t1 = [0, 0, t52, t52, 0; t49, t49, -t56, -t56, 0; t47, t47, t48, t48, 0; 0, 0, -t50, -t50, 0; -t48, -t48, -t47, -t47, 0; -t56, -t56, t49, t49, 0; 0, 0, 0, 0, 0; t51, t51, 0, 0, 0; -t53, -t53, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:16:04
	% EndTime: 2020-01-03 12:16:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (91->17), mult. (20->4), div. (0->0), fcn. (50->4), ass. (0->11)
	t58 = qJ(3) + qJ(4) + qJ(5);
	t54 = sin(t58);
	t59 = qJ(1) + qJ(2);
	t56 = sin(t59);
	t60 = t56 * t54;
	t57 = cos(t59);
	t55 = cos(t58);
	t53 = t57 * t55;
	t52 = t57 * t54;
	t51 = t56 * t55;
	t1 = [0, 0, t55, t55, t55; t53, t53, -t60, -t60, -t60; t51, t51, t52, t52, t52; 0, 0, -t54, -t54, -t54; -t52, -t52, -t51, -t51, -t51; -t60, -t60, t53, t53, t53; 0, 0, 0, 0, 0; t56, t56, 0, 0, 0; -t57, -t57, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end