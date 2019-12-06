% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRPR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:43:51
	% EndTime: 2019-12-05 18:43:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:43:52
	% EndTime: 2019-12-05 18:43:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t4, 0, 0, 0, 0; -t3, 0, 0, 0, 0; 0, 0, 0, 0, 0; t3, 0, 0, 0, 0; -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:43:51
	% EndTime: 2019-12-05 18:43:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (14->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t9 = qJ(1) + qJ(2);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; -t8, -t8, 0, 0, 0; -t7, -t7, 0, 0, 0; 0, 0, 0, 0, 0; t7, t7, 0, 0, 0; -t8, -t8, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:43:52
	% EndTime: 2019-12-05 18:43:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (25->10), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
	t36 = qJ(1) + qJ(2);
	t34 = sin(t36);
	t38 = cos(qJ(3));
	t40 = t34 * t38;
	t35 = cos(t36);
	t39 = t35 * t38;
	t37 = sin(qJ(3));
	t33 = t35 * t37;
	t32 = t34 * t37;
	t1 = [0, 0, t38, 0, 0; -t39, -t39, t32, 0, 0; -t40, -t40, -t33, 0, 0; 0, 0, -t37, 0, 0; t33, t33, t40, 0, 0; t32, t32, -t39, 0, 0; 0, 0, 0, 0, 0; -t34, -t34, 0, 0, 0; t35, t35, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:43:51
	% EndTime: 2019-12-05 18:43:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (39->11), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
	t45 = qJ(3) + pkin(9);
	t42 = cos(t45);
	t46 = qJ(1) + qJ(2);
	t43 = sin(t46);
	t48 = t43 * t42;
	t44 = cos(t46);
	t47 = t44 * t42;
	t41 = sin(t45);
	t40 = t44 * t41;
	t39 = t43 * t41;
	t1 = [0, 0, t42, 0, 0; -t47, -t47, t39, 0, 0; -t48, -t48, -t40, 0, 0; 0, 0, -t41, 0, 0; t40, t40, t48, 0, 0; t39, t39, -t47, 0, 0; 0, 0, 0, 0, 0; -t43, -t43, 0, 0, 0; t44, t44, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:43:52
	% EndTime: 2019-12-05 18:43:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (72->14), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->11)
	t53 = qJ(3) + pkin(9) + qJ(5);
	t52 = cos(t53);
	t56 = qJ(1) + qJ(2);
	t55 = cos(t56);
	t57 = t55 * t52;
	t54 = sin(t56);
	t51 = sin(t53);
	t50 = t55 * t51;
	t49 = t54 * t52;
	t48 = t54 * t51;
	t1 = [0, 0, t52, 0, t52; -t57, -t57, t48, 0, t48; -t49, -t49, -t50, 0, -t50; 0, 0, -t51, 0, -t51; t50, t50, t49, 0, t49; t48, t48, -t57, 0, -t57; 0, 0, 0, 0, 0; -t54, -t54, 0, 0, 0; t55, t55, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end