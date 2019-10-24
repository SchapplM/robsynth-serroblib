% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:42
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRPR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
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
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t6 = qJ(1) + pkin(8);
	t5 = cos(t6);
	t4 = sin(t6);
	t1 = [0, 0, 0, 0, 0; -t5, 0, 0, 0, 0; -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; t4, 0, 0, 0, 0; -t5, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t10 = qJ(1) + pkin(8) + qJ(3);
	t9 = cos(t10);
	t8 = sin(t10);
	t1 = [0, 0, 0, 0, 0; -t9, 0, -t9, 0, 0; -t8, 0, -t8, 0, 0; 0, 0, 0, 0, 0; t8, 0, t8, 0, 0; -t9, 0, -t9, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t28 = qJ(1) + pkin(8) + qJ(3);
	t26 = sin(t28);
	t30 = cos(pkin(9));
	t32 = t26 * t30;
	t27 = cos(t28);
	t31 = t27 * t30;
	t29 = sin(pkin(9));
	t25 = t27 * t29;
	t24 = t26 * t29;
	t1 = [0, 0, 0, 0, 0; -t31, 0, -t31, 0, 0; -t32, 0, -t32, 0, 0; 0, 0, 0, 0, 0; t25, 0, t25, 0, 0; t24, 0, t24, 0, 0; 0, 0, 0, 0, 0; -t26, 0, -t26, 0, 0; t27, 0, t27, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:42:13
	% EndTime: 2019-10-24 10:42:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (74->14), mult. (42->14), div. (0->0), fcn. (72->6), ass. (0->16)
	t52 = qJ(1) + pkin(8) + qJ(3);
	t50 = sin(t52);
	t53 = sin(pkin(9));
	t60 = t50 * t53;
	t51 = cos(t52);
	t59 = t51 * t53;
	t54 = cos(pkin(9));
	t55 = sin(qJ(5));
	t58 = t54 * t55;
	t56 = cos(qJ(5));
	t57 = t54 * t56;
	t48 = -t50 * t55 - t51 * t57;
	t47 = -t50 * t56 + t51 * t58;
	t46 = t50 * t57 - t51 * t55;
	t45 = t50 * t58 + t51 * t56;
	t1 = [0, 0, 0, 0, -t53 * t55; t48, 0, t48, 0, t45; -t46, 0, -t46, 0, -t47; 0, 0, 0, 0, -t53 * t56; t47, 0, t47, 0, t46; t45, 0, t45, 0, t48; 0, 0, 0, 0, 0; -t59, 0, -t59, 0, 0; -t60, 0, -t60, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end