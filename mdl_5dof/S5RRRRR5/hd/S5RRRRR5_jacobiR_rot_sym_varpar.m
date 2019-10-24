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
% Datum: 2019-10-24 10:52
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:52:25
	% EndTime: 2019-10-24 10:52:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:25
	% EndTime: 2019-10-24 10:52:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t4, 0, 0, 0, 0; -t3, 0, 0, 0, 0; 0, 0, 0, 0, 0; t3, 0, 0, 0, 0; -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:26
	% EndTime: 2019-10-24 10:52:26
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
	% StartTime: 2019-10-24 10:52:25
	% EndTime: 2019-10-24 10:52:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->10), mult. (0->0), div. (0->0), fcn. (12->2), ass. (0->4)
	t14 = qJ(1) + qJ(2) + qJ(3);
	t13 = cos(t14);
	t12 = sin(t14);
	t1 = [0, 0, 0, 0, 0; -t13, -t13, -t13, 0, 0; -t12, -t12, -t12, 0, 0; 0, 0, 0, 0, 0; t12, t12, t12, 0, 0; -t13, -t13, -t13, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:26
	% EndTime: 2019-10-24 10:52:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (56->13), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t41 = qJ(1) + qJ(2) + qJ(3);
	t39 = sin(t41);
	t43 = cos(qJ(4));
	t45 = t39 * t43;
	t40 = cos(t41);
	t44 = t40 * t43;
	t42 = sin(qJ(4));
	t38 = t40 * t42;
	t37 = t39 * t42;
	t1 = [0, 0, 0, t43, 0; -t44, -t44, -t44, t37, 0; -t45, -t45, -t45, -t38, 0; 0, 0, 0, -t42, 0; t38, t38, t38, t45, 0; t37, t37, t37, -t44, 0; 0, 0, 0, 0, 0; -t39, -t39, -t39, 0, 0; t40, t40, t40, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:52:26
	% EndTime: 2019-10-24 10:52:26
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (91->17), mult. (20->4), div. (0->0), fcn. (50->4), ass. (0->11)
	t58 = qJ(1) + qJ(2) + qJ(3);
	t55 = cos(t58);
	t59 = qJ(4) + qJ(5);
	t57 = cos(t59);
	t60 = t55 * t57;
	t56 = sin(t59);
	t54 = sin(t58);
	t53 = t55 * t56;
	t52 = t54 * t57;
	t51 = t54 * t56;
	t1 = [0, 0, 0, t57, t57; -t60, -t60, -t60, t51, t51; -t52, -t52, -t52, -t53, -t53; 0, 0, 0, -t56, -t56; t53, t53, t53, t52, t52; t51, t51, t51, -t60, -t60; 0, 0, 0, 0, 0; -t54, -t54, -t54, 0, 0; t55, t55, t55, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end