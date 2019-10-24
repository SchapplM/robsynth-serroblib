% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPPR2
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:22
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRPPR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:30
	% EndTime: 2019-10-24 10:22:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:30
	% EndTime: 2019-10-24 10:22:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:31
	% EndTime: 2019-10-24 10:22:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t9 = cos(qJ(2));
	t8 = sin(qJ(2));
	t7 = cos(pkin(7));
	t6 = sin(pkin(7));
	t1 = [0, -t7 * t8, 0, 0, 0; 0, -t6 * t8, 0, 0, 0; 0, t9, 0, 0, 0; 0, -t7 * t9, 0, 0, 0; 0, -t6 * t9, 0, 0, 0; 0, -t8, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:30
	% EndTime: 2019-10-24 10:22:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->6), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t12 = cos(pkin(7));
	t11 = sin(pkin(7));
	t10 = qJ(2) + pkin(8);
	t9 = cos(t10);
	t8 = sin(t10);
	t1 = [0, -t12 * t8, 0, 0, 0; 0, -t11 * t8, 0, 0, 0; 0, t9, 0, 0, 0; 0, -t12 * t9, 0, 0, 0; 0, -t11 * t9, 0, 0, 0; 0, -t8, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:31
	% EndTime: 2019-10-24 10:22:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->4), mult. (12->10), div. (0->0), fcn. (21->6), ass. (0->10)
	t39 = qJ(2) + pkin(8);
	t37 = sin(t39);
	t41 = sin(pkin(7));
	t45 = t37 * t41;
	t43 = cos(pkin(7));
	t44 = t37 * t43;
	t42 = cos(pkin(9));
	t40 = sin(pkin(9));
	t38 = cos(t39);
	t1 = [0, -t42 * t44, 0, 0, 0; 0, -t42 * t45, 0, 0, 0; 0, t38 * t42, 0, 0, 0; 0, t40 * t44, 0, 0, 0; 0, t40 * t45, 0, 0, 0; 0, -t38 * t40, 0, 0, 0; 0, t43 * t38, 0, 0, 0; 0, t41 * t38, 0, 0, 0; 0, t37, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:31
	% EndTime: 2019-10-24 10:22:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (40->11), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->15)
	t52 = pkin(9) + qJ(5);
	t48 = sin(t52);
	t54 = sin(pkin(7));
	t61 = t54 * t48;
	t50 = cos(t52);
	t60 = t54 * t50;
	t53 = qJ(2) + pkin(8);
	t51 = cos(t53);
	t59 = t54 * t51;
	t55 = cos(pkin(7));
	t58 = t55 * t48;
	t57 = t55 * t50;
	t56 = t55 * t51;
	t49 = sin(t53);
	t1 = [0, -t49 * t57, 0, 0, -t48 * t56 + t60; 0, -t49 * t60, 0, 0, -t48 * t59 - t57; 0, t51 * t50, 0, 0, -t49 * t48; 0, t49 * t58, 0, 0, -t50 * t56 - t61; 0, t49 * t61, 0, 0, -t50 * t59 + t58; 0, -t51 * t48, 0, 0, -t49 * t50; 0, t56, 0, 0, 0; 0, t59, 0, 0, 0; 0, t49, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end