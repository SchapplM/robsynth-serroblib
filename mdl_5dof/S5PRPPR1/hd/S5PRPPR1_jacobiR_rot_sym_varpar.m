% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPPR1
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

function JR_rot = S5PRPPR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:08
	% EndTime: 2019-10-24 10:22:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:08
	% EndTime: 2019-10-24 10:22:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:08
	% EndTime: 2019-10-24 10:22:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = pkin(7) + qJ(2);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [0, -t10, 0, 0, 0; 0, t11, 0, 0, 0; 0, 0, 0, 0, 0; 0, -t11, 0, 0, 0; 0, -t10, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:08
	% EndTime: 2019-10-24 10:22:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->3), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t12 = cos(pkin(8));
	t11 = sin(pkin(8));
	t10 = pkin(7) + qJ(2);
	t9 = cos(t10);
	t8 = sin(t10);
	t1 = [0, -t8 * t12, 0, 0, 0; 0, t9 * t12, 0, 0, 0; 0, 0, 0, 0, 0; 0, t8 * t11, 0, 0, 0; 0, -t9 * t11, 0, 0, 0; 0, 0, 0, 0, 0; 0, t9, 0, 0, 0; 0, t8, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:08
	% EndTime: 2019-10-24 10:22:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (15->6), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->10)
	t40 = sin(pkin(9));
	t43 = cos(pkin(8));
	t45 = t40 * t43;
	t42 = cos(pkin(9));
	t44 = t42 * t43;
	t41 = sin(pkin(8));
	t39 = pkin(7) + qJ(2);
	t38 = cos(t39);
	t37 = sin(t39);
	t1 = [0, -t37 * t44 + t38 * t40, 0, 0, 0; 0, t37 * t40 + t38 * t44, 0, 0, 0; 0, 0, 0, 0, 0; 0, t37 * t45 + t38 * t42, 0, 0, 0; 0, t37 * t42 - t38 * t45, 0, 0, 0; 0, 0, 0, 0, 0; 0, -t37 * t41, 0, 0, 0; 0, t38 * t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:22:08
	% EndTime: 2019-10-24 10:22:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (47->11), mult. (28->14), div. (0->0), fcn. (48->6), ass. (0->15)
	t60 = pkin(7) + qJ(2);
	t56 = sin(t60);
	t62 = cos(pkin(8));
	t64 = t56 * t62;
	t58 = cos(t60);
	t63 = t58 * t62;
	t61 = sin(pkin(8));
	t59 = pkin(9) + qJ(5);
	t57 = cos(t59);
	t55 = sin(t59);
	t54 = t56 * t55 + t57 * t63;
	t53 = -t55 * t63 + t56 * t57;
	t52 = t58 * t55 - t57 * t64;
	t51 = t55 * t64 + t58 * t57;
	t1 = [0, t52, 0, 0, t53; 0, t54, 0, 0, -t51; 0, 0, 0, 0, -t61 * t55; 0, t51, 0, 0, -t54; 0, t53, 0, 0, t52; 0, 0, 0, 0, -t61 * t57; 0, -t56 * t61, 0, 0, 0; 0, t58 * t61, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end