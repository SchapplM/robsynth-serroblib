% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRP3
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRPRP3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:34:00
	% EndTime: 2019-12-05 15:34:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:34:01
	% EndTime: 2019-12-05 15:34:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:34:01
	% EndTime: 2019-12-05 15:34:01
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
	% StartTime: 2019-12-05 15:34:01
	% EndTime: 2019-12-05 15:34:01
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
	% StartTime: 2019-12-05 15:34:01
	% EndTime: 2019-12-05 15:34:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t49 = sin(pkin(7));
	t51 = sin(qJ(4));
	t56 = t49 * t51;
	t52 = cos(qJ(4));
	t55 = t49 * t52;
	t50 = cos(pkin(7));
	t54 = t50 * t51;
	t53 = t50 * t52;
	t48 = qJ(2) + pkin(8);
	t47 = cos(t48);
	t46 = sin(t48);
	t1 = [0, -t46 * t53, 0, -t47 * t54 + t55, 0; 0, -t46 * t55, 0, -t47 * t56 - t53, 0; 0, t47 * t52, 0, -t46 * t51, 0; 0, t46 * t54, 0, -t47 * t53 - t56, 0; 0, t46 * t56, 0, -t47 * t55 + t54, 0; 0, -t47 * t51, 0, -t46 * t52, 0; 0, t50 * t47, 0, 0, 0; 0, t49 * t47, 0, 0, 0; 0, t46, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:34:01
	% EndTime: 2019-12-05 15:34:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t54 = sin(pkin(7));
	t56 = sin(qJ(4));
	t61 = t54 * t56;
	t57 = cos(qJ(4));
	t60 = t54 * t57;
	t55 = cos(pkin(7));
	t59 = t55 * t56;
	t58 = t55 * t57;
	t53 = qJ(2) + pkin(8);
	t52 = cos(t53);
	t51 = sin(t53);
	t1 = [0, -t51 * t58, 0, -t52 * t59 + t60, 0; 0, -t51 * t60, 0, -t52 * t61 - t58, 0; 0, t52 * t57, 0, -t51 * t56, 0; 0, t51 * t59, 0, -t52 * t58 - t61, 0; 0, t51 * t61, 0, -t52 * t60 + t59, 0; 0, -t52 * t56, 0, -t51 * t57, 0; 0, t55 * t52, 0, 0, 0; 0, t54 * t52, 0, 0, 0; 0, t51, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end