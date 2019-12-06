% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRP5
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRPRP5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:39:18
	% EndTime: 2019-12-05 15:39:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:39:19
	% EndTime: 2019-12-05 15:39:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:39:18
	% EndTime: 2019-12-05 15:39:18
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
	% StartTime: 2019-12-05 15:39:19
	% EndTime: 2019-12-05 15:39:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (12->10), div. (0->0), fcn. (21->6), ass. (0->9)
	t35 = sin(pkin(7));
	t38 = sin(qJ(2));
	t41 = t35 * t38;
	t37 = cos(pkin(7));
	t40 = t37 * t38;
	t39 = cos(qJ(2));
	t36 = cos(pkin(8));
	t34 = sin(pkin(8));
	t1 = [0, -t36 * t40, 0, 0, 0; 0, -t36 * t41, 0, 0, 0; 0, t39 * t36, 0, 0, 0; 0, t34 * t40, 0, 0, 0; 0, t34 * t41, 0, 0, 0; 0, -t39 * t34, 0, 0, 0; 0, t37 * t39, 0, 0, 0; 0, t35 * t39, 0, 0, 0; 0, t38, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:39:19
	% EndTime: 2019-12-05 15:39:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (25->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t46 = pkin(8) + qJ(4);
	t44 = sin(t46);
	t49 = sin(qJ(2));
	t54 = t49 * t44;
	t45 = cos(t46);
	t53 = t49 * t45;
	t50 = cos(qJ(2));
	t52 = t50 * t44;
	t51 = t50 * t45;
	t48 = cos(pkin(7));
	t47 = sin(pkin(7));
	t1 = [0, -t48 * t53, 0, t47 * t45 - t48 * t52, 0; 0, -t47 * t53, 0, -t48 * t45 - t47 * t52, 0; 0, t51, 0, -t54, 0; 0, t48 * t54, 0, -t47 * t44 - t48 * t51, 0; 0, t47 * t54, 0, t48 * t44 - t47 * t51, 0; 0, -t52, 0, -t53, 0; 0, t48 * t50, 0, 0, 0; 0, t47 * t50, 0, 0, 0; 0, t49, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:39:19
	% EndTime: 2019-12-05 15:39:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (25->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t66 = pkin(8) + qJ(4);
	t64 = sin(t66);
	t69 = sin(qJ(2));
	t74 = t69 * t64;
	t65 = cos(t66);
	t73 = t69 * t65;
	t70 = cos(qJ(2));
	t72 = t70 * t64;
	t71 = t70 * t65;
	t68 = cos(pkin(7));
	t67 = sin(pkin(7));
	t1 = [0, -t68 * t73, 0, t67 * t65 - t68 * t72, 0; 0, -t67 * t73, 0, -t68 * t65 - t67 * t72, 0; 0, t71, 0, -t74, 0; 0, t68 * t70, 0, 0, 0; 0, t67 * t70, 0, 0, 0; 0, t69, 0, 0, 0; 0, -t68 * t74, 0, t67 * t64 + t68 * t71, 0; 0, -t67 * t74, 0, -t68 * t64 + t67 * t71, 0; 0, t72, 0, t73, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end