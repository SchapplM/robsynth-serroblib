% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRPP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRPP2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:55
	% EndTime: 2019-12-05 16:10:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:55
	% EndTime: 2019-12-05 16:10:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:55
	% EndTime: 2019-12-05 16:10:55
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
	% StartTime: 2019-12-05 16:10:56
	% EndTime: 2019-12-05 16:10:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->9), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->11)
	t42 = sin(qJ(3));
	t43 = sin(qJ(2));
	t49 = t43 * t42;
	t44 = cos(qJ(3));
	t48 = t43 * t44;
	t45 = cos(qJ(2));
	t47 = t45 * t42;
	t46 = t45 * t44;
	t41 = cos(pkin(7));
	t40 = sin(pkin(7));
	t1 = [0, -t41 * t48, t40 * t44 - t41 * t47, 0, 0; 0, -t40 * t48, -t40 * t47 - t41 * t44, 0, 0; 0, t46, -t49, 0, 0; 0, t41 * t49, -t40 * t42 - t41 * t46, 0, 0; 0, t40 * t49, -t40 * t46 + t41 * t42, 0, 0; 0, -t47, -t48, 0, 0; 0, t41 * t45, 0, 0, 0; 0, t40 * t45, 0, 0, 0; 0, t43, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:55
	% EndTime: 2019-12-05 16:10:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (25->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t49 = qJ(3) + pkin(8);
	t47 = sin(t49);
	t52 = sin(qJ(2));
	t57 = t52 * t47;
	t48 = cos(t49);
	t56 = t52 * t48;
	t53 = cos(qJ(2));
	t55 = t53 * t47;
	t54 = t53 * t48;
	t51 = cos(pkin(7));
	t50 = sin(pkin(7));
	t1 = [0, -t51 * t56, t50 * t48 - t51 * t55, 0, 0; 0, -t50 * t56, -t51 * t48 - t50 * t55, 0, 0; 0, t54, -t57, 0, 0; 0, t51 * t57, -t50 * t47 - t51 * t54, 0, 0; 0, t50 * t57, t51 * t47 - t50 * t54, 0, 0; 0, -t55, -t56, 0, 0; 0, t51 * t53, 0, 0, 0; 0, t50 * t53, 0, 0, 0; 0, t52, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:10:56
	% EndTime: 2019-12-05 16:10:56
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (25->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t69 = qJ(3) + pkin(8);
	t67 = sin(t69);
	t72 = sin(qJ(2));
	t77 = t72 * t67;
	t68 = cos(t69);
	t76 = t72 * t68;
	t73 = cos(qJ(2));
	t75 = t73 * t67;
	t74 = t73 * t68;
	t71 = cos(pkin(7));
	t70 = sin(pkin(7));
	t1 = [0, -t71 * t76, t70 * t68 - t71 * t75, 0, 0; 0, -t70 * t76, -t71 * t68 - t70 * t75, 0, 0; 0, t74, -t77, 0, 0; 0, t71 * t73, 0, 0, 0; 0, t70 * t73, 0, 0, 0; 0, t72, 0, 0, 0; 0, -t71 * t77, t70 * t67 + t71 * t74, 0, 0; 0, -t70 * t77, -t71 * t67 + t70 * t74, 0, 0; 0, t75, t76, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end