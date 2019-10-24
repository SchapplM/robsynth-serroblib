% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRP4
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
% Datum: 2019-10-24 10:24
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRPRP4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:26
	% EndTime: 2019-10-24 10:24:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:25
	% EndTime: 2019-10-24 10:24:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:24:26
	% EndTime: 2019-10-24 10:24:26
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
	% StartTime: 2019-10-24 10:24:25
	% EndTime: 2019-10-24 10:24:25
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
	% StartTime: 2019-10-24 10:24:26
	% EndTime: 2019-10-24 10:24:26
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
	% StartTime: 2019-10-24 10:24:26
	% EndTime: 2019-10-24 10:24:26
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t71 = sin(pkin(7));
	t73 = sin(qJ(4));
	t78 = t71 * t73;
	t74 = cos(qJ(4));
	t77 = t71 * t74;
	t72 = cos(pkin(7));
	t76 = t72 * t73;
	t75 = t72 * t74;
	t70 = qJ(2) + pkin(8);
	t69 = cos(t70);
	t68 = sin(t70);
	t1 = [0, -t68 * t75, 0, -t69 * t76 + t77, 0; 0, -t68 * t77, 0, -t69 * t78 - t75, 0; 0, t69 * t74, 0, -t68 * t73, 0; 0, t72 * t69, 0, 0, 0; 0, t71 * t69, 0, 0, 0; 0, t68, 0, 0, 0; 0, -t68 * t76, 0, t69 * t75 + t78, 0; 0, -t68 * t78, 0, t69 * t77 - t76, 0; 0, t69 * t73, 0, t68 * t74, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end