% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR7
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:37
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRRR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:42
	% EndTime: 2019-10-24 10:37:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:42
	% EndTime: 2019-10-24 10:37:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t9 = cos(qJ(2));
	t8 = sin(qJ(2));
	t7 = cos(pkin(9));
	t6 = sin(pkin(9));
	t1 = [0, -t7 * t8, 0, 0, 0; 0, -t6 * t8, 0, 0, 0; 0, t9, 0, 0, 0; 0, -t7 * t9, 0, 0, 0; 0, -t6 * t9, 0, 0, 0; 0, -t8, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
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
	t41 = cos(pkin(9));
	t40 = sin(pkin(9));
	t1 = [0, -t41 * t48, t40 * t44 - t41 * t47, 0, 0; 0, -t40 * t48, -t40 * t47 - t41 * t44, 0, 0; 0, t46, -t49, 0, 0; 0, t41 * t49, -t40 * t42 - t41 * t46, 0, 0; 0, t40 * t49, -t40 * t46 + t41 * t42, 0, 0; 0, -t47, -t48, 0, 0; 0, t41 * t45, 0, 0, 0; 0, t40 * t45, 0, 0, 0; 0, t43, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (41->12), mult. (40->18), div. (0->0), fcn. (69->6), ass. (0->16)
	t64 = qJ(3) + qJ(4);
	t62 = sin(t64);
	t67 = sin(qJ(2));
	t72 = t67 * t62;
	t63 = cos(t64);
	t71 = t67 * t63;
	t68 = cos(qJ(2));
	t70 = t68 * t62;
	t69 = t68 * t63;
	t66 = cos(pkin(9));
	t65 = sin(pkin(9));
	t61 = -t65 * t62 - t66 * t69;
	t60 = t65 * t63 - t66 * t70;
	t59 = t66 * t62 - t65 * t69;
	t58 = -t66 * t63 - t65 * t70;
	t1 = [0, -t66 * t71, t60, t60, 0; 0, -t65 * t71, t58, t58, 0; 0, t69, -t72, -t72, 0; 0, t66 * t72, t61, t61, 0; 0, t65 * t72, t59, t59, 0; 0, -t70, -t71, -t71, 0; 0, t66 * t68, 0, 0, 0; 0, t65 * t68, 0, 0, 0; 0, t67, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:43
	% EndTime: 2019-10-24 10:37:43
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (93->14), mult. (54->18), div. (0->0), fcn. (93->6), ass. (0->16)
	t69 = qJ(3) + qJ(4) + qJ(5);
	t67 = sin(t69);
	t72 = sin(qJ(2));
	t77 = t72 * t67;
	t68 = cos(t69);
	t76 = t72 * t68;
	t73 = cos(qJ(2));
	t75 = t73 * t67;
	t74 = t73 * t68;
	t71 = cos(pkin(9));
	t70 = sin(pkin(9));
	t66 = -t70 * t67 - t71 * t74;
	t65 = t70 * t68 - t71 * t75;
	t64 = t71 * t67 - t70 * t74;
	t63 = -t71 * t68 - t70 * t75;
	t1 = [0, -t71 * t76, t65, t65, t65; 0, -t70 * t76, t63, t63, t63; 0, t74, -t77, -t77, -t77; 0, t71 * t77, t66, t66, t66; 0, t70 * t77, t64, t64, t64; 0, -t75, -t76, -t76, -t76; 0, t71 * t73, 0, 0, 0; 0, t70 * t73, 0, 0, 0; 0, t72, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end