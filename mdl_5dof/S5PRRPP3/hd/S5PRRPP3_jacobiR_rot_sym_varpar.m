% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRPP3
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRPP3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:24
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
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:24
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
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (17->13), mult. (58->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t72 = sin(qJ(3));
	t73 = sin(qJ(2));
	t83 = t72 * t73;
	t68 = sin(pkin(8));
	t82 = t73 * t68;
	t70 = cos(pkin(8));
	t81 = t73 * t70;
	t74 = cos(qJ(3));
	t80 = t73 * t74;
	t75 = cos(qJ(2));
	t79 = t74 * t75;
	t78 = t75 * t72;
	t77 = t68 * t75 - t70 * t80;
	t76 = t68 * t80 + t70 * t75;
	t71 = cos(pkin(7));
	t69 = sin(pkin(7));
	t67 = t69 * t74 - t71 * t78;
	t66 = -t69 * t78 - t71 * t74;
	t1 = [0, t77 * t71, t67 * t70, 0, 0; 0, t77 * t69, t66 * t70, 0, 0; 0, t70 * t79 + t82, -t72 * t81, 0, 0; 0, t76 * t71, -t67 * t68, 0, 0; 0, t76 * t69, -t66 * t68, 0, 0; 0, -t68 * t79 + t81, t72 * t82, 0, 0; 0, -t71 * t83, t69 * t72 + t71 * t79, 0, 0; 0, -t69 * t83, t69 * t79 - t71 * t72, 0, 0; 0, t78, t80, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:14:24
	% EndTime: 2019-12-05 16:14:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (16->12), mult. (58->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t88 = sin(qJ(3));
	t89 = sin(qJ(2));
	t99 = t88 * t89;
	t84 = sin(pkin(8));
	t98 = t89 * t84;
	t86 = cos(pkin(8));
	t97 = t89 * t86;
	t90 = cos(qJ(3));
	t96 = t89 * t90;
	t91 = cos(qJ(2));
	t95 = t90 * t91;
	t94 = t91 * t88;
	t93 = t84 * t91 - t86 * t96;
	t92 = -t84 * t96 - t86 * t91;
	t87 = cos(pkin(7));
	t85 = sin(pkin(7));
	t83 = t85 * t90 - t87 * t94;
	t82 = -t85 * t94 - t87 * t90;
	t1 = [0, t93 * t87, t83 * t86, 0, 0; 0, t93 * t85, t82 * t86, 0, 0; 0, t86 * t95 + t98, -t88 * t97, 0, 0; 0, -t87 * t99, t85 * t88 + t87 * t95, 0, 0; 0, -t85 * t99, t85 * t95 - t87 * t88, 0, 0; 0, t94, t96, 0, 0; 0, t92 * t87, t83 * t84, 0, 0; 0, t92 * t85, t82 * t84, 0, 0; 0, t84 * t95 - t97, -t88 * t98, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end