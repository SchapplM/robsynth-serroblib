% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:25
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0, 0; t13, -t16, 0, 0, 0, 0; 0, t11, 0, 0, 0, 0; t16, -t13, 0, 0, 0, 0; -t15, -t14, 0, 0, 0, 0; 0, -t9, 0, 0, 0, 0; t12, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->8), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t47 = sin(qJ(2));
	t48 = sin(qJ(1));
	t55 = t48 * t47;
	t45 = sin(pkin(9));
	t49 = cos(qJ(2));
	t54 = t49 * t45;
	t46 = cos(pkin(9));
	t53 = t49 * t46;
	t50 = cos(qJ(1));
	t52 = t50 * t47;
	t51 = t50 * t49;
	t1 = [t50 * t45 - t48 * t53, -t46 * t52, 0, 0, 0, 0; t48 * t45 + t46 * t51, -t46 * t55, 0, 0, 0, 0; 0, t53, 0, 0, 0, 0; t50 * t46 + t48 * t54, t45 * t52, 0, 0, 0, 0; -t45 * t51 + t48 * t46, t45 * t55, 0, 0, 0, 0; 0, -t54, 0, 0, 0, 0; -t55, t51, 0, 0, 0, 0; t52, t48 * t49, 0, 0, 0, 0; 0, t47, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->8), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t59 = sin(qJ(2));
	t60 = sin(qJ(1));
	t67 = t60 * t59;
	t57 = sin(pkin(9));
	t61 = cos(qJ(2));
	t66 = t61 * t57;
	t58 = cos(pkin(9));
	t65 = t61 * t58;
	t62 = cos(qJ(1));
	t64 = t62 * t59;
	t63 = t62 * t61;
	t1 = [-t67, t63, 0, 0, 0, 0; t64, t60 * t61, 0, 0, 0, 0; 0, t59, 0, 0, 0, 0; -t62 * t57 + t60 * t65, t58 * t64, 0, 0, 0, 0; -t60 * t57 - t58 * t63, t58 * t67, 0, 0, 0, 0; 0, -t65, 0, 0, 0, 0; -t62 * t58 - t60 * t66, -t57 * t64, 0, 0, 0, 0; t57 * t63 - t60 * t58, -t57 * t67, 0, 0, 0, 0; 0, t66, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (12->12), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t63 = sin(qJ(2));
	t64 = sin(qJ(1));
	t71 = t64 * t63;
	t61 = sin(pkin(9));
	t65 = cos(qJ(2));
	t70 = t65 * t61;
	t62 = cos(pkin(9));
	t69 = t65 * t62;
	t66 = cos(qJ(1));
	t68 = t66 * t63;
	t67 = t66 * t65;
	t1 = [-t66 * t62 - t64 * t70, -t61 * t68, 0, 0, 0, 0; t61 * t67 - t64 * t62, -t61 * t71, 0, 0, 0, 0; 0, t70, 0, 0, 0, 0; t71, -t67, 0, 0, 0, 0; -t68, -t64 * t65, 0, 0, 0, 0; 0, -t63, 0, 0, 0, 0; t66 * t61 - t64 * t69, -t62 * t68, 0, 0, 0, 0; t64 * t61 + t62 * t67, -t62 * t71, 0, 0, 0, 0; 0, t69, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->18), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->23)
	t80 = sin(qJ(2));
	t77 = sin(pkin(9));
	t78 = cos(pkin(9));
	t79 = sin(qJ(6));
	t82 = cos(qJ(6));
	t87 = t77 * t79 - t78 * t82;
	t86 = t80 * t87;
	t81 = sin(qJ(1));
	t83 = cos(qJ(2));
	t92 = t81 * t83;
	t84 = cos(qJ(1));
	t91 = t84 * t83;
	t72 = t77 * t92 + t84 * t78;
	t73 = -t84 * t77 + t78 * t92;
	t90 = -t72 * t82 - t73 * t79;
	t89 = t72 * t79 - t73 * t82;
	t88 = t77 * t82 + t78 * t79;
	t85 = t88 * t80;
	t75 = t81 * t77 + t78 * t91;
	t74 = t77 * t91 - t81 * t78;
	t71 = t74 * t82 + t75 * t79;
	t70 = -t74 * t79 + t75 * t82;
	t1 = [t90, -t84 * t85, 0, 0, 0, t70; t71, -t81 * t85, 0, 0, 0, -t89; 0, t88 * t83, 0, 0, 0, -t86; t89, t84 * t86, 0, 0, 0, -t71; t70, t81 * t86, 0, 0, 0, t90; 0, -t87 * t83, 0, 0, 0, -t85; -t81 * t80, t91, 0, 0, 0, 0; t84 * t80, t92, 0, 0, 0, 0; 0, t80, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end