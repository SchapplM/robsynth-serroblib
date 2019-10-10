% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:23
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:23
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
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
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
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
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t41 = sin(qJ(2));
	t42 = sin(qJ(1));
	t46 = t42 * t41;
	t43 = cos(qJ(2));
	t44 = cos(qJ(1));
	t45 = t44 * t43;
	t40 = t44 * t41;
	t39 = t42 * t43;
	t1 = [t44, 0, 0, 0, 0, 0; t42, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t39, t40, 0, 0, 0, 0; -t45, t46, 0, 0, 0, 0; 0, -t43, 0, 0, 0, 0; -t46, t45, 0, 0, 0, 0; t40, t39, 0, 0, 0, 0; 0, t41, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->11)
	t49 = sin(qJ(2));
	t50 = sin(qJ(1));
	t56 = t50 * t49;
	t51 = cos(qJ(2));
	t55 = t50 * t51;
	t52 = cos(qJ(1));
	t54 = t52 * t49;
	t53 = t52 * t51;
	t48 = cos(pkin(9));
	t47 = sin(pkin(9));
	t1 = [-t47 * t56 + t52 * t48, t47 * t53, 0, 0, 0, 0; t47 * t54 + t50 * t48, t47 * t55, 0, 0, 0, 0; 0, t49 * t47, 0, 0, 0, 0; -t52 * t47 - t48 * t56, t48 * t53, 0, 0, 0, 0; -t50 * t47 + t48 * t54, t48 * t55, 0, 0, 0, 0; 0, t49 * t48, 0, 0, 0, 0; -t55, -t54, 0, 0, 0, 0; t53, -t56, 0, 0, 0, 0; 0, t51, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (10->10), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->11)
	t63 = sin(qJ(2));
	t64 = sin(qJ(1));
	t70 = t64 * t63;
	t65 = cos(qJ(2));
	t69 = t64 * t65;
	t66 = cos(qJ(1));
	t68 = t66 * t63;
	t67 = t66 * t65;
	t62 = cos(pkin(9));
	t61 = sin(pkin(9));
	t1 = [-t61 * t70 + t66 * t62, t61 * t67, 0, 0, 0, 0; t61 * t68 + t64 * t62, t61 * t69, 0, 0, 0, 0; 0, t63 * t61, 0, 0, 0, 0; -t69, -t68, 0, 0, 0, 0; t67, -t70, 0, 0, 0, 0; 0, t65, 0, 0, 0, 0; t66 * t61 + t62 * t70, -t62 * t67, 0, 0, 0, 0; t64 * t61 - t62 * t68, -t62 * t69, 0, 0, 0, 0; 0, -t63 * t62, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (34->17), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->23)
	t85 = cos(qJ(2));
	t79 = sin(pkin(9));
	t80 = cos(pkin(9));
	t81 = sin(qJ(6));
	t84 = cos(qJ(6));
	t89 = t79 * t81 + t80 * t84;
	t87 = t89 * t85;
	t82 = sin(qJ(2));
	t83 = sin(qJ(1));
	t94 = t83 * t82;
	t86 = cos(qJ(1));
	t93 = t86 * t82;
	t76 = t86 * t79 + t80 * t94;
	t77 = -t79 * t94 + t86 * t80;
	t92 = t76 * t84 - t77 * t81;
	t91 = t76 * t81 + t77 * t84;
	t90 = t79 * t84 - t80 * t81;
	t88 = t90 * t85;
	t75 = t79 * t93 + t83 * t80;
	t74 = t83 * t79 - t80 * t93;
	t73 = t74 * t81 + t75 * t84;
	t72 = t74 * t84 - t75 * t81;
	t1 = [t91, t86 * t88, 0, 0, 0, t72; t73, t83 * t88, 0, 0, 0, -t92; 0, t90 * t82, 0, 0, 0, t87; t92, -t86 * t87, 0, 0, 0, -t73; t72, -t83 * t87, 0, 0, 0, t91; 0, -t89 * t82, 0, 0, 0, t88; t83 * t85, t93, 0, 0, 0, 0; -t86 * t85, t94, 0, 0, 0, 0; 0, -t85, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end