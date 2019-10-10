% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:52
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR10_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR10_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
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
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
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
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
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
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
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
	t48 = cos(pkin(10));
	t47 = sin(pkin(10));
	t1 = [-t47 * t56 + t52 * t48, t47 * t53, 0, 0, 0, 0; t47 * t54 + t50 * t48, t47 * t55, 0, 0, 0, 0; 0, t49 * t47, 0, 0, 0, 0; -t52 * t47 - t48 * t56, t48 * t53, 0, 0, 0, 0; -t50 * t47 + t48 * t54, t48 * t55, 0, 0, 0, 0; 0, t49 * t48, 0, 0, 0, 0; -t55, -t54, 0, 0, 0, 0; t53, -t56, 0, 0, 0, 0; 0, t51, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (36->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
	t71 = sin(qJ(2));
	t72 = sin(qJ(1));
	t79 = t72 * t71;
	t70 = pkin(10) + qJ(5);
	t68 = sin(t70);
	t73 = cos(qJ(2));
	t78 = t73 * t68;
	t69 = cos(t70);
	t77 = t73 * t69;
	t74 = cos(qJ(1));
	t76 = t74 * t71;
	t75 = t74 * t73;
	t67 = -t68 * t79 + t74 * t69;
	t66 = t74 * t68 + t69 * t79;
	t65 = t68 * t76 + t72 * t69;
	t64 = -t72 * t68 + t69 * t76;
	t1 = [t67, t68 * t75, 0, 0, t64, 0; t65, t72 * t78, 0, 0, t66, 0; 0, t71 * t68, 0, 0, -t77, 0; -t66, t69 * t75, 0, 0, -t65, 0; t64, t72 * t77, 0, 0, t67, 0; 0, t71 * t69, 0, 0, t78, 0; -t72 * t73, -t76, 0, 0, 0, 0; t75, -t79, 0, 0, 0, 0; 0, t73, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (85->13), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
	t88 = sin(qJ(2));
	t89 = sin(qJ(1));
	t95 = t89 * t88;
	t87 = pkin(10) + qJ(5) + qJ(6);
	t85 = sin(t87);
	t90 = cos(qJ(2));
	t84 = t90 * t85;
	t86 = cos(t87);
	t94 = t90 * t86;
	t91 = cos(qJ(1));
	t93 = t91 * t88;
	t92 = t91 * t90;
	t83 = -t85 * t95 + t91 * t86;
	t82 = t91 * t85 + t86 * t95;
	t81 = t85 * t93 + t89 * t86;
	t80 = -t89 * t85 + t86 * t93;
	t1 = [t83, t85 * t92, 0, 0, t80, t80; t81, t89 * t84, 0, 0, t82, t82; 0, t88 * t85, 0, 0, -t94, -t94; -t82, t86 * t92, 0, 0, -t81, -t81; t80, t89 * t94, 0, 0, t83, t83; 0, t88 * t86, 0, 0, t84, t84; -t89 * t90, -t93, 0, 0, 0, 0; t92, -t95, 0, 0, 0, 0; 0, t90, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end