% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP11_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP11_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
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
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
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
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
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
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (12->10), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t65 = sin(qJ(2));
	t66 = sin(qJ(1));
	t76 = t66 * t65;
	t67 = cos(qJ(4));
	t75 = t66 * t67;
	t64 = sin(qJ(4));
	t68 = cos(qJ(2));
	t74 = t68 * t64;
	t73 = t68 * t67;
	t69 = cos(qJ(1));
	t72 = t69 * t65;
	t71 = t69 * t67;
	t70 = t69 * t68;
	t63 = -t64 * t76 + t71;
	t62 = t69 * t64 + t65 * t75;
	t61 = t64 * t72 + t75;
	t60 = -t66 * t64 + t65 * t71;
	t1 = [t63, t64 * t70, 0, t60, 0, 0; t61, t66 * t74, 0, t62, 0, 0; 0, t65 * t64, 0, -t73, 0, 0; -t62, t67 * t70, 0, -t61, 0, 0; t60, t66 * t73, 0, t63, 0, 0; 0, t65 * t67, 0, t74, 0, 0; -t66 * t68, -t72, 0, 0, 0, 0; t70, -t76, 0, 0, 0, 0; 0, t68, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (51->13), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
	t86 = sin(qJ(2));
	t87 = sin(qJ(1));
	t93 = t87 * t86;
	t85 = qJ(4) + qJ(5);
	t83 = sin(t85);
	t88 = cos(qJ(2));
	t82 = t88 * t83;
	t84 = cos(t85);
	t92 = t88 * t84;
	t89 = cos(qJ(1));
	t91 = t89 * t86;
	t90 = t89 * t88;
	t81 = -t83 * t93 + t89 * t84;
	t80 = t89 * t83 + t84 * t93;
	t79 = t83 * t91 + t87 * t84;
	t78 = -t87 * t83 + t84 * t91;
	t1 = [t81, t83 * t90, 0, t78, t78, 0; t79, t87 * t82, 0, t80, t80, 0; 0, t86 * t83, 0, -t92, -t92, 0; -t80, t84 * t90, 0, -t79, -t79, 0; t78, t87 * t92, 0, t81, t81, 0; 0, t86 * t84, 0, t82, t82, 0; -t87 * t88, -t91, 0, 0, 0, 0; t90, -t93, 0, 0, 0, 0; 0, t88, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (51->13), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
	t86 = sin(qJ(2));
	t87 = sin(qJ(1));
	t93 = t87 * t86;
	t85 = qJ(4) + qJ(5);
	t83 = sin(t85);
	t88 = cos(qJ(2));
	t82 = t88 * t83;
	t84 = cos(t85);
	t92 = t88 * t84;
	t89 = cos(qJ(1));
	t91 = t89 * t86;
	t90 = t89 * t88;
	t81 = -t83 * t93 + t89 * t84;
	t80 = t89 * t83 + t84 * t93;
	t79 = t83 * t91 + t87 * t84;
	t78 = -t87 * t83 + t84 * t91;
	t1 = [t81, t83 * t90, 0, t78, t78, 0; t79, t87 * t82, 0, t80, t80, 0; 0, t86 * t83, 0, -t92, -t92, 0; -t80, t84 * t90, 0, -t79, -t79, 0; t78, t87 * t92, 0, t81, t81, 0; 0, t86 * t84, 0, t82, t82, 0; -t87 * t88, -t91, 0, 0, 0, 0; t90, -t93, 0, 0, 0, 0; 0, t88, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end