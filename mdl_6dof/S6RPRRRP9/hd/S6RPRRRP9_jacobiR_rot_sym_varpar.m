% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
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
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; -t6, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t6, 0, 0, 0, 0, 0; t5, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (5->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t11 = sin(qJ(3));
	t12 = sin(qJ(1));
	t16 = t12 * t11;
	t13 = cos(qJ(3));
	t14 = cos(qJ(1));
	t15 = t14 * t13;
	t10 = t14 * t11;
	t9 = t12 * t13;
	t1 = [t10, 0, t9, 0, 0, 0; t16, 0, -t15, 0, 0, 0; 0, 0, -t11, 0, 0, 0; t15, 0, -t16, 0, 0, 0; t9, 0, t10, 0, 0, 0; 0, 0, -t13, 0, 0, 0; -t12, 0, 0, 0, 0, 0; t14, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (16->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t65 = sin(qJ(3));
	t66 = sin(qJ(1));
	t76 = t66 * t65;
	t67 = cos(qJ(4));
	t75 = t66 * t67;
	t64 = sin(qJ(4));
	t68 = cos(qJ(3));
	t74 = t68 * t64;
	t73 = t68 * t67;
	t69 = cos(qJ(1));
	t72 = t69 * t65;
	t71 = t69 * t67;
	t70 = t69 * t68;
	t63 = -t66 * t64 + t65 * t71;
	t62 = t64 * t72 + t75;
	t61 = t69 * t64 + t65 * t75;
	t60 = -t64 * t76 + t71;
	t1 = [t63, 0, t66 * t73, t60, 0, 0; t61, 0, -t67 * t70, t62, 0, 0; 0, 0, -t65 * t67, -t74, 0, 0; -t62, 0, -t66 * t74, -t61, 0, 0; t60, 0, t64 * t70, t63, 0, 0; 0, 0, t65 * t64, -t73, 0, 0; -t70, 0, t76, 0, 0, 0; -t66 * t68, 0, -t72, 0, 0, 0; 0, 0, t68, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (56->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
	t83 = sin(qJ(3));
	t84 = sin(qJ(1));
	t91 = t84 * t83;
	t82 = qJ(4) + qJ(5);
	t80 = sin(t82);
	t85 = cos(qJ(3));
	t90 = t85 * t80;
	t81 = cos(t82);
	t89 = t85 * t81;
	t86 = cos(qJ(1));
	t88 = t86 * t83;
	t87 = t86 * t85;
	t79 = -t84 * t80 + t81 * t88;
	t78 = t80 * t88 + t84 * t81;
	t77 = t86 * t80 + t81 * t91;
	t76 = -t80 * t91 + t86 * t81;
	t1 = [t79, 0, t84 * t89, t76, t76, 0; t77, 0, -t81 * t87, t78, t78, 0; 0, 0, -t83 * t81, -t90, -t90, 0; -t78, 0, -t84 * t90, -t77, -t77, 0; t76, 0, t80 * t87, t79, t79, 0; 0, 0, t83 * t80, -t89, -t89, 0; -t87, 0, t91, 0, 0, 0; -t84 * t85, 0, -t88, 0, 0, 0; 0, 0, t85, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (56->18), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
	t85 = sin(qJ(3));
	t86 = sin(qJ(1));
	t93 = t86 * t85;
	t84 = qJ(4) + qJ(5);
	t82 = sin(t84);
	t87 = cos(qJ(3));
	t92 = t87 * t82;
	t83 = cos(t84);
	t91 = t87 * t83;
	t88 = cos(qJ(1));
	t90 = t88 * t85;
	t89 = t88 * t87;
	t81 = -t86 * t82 + t83 * t90;
	t80 = t82 * t90 + t86 * t83;
	t79 = t88 * t82 + t83 * t93;
	t78 = -t82 * t93 + t88 * t83;
	t1 = [t81, 0, t86 * t91, t78, t78, 0; t79, 0, -t83 * t89, t80, t80, 0; 0, 0, -t85 * t83, -t92, -t92, 0; -t80, 0, -t86 * t92, -t79, -t79, 0; t78, 0, t82 * t89, t81, t81, 0; 0, 0, t85 * t82, -t91, -t91, 0; -t89, 0, t93, 0, 0, 0; -t86 * t87, 0, -t90, 0, 0, 0; 0, 0, t87, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end