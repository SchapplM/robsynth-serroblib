% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RRRRRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:02
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:02
	% EndTime: 2019-10-10 13:18:03
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (741->21), mult. (440->54), div. (106->9), fcn. (646->9), ass. (0->38)
	t74 = qJ(2) + qJ(3) + qJ(4);
	t73 = cos(t74);
	t72 = sin(t74);
	t76 = sin(qJ(1));
	t84 = t76 * t72;
	t67 = atan2(-t84, -t73);
	t61 = sin(t67);
	t62 = cos(t67);
	t58 = -t61 * t84 - t62 * t73;
	t57 = 0.1e1 / t58 ^ 2;
	t78 = cos(qJ(1));
	t90 = t57 * t78 ^ 2;
	t89 = t61 * t73;
	t77 = cos(qJ(5));
	t80 = t78 * t77;
	t75 = sin(qJ(5));
	t83 = t76 * t75;
	t66 = t73 * t80 + t83;
	t64 = 0.1e1 / t66 ^ 2;
	t81 = t78 * t75;
	t82 = t76 * t77;
	t65 = t73 * t81 - t82;
	t88 = t64 * t65;
	t69 = t72 ^ 2;
	t87 = t69 / t73 ^ 2;
	t86 = t72 * t78;
	t68 = 0.1e1 / (t76 ^ 2 * t87 + 0.1e1);
	t85 = t76 * t68;
	t79 = t65 ^ 2 * t64 + 0.1e1;
	t70 = 0.1e1 / t73;
	t63 = 0.1e1 / t66;
	t60 = 0.1e1 / t79;
	t59 = (0.1e1 + t87) * t85;
	t56 = 0.1e1 / t58;
	t55 = 0.1e1 / (t69 * t90 + 0.1e1);
	t54 = (-t63 * t75 + t77 * t88) * t60 * t86;
	t53 = (t73 * t56 - (-t76 * t89 + t62 * t72 + (-t62 * t84 + t89) * t59) * t72 * t57) * t78 * t55;
	t1 = [t70 * t68 * t86, t59, t59, t59, 0, 0; (-t56 * t84 - (-t62 * t69 * t70 * t85 + (t68 - 0.1e1) * t72 * t61) * t72 * t90) * t55, t53, t53, t53, 0, 0; ((-t73 * t83 - t80) * t63 - (-t73 * t82 + t81) * t88) * t60, t54, t54, t54, t79 * t60, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:18:03
	% EndTime: 2019-10-10 13:18:03
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (859->22), mult. (467->54), div. (111->9), fcn. (681->9), ass. (0->40)
	t85 = qJ(2) + qJ(3) + qJ(4);
	t82 = cos(t85);
	t81 = sin(t85);
	t87 = sin(qJ(1));
	t94 = t87 * t81;
	t76 = atan2(-t94, -t82);
	t74 = sin(t76);
	t75 = cos(t76);
	t67 = -t74 * t94 - t75 * t82;
	t66 = 0.1e1 / t67 ^ 2;
	t88 = cos(qJ(1));
	t100 = t66 * t88 ^ 2;
	t86 = qJ(5) + qJ(6);
	t84 = cos(t86);
	t90 = t88 * t84;
	t83 = sin(t86);
	t93 = t87 * t83;
	t73 = t82 * t90 + t93;
	t71 = 0.1e1 / t73 ^ 2;
	t91 = t88 * t83;
	t92 = t87 * t84;
	t72 = t82 * t91 - t92;
	t99 = t71 * t72;
	t98 = t74 * t82;
	t78 = t81 ^ 2;
	t97 = t78 / t82 ^ 2;
	t96 = t81 * t88;
	t77 = 0.1e1 / (t87 ^ 2 * t97 + 0.1e1);
	t95 = t87 * t77;
	t89 = t72 ^ 2 * t71 + 0.1e1;
	t79 = 0.1e1 / t82;
	t70 = 0.1e1 / t73;
	t69 = 0.1e1 / t89;
	t68 = (0.1e1 + t97) * t95;
	t65 = 0.1e1 / t67;
	t64 = 0.1e1 / (t78 * t100 + 0.1e1);
	t63 = t89 * t69;
	t62 = (-t70 * t83 + t84 * t99) * t69 * t96;
	t61 = (t82 * t65 - (-t87 * t98 + t75 * t81 + (-t75 * t94 + t98) * t68) * t81 * t66) * t88 * t64;
	t1 = [t79 * t77 * t96, t68, t68, t68, 0, 0; (-t65 * t94 - (-t75 * t78 * t79 * t95 + (t77 - 0.1e1) * t81 * t74) * t81 * t100) * t64, t61, t61, t61, 0, 0; ((-t82 * t93 - t90) * t70 - (-t82 * t92 + t91) * t99) * t69, t62, t62, t62, t63, t63;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end