% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP3
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
%   Wie in S6RRRRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (358->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t68 = qJ(2) + qJ(3);
	t67 = cos(t68);
	t66 = sin(t68);
	t70 = sin(qJ(1));
	t78 = t70 * t66;
	t61 = atan2(-t78, -t67);
	t59 = sin(t61);
	t60 = cos(t61);
	t52 = -t59 * t78 - t60 * t67;
	t51 = 0.1e1 / t52 ^ 2;
	t72 = cos(qJ(1));
	t84 = t51 * t72 ^ 2;
	t71 = cos(qJ(4));
	t74 = t72 * t71;
	t69 = sin(qJ(4));
	t77 = t70 * t69;
	t58 = t67 * t74 + t77;
	t56 = 0.1e1 / t58 ^ 2;
	t75 = t72 * t69;
	t76 = t70 * t71;
	t57 = t67 * t75 - t76;
	t83 = t56 * t57;
	t82 = t59 * t67;
	t63 = t66 ^ 2;
	t81 = t63 / t67 ^ 2;
	t80 = t66 * t72;
	t62 = 0.1e1 / (t70 ^ 2 * t81 + 0.1e1);
	t79 = t70 * t62;
	t73 = t57 ^ 2 * t56 + 0.1e1;
	t64 = 0.1e1 / t67;
	t55 = 0.1e1 / t58;
	t54 = 0.1e1 / t73;
	t53 = (0.1e1 + t81) * t79;
	t50 = 0.1e1 / t52;
	t49 = 0.1e1 / (t63 * t84 + 0.1e1);
	t48 = (-t55 * t69 + t71 * t83) * t54 * t80;
	t47 = (t67 * t50 - (-t70 * t82 + t60 * t66 + (-t60 * t78 + t82) * t53) * t66 * t51) * t72 * t49;
	t1 = [t64 * t62 * t80, t53, t53, 0, 0, 0; (-t50 * t78 - (-t60 * t63 * t64 * t79 + (t62 - 0.1e1) * t66 * t59) * t66 * t84) * t49, t47, t47, 0, 0, 0; ((-t67 * t77 - t74) * t55 - (-t67 * t76 + t75) * t83) * t54, t48, t48, t73 * t54, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (469->27), mult. (656->69), div. (146->11), fcn. (1013->9), ass. (0->42)
	t82 = qJ(2) + qJ(3);
	t77 = cos(t82);
	t85 = cos(qJ(4));
	t86 = cos(qJ(1));
	t87 = t86 * t85;
	t83 = sin(qJ(4));
	t84 = sin(qJ(1));
	t90 = t84 * t83;
	t67 = t77 * t90 + t87;
	t76 = sin(t82);
	t91 = t76 * t83;
	t66 = atan2(-t67, t91);
	t62 = sin(t66);
	t63 = cos(t66);
	t61 = -t62 * t67 + t63 * t91;
	t60 = 0.1e1 / t61 ^ 2;
	t88 = t86 * t83;
	t89 = t84 * t85;
	t70 = t77 * t88 - t89;
	t98 = t60 * t70;
	t96 = t63 * t67;
	t71 = t77 * t87 + t90;
	t75 = 0.1e1 / t76 ^ 2;
	t81 = 0.1e1 / t86 ^ 2;
	t65 = 0.1e1 / (t71 ^ 2 * t75 * t81 + 0.1e1);
	t74 = 0.1e1 / t76;
	t95 = t65 * t74;
	t94 = t70 ^ 2 * t60;
	t78 = 0.1e1 / t83;
	t93 = t74 * t78;
	t92 = t75 * t77;
	t80 = 0.1e1 / t86;
	t79 = 0.1e1 / t83 ^ 2;
	t69 = t77 * t89 - t88;
	t64 = 0.1e1 / (t67 ^ 2 * t75 * t79 + 0.1e1);
	t59 = 0.1e1 / t61;
	t58 = (-t71 * t80 * t92 - t85) * t65;
	t57 = (t67 * t78 * t92 + t84) * t64;
	t56 = 0.1e1 / (0.1e1 + t94);
	t55 = (t67 * t79 * t85 - t69 * t78) * t74 * t64;
	t54 = (t57 * t96 * t98 + (-t86 * t76 * t59 - (t63 * t77 + (-t57 + t84) * t76 * t62) * t98) * t83) * t56;
	t1 = [-t70 * t64 * t93, t57, t57, t55, 0, 0; (-t67 * t59 - (-t62 + (t93 * t96 + t62) * t64) * t94) * t56, t54, t54, (t71 * t59 - (t63 * t76 * t85 - t62 * t69 + (-t62 * t91 - t96) * t55) * t98) * t56, 0, 0; (t71 * t81 * t84 - t69 * t80) * t95, t58, t58, -t70 * t80 * t95, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:24:00
	% EndTime: 2019-10-10 12:24:00
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (469->27), mult. (656->69), div. (146->11), fcn. (1013->9), ass. (0->42)
	t84 = qJ(2) + qJ(3);
	t79 = cos(t84);
	t85 = sin(qJ(4));
	t88 = cos(qJ(1));
	t90 = t88 * t85;
	t86 = sin(qJ(1));
	t87 = cos(qJ(4));
	t91 = t86 * t87;
	t71 = t79 * t91 - t90;
	t78 = sin(t84);
	t93 = t78 * t87;
	t69 = atan2(-t71, t93);
	t65 = sin(t69);
	t66 = cos(t69);
	t64 = -t65 * t71 + t66 * t93;
	t63 = 0.1e1 / t64 ^ 2;
	t89 = t88 * t87;
	t92 = t86 * t85;
	t74 = t79 * t89 + t92;
	t100 = t63 * t74;
	t98 = t66 * t71;
	t73 = -t79 * t90 + t91;
	t77 = 0.1e1 / t78 ^ 2;
	t83 = 0.1e1 / t88 ^ 2;
	t68 = 0.1e1 / (t73 ^ 2 * t83 * t77 + 0.1e1);
	t76 = 0.1e1 / t78;
	t97 = t68 * t76;
	t96 = t74 ^ 2 * t63;
	t80 = 0.1e1 / t87;
	t95 = t76 * t80;
	t94 = t77 * t79;
	t82 = 0.1e1 / t88;
	t81 = 0.1e1 / t87 ^ 2;
	t70 = t79 * t92 + t89;
	t67 = 0.1e1 / (t71 ^ 2 * t77 * t81 + 0.1e1);
	t62 = 0.1e1 / t64;
	t61 = (t71 * t80 * t94 + t86) * t67;
	t60 = (-t73 * t82 * t94 + t85) * t68;
	t59 = 0.1e1 / (0.1e1 + t96);
	t58 = (-t71 * t81 * t85 + t70 * t80) * t76 * t67;
	t57 = (t61 * t98 * t100 + (-t88 * t78 * t62 - (t66 * t79 + (-t61 + t86) * t78 * t65) * t100) * t87) * t59;
	t1 = [-t74 * t67 * t95, t61, t61, t58, 0, 0; (-t71 * t62 - (-t65 + (t95 * t98 + t65) * t67) * t96) * t59, t57, t57, (t73 * t62 - (-t66 * t78 * t85 + t65 * t70 + (-t65 * t93 - t98) * t58) * t100) * t59, 0, 0; (t73 * t83 * t86 + t70 * t82) * t97, t60, t60, -t74 * t82 * t97, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end