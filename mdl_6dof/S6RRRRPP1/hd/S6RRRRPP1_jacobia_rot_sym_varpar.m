% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP1
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
%   Wie in S6RRRRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
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
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (422->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
	t72 = qJ(2) + qJ(3);
	t70 = cos(t72);
	t69 = sin(t72);
	t73 = sin(qJ(1));
	t78 = t73 * t69;
	t62 = atan2(-t78, -t70);
	t60 = sin(t62);
	t61 = cos(t62);
	t53 = -t60 * t78 - t61 * t70;
	t52 = 0.1e1 / t53 ^ 2;
	t74 = cos(qJ(1));
	t86 = t52 * t74 ^ 2;
	t71 = qJ(4) + pkin(10);
	t68 = cos(t71);
	t76 = t74 * t68;
	t67 = sin(t71);
	t80 = t73 * t67;
	t59 = t70 * t76 + t80;
	t57 = 0.1e1 / t59 ^ 2;
	t77 = t74 * t67;
	t79 = t73 * t68;
	t58 = t70 * t77 - t79;
	t85 = t57 * t58;
	t84 = t60 * t70;
	t64 = t69 ^ 2;
	t83 = t64 / t70 ^ 2;
	t82 = t69 * t74;
	t63 = 0.1e1 / (t73 ^ 2 * t83 + 0.1e1);
	t81 = t73 * t63;
	t75 = t58 ^ 2 * t57 + 0.1e1;
	t65 = 0.1e1 / t70;
	t56 = 0.1e1 / t59;
	t55 = (0.1e1 + t83) * t81;
	t54 = 0.1e1 / t75;
	t51 = 0.1e1 / t53;
	t50 = 0.1e1 / (t64 * t86 + 0.1e1);
	t49 = (-t56 * t67 + t68 * t85) * t54 * t82;
	t48 = (t70 * t51 - (-t73 * t84 + t61 * t69 + (-t61 * t78 + t84) * t55) * t69 * t52) * t74 * t50;
	t1 = [t65 * t63 * t82, t55, t55, 0, 0, 0; (-t51 * t78 - (-t61 * t64 * t65 * t81 + (t63 - 0.1e1) * t69 * t60) * t69 * t86) * t50, t48, t48, 0, 0, 0; ((-t70 * t80 - t76) * t56 - (-t70 * t79 + t77) * t85) * t54, t49, t49, t75 * t54, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:20:28
	% EndTime: 2019-10-10 12:20:28
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (860->28), mult. (688->69), div. (136->11), fcn. (1034->9), ass. (0->44)
	t85 = qJ(2) + qJ(3);
	t81 = sin(t85);
	t100 = t81 ^ 2;
	t82 = cos(t85);
	t83 = qJ(4) + pkin(10);
	t80 = cos(t83);
	t87 = cos(qJ(1));
	t89 = t87 * t80;
	t79 = sin(t83);
	t86 = sin(qJ(1));
	t92 = t86 * t79;
	t67 = t82 * t92 + t89;
	t94 = t81 * t79;
	t63 = atan2(-t67, t94);
	t60 = sin(t63);
	t61 = cos(t63);
	t59 = -t60 * t67 + t61 * t94;
	t58 = 0.1e1 / t59 ^ 2;
	t90 = t87 * t79;
	t91 = t86 * t80;
	t70 = t82 * t90 - t91;
	t99 = t58 * t70;
	t97 = t61 * t67;
	t96 = t70 ^ 2 * t58;
	t74 = 0.1e1 / t79;
	t77 = 0.1e1 / t81;
	t95 = t74 * t77;
	t93 = t81 * t87;
	t71 = t82 * t89 + t92;
	t66 = 0.1e1 / t71 ^ 2;
	t88 = t87 ^ 2 * t100 * t66;
	t78 = 0.1e1 / t100;
	t75 = 0.1e1 / t79 ^ 2;
	t69 = t82 * t91 - t90;
	t65 = 0.1e1 / t71;
	t64 = 0.1e1 / (0.1e1 + t88);
	t62 = 0.1e1 / (t67 ^ 2 * t78 * t75 + 0.1e1);
	t57 = 0.1e1 / t59;
	t56 = (t67 * t74 * t78 * t82 + t86) * t62;
	t55 = (-t65 * t82 * t87 - t80 * t88) * t64;
	t54 = 0.1e1 / (0.1e1 + t96);
	t53 = (t67 * t75 * t80 - t69 * t74) * t77 * t62;
	t52 = (t56 * t97 * t99 + (-t57 * t93 - (t61 * t82 + (-t56 + t86) * t81 * t60) * t99) * t79) * t54;
	t1 = [-t70 * t62 * t95, t56, t56, t53, 0, 0; (-t67 * t57 - (-t60 + (t95 * t97 + t60) * t62) * t96) * t54, t52, t52, (t71 * t57 - (t61 * t81 * t80 - t60 * t69 + (-t60 * t94 - t97) * t53) * t99) * t54, 0, 0; (-t66 * t69 * t87 + t65 * t86) * t81 * t64, t55, t55, -t70 * t66 * t64 * t93, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end