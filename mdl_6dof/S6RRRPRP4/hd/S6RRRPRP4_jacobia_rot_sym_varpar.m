% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP4
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
%   Wie in S6RRRPRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:40
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (297->18), mult. (231->39), div. (67->9), fcn. (349->7), ass. (0->30)
	t55 = cos(qJ(1));
	t52 = t55 ^ 2;
	t53 = qJ(2) + qJ(3);
	t49 = cos(t53);
	t48 = sin(t53);
	t54 = sin(qJ(1));
	t59 = t54 * t48;
	t42 = atan2(-t59, -t49);
	t40 = sin(t42);
	t41 = cos(t42);
	t38 = -t40 * t59 - t41 * t49;
	t37 = 0.1e1 / t38 ^ 2;
	t65 = t37 * t48;
	t64 = t40 * t49;
	t45 = t48 ^ 2;
	t56 = t49 ^ 2;
	t63 = t45 / t56;
	t62 = t48 * t55;
	t57 = t54 ^ 2;
	t61 = 0.1e1 / t57 * t52;
	t43 = 0.1e1 / (t57 * t63 + 0.1e1);
	t60 = t54 * t43;
	t44 = 0.1e1 / (t56 * t61 + 0.1e1);
	t58 = 0.1e1 / t54 * t44 * t62;
	t46 = 0.1e1 / t49;
	t39 = (0.1e1 + t63) * t60;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / (t52 * t45 * t37 + 0.1e1);
	t34 = (t49 * t36 - (-t54 * t64 + t41 * t48 + (-t41 * t59 + t64) * t39) * t65) * t55 * t35;
	t1 = [t46 * t43 * t62, t39, t39, 0, 0, 0; (-t36 * t59 - (-t41 * t45 * t46 * t60 + (t43 - 0.1e1) * t48 * t40) * t52 * t65) * t35, t34, t34, 0, 0, 0; (-0.1e1 - t61) * t49 * t44, -t58, -t58, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (324->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t70 = qJ(2) + qJ(3);
	t68 = sin(t70);
	t69 = cos(t70);
	t72 = sin(qJ(1));
	t80 = t72 * t69;
	t63 = atan2(-t80, t68);
	t61 = sin(t63);
	t62 = cos(t63);
	t54 = -t61 * t80 + t62 * t68;
	t53 = 0.1e1 / t54 ^ 2;
	t74 = cos(qJ(1));
	t86 = t53 * t74 ^ 2;
	t71 = sin(qJ(5));
	t77 = t74 * t71;
	t73 = cos(qJ(5));
	t78 = t72 * t73;
	t60 = t68 * t77 + t78;
	t58 = 0.1e1 / t60 ^ 2;
	t76 = t74 * t73;
	t79 = t72 * t71;
	t59 = -t68 * t76 + t79;
	t85 = t58 * t59;
	t84 = t61 * t68;
	t67 = t69 ^ 2;
	t83 = 0.1e1 / t68 ^ 2 * t67;
	t82 = t69 * t74;
	t64 = 0.1e1 / (t72 ^ 2 * t83 + 0.1e1);
	t81 = t72 * t64;
	t75 = t59 ^ 2 * t58 + 0.1e1;
	t65 = 0.1e1 / t68;
	t57 = 0.1e1 / t60;
	t56 = 0.1e1 / t75;
	t55 = (0.1e1 + t83) * t81;
	t52 = 0.1e1 / t54;
	t51 = 0.1e1 / (t67 * t86 + 0.1e1);
	t50 = (-t57 * t73 - t71 * t85) * t56 * t82;
	t49 = (-t68 * t52 - (t72 * t84 + t62 * t69 + (-t62 * t80 - t84) * t55) * t69 * t53) * t74 * t51;
	t1 = [-t65 * t64 * t82, t55, t55, 0, 0, 0; (-t52 * t80 - (t62 * t65 * t67 * t81 + (t64 - 0.1e1) * t69 * t61) * t69 * t86) * t51, t49, t49, 0, 0, 0; ((t68 * t78 + t77) * t57 - (-t68 * t79 + t76) * t85) * t56, t50, t50, 0, t75 * t56, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (477->26), mult. (688->69), div. (136->11), fcn. (1034->9), ass. (0->43)
	t80 = qJ(2) + qJ(3);
	t76 = cos(t80);
	t97 = t76 ^ 2;
	t75 = sin(t80);
	t81 = sin(qJ(5));
	t84 = cos(qJ(1));
	t87 = t84 * t81;
	t82 = sin(qJ(1));
	t83 = cos(qJ(5));
	t88 = t82 * t83;
	t69 = t75 * t88 + t87;
	t91 = t76 * t83;
	t64 = atan2(t69, t91);
	t60 = sin(t64);
	t61 = cos(t64);
	t59 = t60 * t69 + t61 * t91;
	t58 = 0.1e1 / t59 ^ 2;
	t86 = t84 * t83;
	t89 = t82 * t81;
	t67 = -t75 * t86 + t89;
	t96 = t58 * t67;
	t94 = t61 * t69;
	t93 = t67 ^ 2 * t58;
	t73 = 0.1e1 / t76;
	t77 = 0.1e1 / t83;
	t92 = t73 * t77;
	t90 = t76 * t84;
	t68 = t75 * t87 + t88;
	t66 = 0.1e1 / t68 ^ 2;
	t85 = t84 ^ 2 * t97 * t66;
	t78 = 0.1e1 / t83 ^ 2;
	t74 = 0.1e1 / t97;
	t70 = -t75 * t89 + t86;
	t65 = 0.1e1 / t68;
	t63 = 0.1e1 / (t69 ^ 2 * t74 * t78 + 0.1e1);
	t62 = 0.1e1 / (0.1e1 + t85);
	t57 = 0.1e1 / t59;
	t56 = (t69 * t74 * t75 * t77 + t82) * t63;
	t55 = 0.1e1 / (0.1e1 + t93);
	t54 = (t69 * t78 * t81 + t70 * t77) * t73 * t63;
	t53 = (t65 * t75 * t84 + t81 * t85) * t62;
	t52 = (-t56 * t94 * t96 + (-t57 * t90 - (-t61 * t75 + (-t56 + t82) * t60 * t76) * t96) * t83) * t55;
	t1 = [-t67 * t63 * t92, t56, t56, 0, t54, 0; (t69 * t57 - (-t60 + (-t92 * t94 + t60) * t63) * t93) * t55, t52, t52, 0, (t68 * t57 - (-t61 * t76 * t81 + t60 * t70 + (-t60 * t91 + t94) * t54) * t96) * t55, 0; (t66 * t70 * t84 + t65 * t82) * t76 * t62, t53, t53, 0, -t67 * t66 * t62 * t90, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end