% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP3
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
%   Wie in S6RRRPRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:41
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (341->21), mult. (305->53), div. (74->9), fcn. (454->9), ass. (0->37)
	t62 = qJ(2) + qJ(3);
	t61 = cos(t62);
	t60 = sin(t62);
	t65 = sin(qJ(1));
	t71 = t65 * t60;
	t55 = atan2(-t71, -t61);
	t53 = sin(t55);
	t54 = cos(t55);
	t46 = -t53 * t71 - t54 * t61;
	t45 = 0.1e1 / t46 ^ 2;
	t66 = cos(qJ(1));
	t77 = t45 * t66 ^ 2;
	t64 = cos(pkin(10));
	t67 = t66 * t64;
	t63 = sin(pkin(10));
	t70 = t65 * t63;
	t52 = t61 * t67 + t70;
	t50 = 0.1e1 / t52 ^ 2;
	t68 = t66 * t63;
	t69 = t65 * t64;
	t51 = t61 * t68 - t69;
	t76 = t50 * t51;
	t75 = t53 * t61;
	t57 = t60 ^ 2;
	t74 = t57 / t61 ^ 2;
	t73 = t60 * t66;
	t56 = 0.1e1 / (t65 ^ 2 * t74 + 0.1e1);
	t72 = t65 * t56;
	t58 = 0.1e1 / t61;
	t49 = 0.1e1 / t52;
	t48 = 0.1e1 / (t51 ^ 2 * t50 + 0.1e1);
	t47 = (0.1e1 + t74) * t72;
	t44 = 0.1e1 / t46;
	t43 = 0.1e1 / (t57 * t77 + 0.1e1);
	t42 = (-t49 * t63 + t64 * t76) * t48 * t73;
	t41 = (t61 * t44 - (-t65 * t75 + t54 * t60 + (-t54 * t71 + t75) * t47) * t60 * t45) * t66 * t43;
	t1 = [t58 * t56 * t73, t47, t47, 0, 0, 0; (-t44 * t71 - (-t54 * t57 * t58 * t72 + (t56 - 0.1e1) * t60 * t53) * t60 * t77) * t43, t41, t41, 0, 0, 0; ((-t61 * t70 - t67) * t49 - (-t61 * t69 + t68) * t76) * t48, t42, t42, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:41
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (422->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
	t66 = qJ(2) + qJ(3);
	t64 = cos(t66);
	t63 = sin(t66);
	t67 = sin(qJ(1));
	t72 = t67 * t63;
	t56 = atan2(-t72, -t64);
	t54 = sin(t56);
	t55 = cos(t56);
	t47 = -t54 * t72 - t55 * t64;
	t46 = 0.1e1 / t47 ^ 2;
	t68 = cos(qJ(1));
	t80 = t46 * t68 ^ 2;
	t65 = pkin(10) + qJ(5);
	t62 = cos(t65);
	t70 = t68 * t62;
	t61 = sin(t65);
	t74 = t67 * t61;
	t53 = t64 * t70 + t74;
	t51 = 0.1e1 / t53 ^ 2;
	t71 = t68 * t61;
	t73 = t67 * t62;
	t52 = t64 * t71 - t73;
	t79 = t51 * t52;
	t78 = t54 * t64;
	t58 = t63 ^ 2;
	t76 = t58 / t64 ^ 2;
	t57 = 0.1e1 / (t67 ^ 2 * t76 + 0.1e1);
	t77 = t57 * t67;
	t75 = t63 * t68;
	t69 = t52 ^ 2 * t51 + 0.1e1;
	t59 = 0.1e1 / t64;
	t50 = 0.1e1 / t53;
	t49 = (0.1e1 + t76) * t77;
	t48 = 0.1e1 / t69;
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / (t58 * t80 + 0.1e1);
	t43 = (-t50 * t61 + t62 * t79) * t48 * t75;
	t42 = (t64 * t45 - (-t67 * t78 + t55 * t63 + (-t55 * t72 + t78) * t49) * t63 * t46) * t68 * t44;
	t1 = [t59 * t57 * t75, t49, t49, 0, 0, 0; (-t45 * t72 - (-t55 * t58 * t59 * t77 + (t57 - 0.1e1) * t63 * t54) * t63 * t80) * t44, t42, t42, 0, 0, 0; ((-t64 * t74 - t70) * t50 - (-t64 * t73 + t71) * t79) * t48, t43, t43, 0, t69 * t48, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:41
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (860->28), mult. (688->69), div. (136->11), fcn. (1034->9), ass. (0->44)
	t82 = qJ(2) + qJ(3);
	t78 = sin(t82);
	t97 = t78 ^ 2;
	t80 = pkin(10) + qJ(5);
	t77 = cos(t80);
	t79 = cos(t82);
	t84 = cos(qJ(1));
	t76 = sin(t80);
	t83 = sin(qJ(1));
	t88 = t83 * t76;
	t64 = t77 * t84 + t79 * t88;
	t91 = t78 * t76;
	t60 = atan2(-t64, t91);
	t57 = sin(t60);
	t58 = cos(t60);
	t56 = -t57 * t64 + t58 * t91;
	t55 = 0.1e1 / t56 ^ 2;
	t86 = t84 * t76;
	t87 = t83 * t77;
	t67 = t79 * t86 - t87;
	t96 = t55 * t67;
	t95 = t55 * t67 ^ 2;
	t93 = t58 * t64;
	t71 = 0.1e1 / t76;
	t74 = 0.1e1 / t78;
	t92 = t71 * t74;
	t90 = t78 * t84;
	t89 = t79 * t84;
	t68 = t77 * t89 + t88;
	t63 = 0.1e1 / t68 ^ 2;
	t85 = t63 * t97 * t84 ^ 2;
	t75 = 0.1e1 / t97;
	t72 = 0.1e1 / t76 ^ 2;
	t66 = t79 * t87 - t86;
	t62 = 0.1e1 / t68;
	t61 = 0.1e1 / (0.1e1 + t85);
	t59 = 0.1e1 / (t64 ^ 2 * t72 * t75 + 0.1e1);
	t54 = 0.1e1 / t56;
	t53 = (t64 * t71 * t75 * t79 + t83) * t59;
	t52 = (-t62 * t89 - t77 * t85) * t61;
	t51 = 0.1e1 / (0.1e1 + t95);
	t50 = (t64 * t72 * t77 - t66 * t71) * t74 * t59;
	t49 = (t53 * t93 * t96 + (-t54 * t90 - (t58 * t79 + (-t53 + t83) * t78 * t57) * t96) * t76) * t51;
	t1 = [-t67 * t59 * t92, t53, t53, 0, t50, 0; (-t64 * t54 - (-t57 + (t92 * t93 + t57) * t59) * t95) * t51, t49, t49, 0, (t68 * t54 - (t58 * t77 * t78 - t57 * t66 + (-t57 * t91 - t93) * t50) * t96) * t51, 0; (-t63 * t66 * t84 + t62 * t83) * t78 * t61, t52, t52, 0, -t67 * t63 * t61 * t90, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end