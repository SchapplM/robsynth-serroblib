% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRPPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPPP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:41
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:41
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:40
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:41
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (151->26), mult. (410->65), div. (61->11), fcn. (606->11), ass. (0->42)
	t48 = sin(qJ(2));
	t47 = cos(pkin(5));
	t51 = cos(qJ(1));
	t54 = t51 * t47;
	t45 = sin(pkin(5));
	t49 = sin(qJ(1));
	t61 = t45 * t49;
	t37 = t48 * t61 - t54;
	t50 = cos(qJ(2));
	t56 = t50 * t45;
	t35 = atan2(-t37, -t56);
	t33 = sin(t35);
	t34 = cos(t35);
	t27 = -t33 * t37 - t34 * t56;
	t25 = 0.1e1 / t27 ^ 2;
	t58 = t49 * t47;
	t60 = t45 * t51;
	t39 = t48 * t60 + t58;
	t66 = t25 * t39;
	t44 = sin(pkin(8));
	t46 = cos(pkin(8));
	t53 = t48 * t54 - t61;
	t55 = t50 * t51;
	t32 = -t53 * t44 + t46 * t55;
	t30 = 0.1e1 / t32 ^ 2;
	t31 = t44 * t55 + t53 * t46;
	t65 = t30 * t31;
	t64 = t34 * t37;
	t63 = t39 ^ 2 * t25;
	t41 = 0.1e1 / t45;
	t62 = t41 / t50;
	t59 = t47 * t50;
	t57 = t49 * t50;
	t52 = t48 * t58 + t60;
	t43 = 0.1e1 / t50 ^ 2;
	t36 = 0.1e1 / (0.1e1 + t37 ^ 2 * t43 / t45 ^ 2);
	t29 = 0.1e1 / t32;
	t28 = 0.1e1 / (t31 ^ 2 * t30 + 0.1e1);
	t26 = (t37 * t41 * t43 * t48 + t49) * t36;
	t24 = 0.1e1 / t27;
	t23 = 0.1e1 / (0.1e1 + t63);
	t1 = [t39 * t36 * t62, t26, 0, 0, 0; (-t37 * t24 - (-t33 + (-t62 * t64 + t33) * t36) * t63) * t23, (t26 * t64 * t66 + (t24 * t55 - (t34 * t48 + (t26 * t50 - t57) * t33) * t66) * t45) * t23, 0, 0, 0; ((-t44 * t57 - t52 * t46) * t29 - (t52 * t44 - t46 * t57) * t65) * t28, ((-t44 * t48 + t46 * t59) * t29 - (-t44 * t59 - t46 * t48) * t65) * t28 * t51, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:41
	% EndTime: 2019-12-29 18:08:42
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (231->27), mult. (733->64), div. (50->9), fcn. (1031->11), ass. (0->41)
	t68 = cos(pkin(8));
	t71 = sin(qJ(1));
	t70 = sin(qJ(2));
	t69 = cos(pkin(5));
	t78 = t71 * t69;
	t67 = sin(pkin(5));
	t73 = cos(qJ(1));
	t81 = t67 * t73;
	t74 = t70 * t78 + t81;
	t66 = sin(pkin(8));
	t72 = cos(qJ(2));
	t77 = t72 * t66;
	t53 = t74 * t68 + t71 * t77;
	t79 = t68 * t72;
	t75 = -t70 * t66 + t69 * t79;
	t52 = atan2(-t53, -t75);
	t49 = cos(t52);
	t85 = t49 * t53;
	t48 = sin(t52);
	t47 = -t48 * t53 - t49 * t75;
	t46 = 0.1e1 / t47 ^ 2;
	t80 = t68 * t70;
	t63 = t69 * t80 + t77;
	t82 = t67 * t71;
	t55 = t63 * t73 - t68 * t82;
	t84 = t55 ^ 2 * t46;
	t76 = t73 * t69;
	t56 = t73 * t79 + (-t70 * t76 + t82) * t66;
	t64 = t70 * t81 + t78;
	t61 = 0.1e1 / t64 ^ 2;
	t83 = t56 * t61;
	t60 = 0.1e1 / t64;
	t59 = 0.1e1 / t75 ^ 2;
	t58 = 0.1e1 / t75;
	t57 = t75 * t71;
	t51 = 0.1e1 / (t56 ^ 2 * t61 + 0.1e1);
	t50 = 0.1e1 / (t53 ^ 2 * t59 + 0.1e1);
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / (0.1e1 + t84);
	t43 = (t53 * t59 * t63 + t57 * t58) * t50;
	t1 = [t55 * t58 * t50, t43, 0, 0, 0; (-t53 * t45 - (-t48 + (-t58 * t85 + t48) * t50) * t84) * t44, (-(-t48 * t57 + t49 * t63 + (t48 * t75 - t85) * t43) * t55 * t46 + t75 * t45 * t73) * t44, 0, 0, 0; ((t74 * t66 - t71 * t79) * t60 - (-t70 * t82 + t76) * t83) * t51, ((-t69 * t77 - t80) * t60 - t72 * t67 * t83) * t51 * t73, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:42
	% EndTime: 2019-12-29 18:08:42
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (231->28), mult. (733->66), div. (50->9), fcn. (1031->11), ass. (0->42)
	t77 = cos(pkin(5));
	t74 = sin(pkin(8));
	t78 = sin(qJ(2));
	t90 = t74 * t78;
	t82 = t77 * t90;
	t76 = cos(pkin(8));
	t80 = cos(qJ(2));
	t85 = t80 * t76;
	t69 = -t82 + t85;
	t79 = sin(qJ(1));
	t75 = sin(pkin(5));
	t81 = cos(qJ(1));
	t87 = t75 * t81;
	t59 = t69 * t79 - t74 * t87;
	t89 = t74 * t80;
	t68 = t78 * t76 + t77 * t89;
	t58 = atan2(-t59, t68);
	t55 = cos(t58);
	t93 = t55 * t59;
	t83 = t81 * t77;
	t84 = t80 * t81;
	t88 = t75 * t79;
	t61 = -t74 * t84 + (-t78 * t83 + t88) * t76;
	t86 = t79 * t77;
	t70 = t78 * t87 + t86;
	t67 = 0.1e1 / t70 ^ 2;
	t92 = t61 * t67;
	t54 = sin(t58);
	t53 = -t54 * t59 + t55 * t68;
	t52 = 0.1e1 / t53 ^ 2;
	t62 = t74 * t88 + t76 * t84 - t81 * t82;
	t91 = t62 ^ 2 * t52;
	t66 = 0.1e1 / t70;
	t65 = 0.1e1 / t68 ^ 2;
	t64 = 0.1e1 / t68;
	t63 = t68 * t79;
	t57 = 0.1e1 / (t61 ^ 2 * t67 + 0.1e1);
	t56 = 0.1e1 / (t59 ^ 2 * t65 + 0.1e1);
	t51 = 0.1e1 / t53;
	t50 = 0.1e1 / (0.1e1 + t91);
	t49 = (t59 * t65 * t69 + t63 * t64) * t56;
	t1 = [-t62 * t64 * t56, t49, 0, 0, 0; (-t59 * t51 - (-t54 + (t64 * t93 + t54) * t56) * t91) * t50, (-(t54 * t63 + t55 * t69 + (-t54 * t68 - t93) * t49) * t62 * t52 - t68 * t51 * t81) * t50, 0, 0, 0; ((t79 * t89 + (t78 * t86 + t87) * t76) * t66 - (-t78 * t88 + t83) * t92) * t57, ((-t77 * t85 + t90) * t66 - t80 * t75 * t92) * t57 * t81, 0, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end