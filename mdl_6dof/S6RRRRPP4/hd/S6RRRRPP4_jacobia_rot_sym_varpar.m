% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPP4
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
%   Wie in S6RRRRPP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:25
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPP4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:46
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (109->20), mult. (224->55), div. (52->9), fcn. (332->9), ass. (0->33)
	t40 = cos(qJ(2));
	t37 = sin(qJ(2));
	t38 = sin(qJ(1));
	t46 = t38 * t37;
	t31 = atan2(-t46, -t40);
	t29 = sin(t31);
	t30 = cos(t31);
	t22 = -t29 * t46 - t30 * t40;
	t21 = 0.1e1 / t22 ^ 2;
	t41 = cos(qJ(1));
	t51 = t21 * t41 ^ 2;
	t36 = sin(qJ(3));
	t39 = cos(qJ(3));
	t43 = t41 * t39;
	t28 = t38 * t36 + t40 * t43;
	t26 = 0.1e1 / t28 ^ 2;
	t44 = t41 * t36;
	t27 = -t38 * t39 + t40 * t44;
	t50 = t26 * t27;
	t33 = t37 ^ 2;
	t49 = t33 / t40 ^ 2;
	t48 = t37 * t41;
	t32 = 0.1e1 / (t38 ^ 2 * t49 + 0.1e1);
	t47 = t38 * t32;
	t45 = t38 * t40;
	t42 = t27 ^ 2 * t26 + 0.1e1;
	t34 = 0.1e1 / t40;
	t25 = 0.1e1 / t28;
	t24 = (0.1e1 + t49) * t47;
	t23 = 0.1e1 / t42;
	t20 = 0.1e1 / t22;
	t19 = 0.1e1 / (t33 * t51 + 0.1e1);
	t1 = [t34 * t32 * t48, t24, 0, 0, 0, 0; (-t20 * t46 - (-t30 * t33 * t34 * t47 + (t32 - 0.1e1) * t37 * t29) * t37 * t51) * t19, (t40 * t20 - (-t29 * t45 + t30 * t37 + (t29 * t40 - t30 * t46) * t24) * t37 * t21) * t41 * t19, 0, 0, 0, 0; ((-t36 * t45 - t43) * t25 - (-t39 * t45 + t44) * t50) * t23, (-t25 * t36 + t39 * t50) * t23 * t48, t42 * t23, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:46
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (181->21), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->35)
	t56 = cos(qJ(2));
	t54 = sin(qJ(2));
	t55 = sin(qJ(1));
	t62 = t55 * t54;
	t46 = atan2(-t62, -t56);
	t44 = sin(t46);
	t45 = cos(t46);
	t38 = -t44 * t62 - t45 * t56;
	t37 = 0.1e1 / t38 ^ 2;
	t57 = cos(qJ(1));
	t67 = t37 * t57 ^ 2;
	t53 = qJ(3) + qJ(4);
	t48 = sin(t53);
	t49 = cos(t53);
	t59 = t57 * t49;
	t43 = t55 * t48 + t56 * t59;
	t41 = 0.1e1 / t43 ^ 2;
	t60 = t57 * t48;
	t42 = -t55 * t49 + t56 * t60;
	t66 = t41 * t42;
	t50 = t54 ^ 2;
	t65 = t50 / t56 ^ 2;
	t64 = t54 * t57;
	t47 = 0.1e1 / (t55 ^ 2 * t65 + 0.1e1);
	t63 = t55 * t47;
	t61 = t55 * t56;
	t58 = t42 ^ 2 * t41 + 0.1e1;
	t51 = 0.1e1 / t56;
	t40 = 0.1e1 / t43;
	t39 = (0.1e1 + t65) * t63;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / t58;
	t34 = 0.1e1 / (t50 * t67 + 0.1e1);
	t33 = t58 * t35;
	t1 = [t51 * t47 * t64, t39, 0, 0, 0, 0; (-t36 * t62 - (-t45 * t50 * t51 * t63 + (t47 - 0.1e1) * t54 * t44) * t54 * t67) * t34, (t56 * t36 - (-t44 * t61 + t45 * t54 + (t44 * t56 - t45 * t62) * t39) * t54 * t37) * t57 * t34, 0, 0, 0, 0; ((-t48 * t61 - t59) * t40 - (-t49 * t61 + t60) * t66) * t35, (-t40 * t48 + t49 * t66) * t35 * t64, t33, t33, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:46
	% EndTime: 2019-10-10 12:25:46
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (243->21), mult. (251->55), div. (57->9), fcn. (367->9), ass. (0->35)
	t63 = cos(qJ(2));
	t61 = sin(qJ(2));
	t62 = sin(qJ(1));
	t69 = t62 * t61;
	t53 = atan2(-t69, -t63);
	t51 = sin(t53);
	t52 = cos(t53);
	t45 = -t51 * t69 - t52 * t63;
	t44 = 0.1e1 / t45 ^ 2;
	t64 = cos(qJ(1));
	t74 = t44 * t64 ^ 2;
	t57 = qJ(3) + qJ(4) + pkin(10);
	t55 = sin(t57);
	t56 = cos(t57);
	t66 = t64 * t56;
	t50 = t62 * t55 + t63 * t66;
	t48 = 0.1e1 / t50 ^ 2;
	t67 = t64 * t55;
	t49 = -t62 * t56 + t63 * t67;
	t73 = t48 * t49;
	t58 = t61 ^ 2;
	t72 = t58 / t63 ^ 2;
	t71 = t61 * t64;
	t54 = 0.1e1 / (t62 ^ 2 * t72 + 0.1e1);
	t70 = t62 * t54;
	t68 = t62 * t63;
	t65 = t49 ^ 2 * t48 + 0.1e1;
	t59 = 0.1e1 / t63;
	t47 = 0.1e1 / t50;
	t46 = (0.1e1 + t72) * t70;
	t43 = 0.1e1 / t45;
	t42 = 0.1e1 / (t58 * t74 + 0.1e1);
	t41 = 0.1e1 / t65;
	t40 = t65 * t41;
	t1 = [t59 * t54 * t71, t46, 0, 0, 0, 0; (-t43 * t69 - (-t52 * t58 * t59 * t70 + (t54 - 0.1e1) * t61 * t51) * t61 * t74) * t42, (t63 * t43 - (-t51 * t68 + t52 * t61 + (t51 * t63 - t52 * t69) * t46) * t61 * t44) * t64 * t42, 0, 0, 0, 0; ((-t55 * t68 - t66) * t47 - (-t56 * t68 + t67) * t73) * t41, (-t47 * t55 + t56 * t73) * t41 * t71, t40, t40, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:25:45
	% EndTime: 2019-10-10 12:25:46
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (1004->28), mult. (689->68), div. (139->11), fcn. (1046->9), ass. (0->42)
	t81 = sin(qJ(2));
	t96 = t81 ^ 2;
	t76 = qJ(3) + qJ(4) + pkin(10);
	t74 = sin(t76);
	t75 = cos(t76);
	t84 = cos(qJ(1));
	t86 = t84 * t75;
	t82 = sin(qJ(1));
	t83 = cos(qJ(2));
	t88 = t82 * t83;
	t65 = t74 * t88 + t86;
	t90 = t81 * t74;
	t61 = atan2(-t65, t90);
	t58 = sin(t61);
	t59 = cos(t61);
	t56 = -t58 * t65 + t59 * t90;
	t55 = 0.1e1 / t56 ^ 2;
	t87 = t84 * t74;
	t68 = -t82 * t75 + t83 * t87;
	t95 = t55 * t68;
	t93 = t59 * t65;
	t92 = t68 ^ 2 * t55;
	t72 = 0.1e1 / t74;
	t78 = 0.1e1 / t81;
	t91 = t72 * t78;
	t89 = t81 * t84;
	t69 = t82 * t74 + t83 * t86;
	t64 = 0.1e1 / t69 ^ 2;
	t85 = t84 ^ 2 * t96 * t64;
	t79 = 0.1e1 / t96;
	t73 = 0.1e1 / t74 ^ 2;
	t67 = t75 * t88 - t87;
	t63 = 0.1e1 / t69;
	t62 = 0.1e1 / (0.1e1 + t85);
	t60 = 0.1e1 / (t65 ^ 2 * t79 * t73 + 0.1e1);
	t57 = t68 * t64 * t62 * t89;
	t54 = 0.1e1 / t56;
	t53 = (t65 * t72 * t79 * t83 + t82) * t60;
	t52 = 0.1e1 / (0.1e1 + t92);
	t51 = (t65 * t73 * t75 - t67 * t72) * t78 * t60;
	t50 = (t69 * t54 - (t59 * t81 * t75 - t58 * t67 + (-t58 * t90 - t93) * t51) * t95) * t52;
	t1 = [-t68 * t60 * t91, t53, t51, t51, 0, 0; (-t65 * t54 - (-t58 + (t91 * t93 + t58) * t60) * t92) * t52, (t53 * t93 * t95 + (-t54 * t89 - (t59 * t83 + (-t53 + t82) * t81 * t58) * t95) * t74) * t52, t50, t50, 0, 0; (-t64 * t67 * t84 + t63 * t82) * t81 * t62, (-t63 * t83 * t84 - t75 * t85) * t62, -t57, -t57, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end