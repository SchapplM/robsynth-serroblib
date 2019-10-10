% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR4
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
%   Wie in S6PRPRPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (84->18), mult. (220->41), div. (42->11), fcn. (332->11), ass. (0->27)
	t39 = sin(pkin(10));
	t40 = sin(pkin(6));
	t49 = t39 * t40;
	t45 = cos(qJ(2));
	t48 = t40 * t45;
	t43 = cos(pkin(6));
	t44 = sin(qJ(2));
	t47 = t43 * t44;
	t46 = t43 * t45;
	t42 = cos(pkin(10));
	t41 = cos(pkin(11));
	t38 = sin(pkin(11));
	t37 = 0.1e1 / t45 ^ 2;
	t34 = -t39 * t47 + t42 * t45;
	t33 = t39 * t46 + t42 * t44;
	t32 = t39 * t45 + t42 * t47;
	t30 = t39 * t44 - t42 * t46;
	t28 = atan2(-t30, -t48);
	t27 = cos(t28);
	t26 = sin(t28);
	t25 = t34 * t41 + t38 * t49;
	t24 = t34 * t38 - t41 * t49;
	t23 = 0.1e1 / t25 ^ 2;
	t21 = -t26 * t30 - t27 * t48;
	t20 = 0.1e1 / t21 ^ 2;
	t18 = (t32 / t45 + t44 * t30 * t37) / t40 / (0.1e1 + t30 ^ 2 / t40 ^ 2 * t37);
	t1 = [0, t18, 0, 0, 0, 0; 0, (t34 / t21 - (t27 * t40 * t44 - t26 * t32 + (t26 * t48 - t27 * t30) * t18) * t33 * t20) / (t33 ^ 2 * t20 + 0.1e1), 0, 0, 0, 0; 0, (-t38 / t25 + t41 * t24 * t23) * t33 / (t24 ^ 2 * t23 + 0.1e1), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (131->19), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->30)
	t47 = sin(pkin(10));
	t48 = sin(pkin(6));
	t57 = t47 * t48;
	t52 = cos(qJ(2));
	t56 = t48 * t52;
	t50 = cos(pkin(6));
	t51 = sin(qJ(2));
	t55 = t50 * t51;
	t54 = t50 * t52;
	t49 = cos(pkin(10));
	t40 = -t47 * t55 + t49 * t52;
	t45 = pkin(11) + qJ(4);
	t42 = sin(t45);
	t43 = cos(t45);
	t31 = t40 * t43 + t42 * t57;
	t29 = 0.1e1 / t31 ^ 2;
	t30 = t40 * t42 - t43 * t57;
	t53 = t30 ^ 2 * t29 + 0.1e1;
	t46 = 0.1e1 / t52 ^ 2;
	t39 = t47 * t54 + t49 * t51;
	t38 = t47 * t52 + t49 * t55;
	t36 = t47 * t51 - t49 * t54;
	t34 = atan2(-t36, -t56);
	t33 = cos(t34);
	t32 = sin(t34);
	t28 = 0.1e1 / t53;
	t27 = -t32 * t36 - t33 * t56;
	t26 = 0.1e1 / t27 ^ 2;
	t24 = (t38 / t52 + t51 * t36 * t46) / t48 / (0.1e1 + t36 ^ 2 / t48 ^ 2 * t46);
	t1 = [0, t24, 0, 0, 0, 0; 0, (t40 / t27 - (t33 * t48 * t51 - t32 * t38 + (t32 * t56 - t33 * t36) * t24) * t39 * t26) / (t39 ^ 2 * t26 + 0.1e1), 0, 0, 0, 0; 0, (-t42 / t31 + t43 * t30 * t29) * t39 * t28, 0, t53 * t28, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (596->31), mult. (866->78), div. (60->9), fcn. (1237->13), ass. (0->51)
	t70 = sin(pkin(10));
	t73 = cos(pkin(10));
	t76 = cos(qJ(2));
	t74 = cos(pkin(6));
	t75 = sin(qJ(2));
	t79 = t74 * t75;
	t62 = t70 * t76 + t73 * t79;
	t68 = pkin(11) + qJ(4);
	t66 = sin(t68);
	t67 = cos(t68);
	t71 = sin(pkin(6));
	t82 = t71 * t73;
	t52 = t62 * t66 + t67 * t82;
	t81 = t71 * t75;
	t59 = t66 * t81 - t74 * t67;
	t51 = atan2(-t52, t59);
	t48 = sin(t51);
	t49 = cos(t51);
	t42 = -t48 * t52 + t49 * t59;
	t41 = 0.1e1 / t42 ^ 2;
	t64 = -t70 * t79 + t73 * t76;
	t83 = t70 * t71;
	t55 = t64 * t66 - t67 * t83;
	t88 = t41 * t55;
	t56 = t64 * t67 + t66 * t83;
	t72 = cos(pkin(12));
	t78 = t74 * t76;
	t63 = t70 * t78 + t73 * t75;
	t69 = sin(pkin(12));
	t85 = t63 * t69;
	t47 = t56 * t72 + t85;
	t45 = 0.1e1 / t47 ^ 2;
	t84 = t63 * t72;
	t46 = t56 * t69 - t84;
	t87 = t45 * t46;
	t58 = 0.1e1 / t59 ^ 2;
	t86 = t52 * t58;
	t80 = t71 * t76;
	t77 = -t48 * t59 - t49 * t52;
	t61 = -t70 * t75 + t73 * t78;
	t60 = t74 * t66 + t67 * t81;
	t57 = 0.1e1 / t59;
	t54 = t62 * t67 - t66 * t82;
	t50 = 0.1e1 / (t52 ^ 2 * t58 + 0.1e1);
	t44 = 0.1e1 / t47;
	t43 = 0.1e1 / (t46 ^ 2 * t45 + 0.1e1);
	t40 = 0.1e1 / t42;
	t39 = 0.1e1 / (t55 ^ 2 * t41 + 0.1e1);
	t38 = (-t57 * t61 + t80 * t86) * t66 * t50;
	t37 = (-t54 * t57 + t60 * t86) * t50;
	t1 = [0, t38, 0, t37, 0, 0; 0, (-t63 * t66 * t40 - ((-t48 * t61 + t49 * t80) * t66 + t77 * t38) * t88) * t39, 0, (t56 * t40 - (t77 * t37 - t48 * t54 + t49 * t60) * t88) * t39, 0, 0; 0, ((-t64 * t72 - t67 * t85) * t44 - (t64 * t69 - t67 * t84) * t87) * t43, 0, (-t44 * t69 + t72 * t87) * t55 * t43, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (689->32), mult. (949->80), div. (65->9), fcn. (1349->13), ass. (0->52)
	t77 = sin(pkin(10));
	t79 = cos(pkin(10));
	t82 = cos(qJ(2));
	t80 = cos(pkin(6));
	t81 = sin(qJ(2));
	t86 = t80 * t81;
	t67 = t77 * t82 + t79 * t86;
	t76 = pkin(11) + qJ(4);
	t72 = sin(t76);
	t74 = cos(t76);
	t78 = sin(pkin(6));
	t89 = t78 * t79;
	t57 = t67 * t72 + t74 * t89;
	t88 = t78 * t81;
	t64 = t72 * t88 - t80 * t74;
	t56 = atan2(-t57, t64);
	t53 = sin(t56);
	t54 = cos(t56);
	t47 = -t53 * t57 + t54 * t64;
	t46 = 0.1e1 / t47 ^ 2;
	t69 = -t77 * t86 + t79 * t82;
	t90 = t77 * t78;
	t60 = t69 * t72 - t74 * t90;
	t94 = t46 * t60;
	t61 = t69 * t74 + t72 * t90;
	t85 = t80 * t82;
	t68 = t77 * t85 + t79 * t81;
	t75 = pkin(12) + qJ(6);
	t71 = sin(t75);
	t73 = cos(t75);
	t52 = t61 * t73 + t68 * t71;
	t50 = 0.1e1 / t52 ^ 2;
	t51 = t61 * t71 - t68 * t73;
	t93 = t50 * t51;
	t63 = 0.1e1 / t64 ^ 2;
	t92 = t57 * t63;
	t91 = t68 * t74;
	t87 = t78 * t82;
	t84 = t51 ^ 2 * t50 + 0.1e1;
	t83 = -t53 * t64 - t54 * t57;
	t66 = -t77 * t81 + t79 * t85;
	t65 = t80 * t72 + t74 * t88;
	t62 = 0.1e1 / t64;
	t59 = t67 * t74 - t72 * t89;
	t55 = 0.1e1 / (t57 ^ 2 * t63 + 0.1e1);
	t49 = 0.1e1 / t52;
	t48 = 0.1e1 / t84;
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / (t60 ^ 2 * t46 + 0.1e1);
	t43 = (-t62 * t66 + t87 * t92) * t72 * t55;
	t42 = (-t59 * t62 + t65 * t92) * t55;
	t1 = [0, t43, 0, t42, 0, 0; 0, (-t68 * t72 * t45 - ((-t53 * t66 + t54 * t87) * t72 + t83 * t43) * t94) * t44, 0, (t61 * t45 - (t83 * t42 - t53 * t59 + t54 * t65) * t94) * t44, 0, 0; 0, ((-t69 * t73 - t71 * t91) * t49 - (t69 * t71 - t73 * t91) * t93) * t48, 0, (-t49 * t71 + t73 * t93) * t60 * t48, 0, t84 * t48;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end