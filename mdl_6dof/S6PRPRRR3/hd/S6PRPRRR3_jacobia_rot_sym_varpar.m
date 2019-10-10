% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR3
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
%   Wie in S6PRPRRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:57
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (84->18), mult. (220->41), div. (42->11), fcn. (332->11), ass. (0->27)
	t39 = sin(pkin(11));
	t40 = sin(pkin(6));
	t49 = t39 * t40;
	t45 = cos(qJ(2));
	t48 = t40 * t45;
	t43 = cos(pkin(6));
	t44 = sin(qJ(2));
	t47 = t43 * t44;
	t46 = t43 * t45;
	t42 = cos(pkin(11));
	t41 = cos(pkin(12));
	t38 = sin(pkin(12));
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
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (131->19), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->30)
	t47 = sin(pkin(11));
	t48 = sin(pkin(6));
	t57 = t47 * t48;
	t52 = cos(qJ(2));
	t56 = t48 * t52;
	t50 = cos(pkin(6));
	t51 = sin(qJ(2));
	t55 = t50 * t51;
	t54 = t50 * t52;
	t49 = cos(pkin(11));
	t40 = -t47 * t55 + t49 * t52;
	t45 = pkin(12) + qJ(4);
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
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (206->19), mult. (316->43), div. (52->11), fcn. (458->11), ass. (0->31)
	t60 = sin(pkin(11));
	t61 = sin(pkin(6));
	t70 = t60 * t61;
	t65 = cos(qJ(2));
	t69 = t61 * t65;
	t63 = cos(pkin(6));
	t64 = sin(qJ(2));
	t68 = t63 * t64;
	t67 = t63 * t65;
	t62 = cos(pkin(11));
	t53 = -t60 * t68 + t62 * t65;
	t57 = pkin(12) + qJ(4) + qJ(5);
	t55 = sin(t57);
	t56 = cos(t57);
	t44 = t53 * t56 + t55 * t70;
	t42 = 0.1e1 / t44 ^ 2;
	t43 = t53 * t55 - t56 * t70;
	t66 = t43 ^ 2 * t42 + 0.1e1;
	t59 = 0.1e1 / t65 ^ 2;
	t52 = t60 * t67 + t62 * t64;
	t51 = t60 * t65 + t62 * t68;
	t49 = t60 * t64 - t62 * t67;
	t47 = atan2(-t49, -t69);
	t46 = cos(t47);
	t45 = sin(t47);
	t41 = 0.1e1 / t66;
	t40 = -t45 * t49 - t46 * t69;
	t39 = 0.1e1 / t40 ^ 2;
	t37 = (t51 / t65 + t64 * t49 * t59) / t61 / (0.1e1 + t49 ^ 2 / t61 ^ 2 * t59);
	t36 = t66 * t41;
	t1 = [0, t37, 0, 0, 0, 0; 0, (t53 / t40 - (t46 * t61 * t64 - t45 * t51 + (t45 * t69 - t46 * t49) * t37) * t52 * t39) / (t52 ^ 2 * t39 + 0.1e1), 0, 0, 0, 0; 0, (-t55 / t44 + t56 * t43 * t42) * t52 * t41, 0, t36, t36, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (1407->31), mult. (1387->79), div. (95->9), fcn. (1976->13), ass. (0->54)
	t93 = sin(pkin(6));
	t94 = cos(pkin(11));
	t106 = t93 * t94;
	t95 = cos(pkin(6));
	t97 = sin(qJ(2));
	t103 = t95 * t97;
	t92 = sin(pkin(11));
	t99 = cos(qJ(2));
	t85 = t94 * t103 + t92 * t99;
	t91 = pkin(12) + qJ(4) + qJ(5);
	t89 = sin(t91);
	t90 = cos(t91);
	t75 = t90 * t106 + t85 * t89;
	t105 = t93 * t97;
	t82 = t89 * t105 - t95 * t90;
	t70 = atan2(-t75, t82);
	t67 = sin(t70);
	t68 = cos(t70);
	t65 = -t67 * t75 + t68 * t82;
	t64 = 0.1e1 / t65 ^ 2;
	t107 = t92 * t93;
	t87 = -t92 * t103 + t94 * t99;
	t78 = -t90 * t107 + t87 * t89;
	t112 = t64 * t78;
	t102 = t95 * t99;
	t86 = t92 * t102 + t94 * t97;
	t96 = sin(qJ(6));
	t109 = t86 * t96;
	t79 = t89 * t107 + t87 * t90;
	t98 = cos(qJ(6));
	t74 = t79 * t98 + t109;
	t72 = 0.1e1 / t74 ^ 2;
	t108 = t86 * t98;
	t73 = t79 * t96 - t108;
	t111 = t72 * t73;
	t81 = 0.1e1 / t82 ^ 2;
	t110 = t75 * t81;
	t104 = t93 * t99;
	t101 = t73 ^ 2 * t72 + 0.1e1;
	t100 = -t67 * t82 - t68 * t75;
	t84 = t94 * t102 - t92 * t97;
	t83 = t90 * t105 + t95 * t89;
	t80 = 0.1e1 / t82;
	t77 = -t89 * t106 + t85 * t90;
	t71 = 0.1e1 / t74;
	t69 = 0.1e1 / (t75 ^ 2 * t81 + 0.1e1);
	t66 = 0.1e1 / t101;
	t63 = 0.1e1 / t65;
	t62 = 0.1e1 / (t78 ^ 2 * t64 + 0.1e1);
	t61 = (t104 * t110 - t80 * t84) * t89 * t69;
	t60 = (t83 * t110 - t77 * t80) * t69;
	t59 = (t98 * t111 - t71 * t96) * t78 * t66;
	t58 = (t79 * t63 - (t100 * t60 - t67 * t77 + t68 * t83) * t112) * t62;
	t1 = [0, t61, 0, t60, t60, 0; 0, (-t86 * t89 * t63 - ((t68 * t104 - t67 * t84) * t89 + t100 * t61) * t112) * t62, 0, t58, t58, 0; 0, ((-t90 * t109 - t87 * t98) * t71 - (-t90 * t108 + t87 * t96) * t111) * t66, 0, t59, t59, t101 * t66;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end