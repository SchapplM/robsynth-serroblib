% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR2
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
%   Wie in S6PRRRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (101->18), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->29)
	t41 = sin(pkin(11));
	t42 = sin(pkin(6));
	t53 = t41 * t42;
	t48 = cos(qJ(2));
	t52 = t42 * t48;
	t44 = cos(pkin(6));
	t46 = sin(qJ(2));
	t51 = t44 * t46;
	t50 = t44 * t48;
	t43 = cos(pkin(11));
	t37 = -t41 * t51 + t43 * t48;
	t45 = sin(qJ(3));
	t47 = cos(qJ(3));
	t28 = t37 * t47 + t45 * t53;
	t26 = 0.1e1 / t28 ^ 2;
	t27 = t37 * t45 - t47 * t53;
	t49 = t27 ^ 2 * t26 + 0.1e1;
	t40 = 0.1e1 / t48 ^ 2;
	t36 = t41 * t50 + t43 * t46;
	t35 = t41 * t48 + t43 * t51;
	t33 = t41 * t46 - t43 * t50;
	t31 = atan2(-t33, -t52);
	t30 = cos(t31);
	t29 = sin(t31);
	t25 = 0.1e1 / t49;
	t24 = -t29 * t33 - t30 * t52;
	t23 = 0.1e1 / t24 ^ 2;
	t21 = (t35 / t48 + t46 * t33 * t40) / t42 / (0.1e1 + t33 ^ 2 / t42 ^ 2 * t40);
	t1 = [0, t21, 0, 0, 0, 0; 0, (t37 / t24 - (t30 * t42 * t46 - t29 * t35 + (t29 * t52 - t30 * t33) * t21) * t36 * t23) / (t36 ^ 2 * t23 + 0.1e1), 0, 0, 0, 0; 0, (-t45 / t28 + t47 * t27 * t26) * t36 * t25, t49 * t25, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (162->19), mult. (316->43), div. (52->11), fcn. (458->11), ass. (0->31)
	t59 = sin(pkin(11));
	t60 = sin(pkin(6));
	t69 = t59 * t60;
	t64 = cos(qJ(2));
	t68 = t60 * t64;
	t62 = cos(pkin(6));
	t63 = sin(qJ(2));
	t67 = t62 * t63;
	t66 = t62 * t64;
	t61 = cos(pkin(11));
	t52 = -t59 * t67 + t61 * t64;
	t58 = qJ(3) + qJ(4);
	t54 = sin(t58);
	t55 = cos(t58);
	t43 = t52 * t55 + t54 * t69;
	t41 = 0.1e1 / t43 ^ 2;
	t42 = t52 * t54 - t55 * t69;
	t65 = t42 ^ 2 * t41 + 0.1e1;
	t57 = 0.1e1 / t64 ^ 2;
	t51 = t59 * t66 + t61 * t63;
	t50 = t59 * t64 + t61 * t67;
	t48 = t59 * t63 - t61 * t66;
	t46 = atan2(-t48, -t68);
	t45 = cos(t46);
	t44 = sin(t46);
	t40 = 0.1e1 / t65;
	t39 = -t44 * t48 - t45 * t68;
	t38 = 0.1e1 / t39 ^ 2;
	t36 = (t50 / t64 + t63 * t48 * t57) / t60 / (0.1e1 + t48 ^ 2 / t60 ^ 2 * t57);
	t35 = t65 * t40;
	t1 = [0, t36, 0, 0, 0, 0; 0, (t52 / t39 - (t45 * t60 * t63 - t44 * t50 + (t44 * t68 - t45 * t48) * t36) * t51 * t38) / (t51 ^ 2 * t38 + 0.1e1), 0, 0, 0, 0; 0, (-t54 / t43 + t55 * t42 * t41) * t51 * t40, t35, t35, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (903->31), mult. (1304->78), div. (90->9), fcn. (1864->13), ass. (0->53)
	t89 = sin(pkin(6));
	t91 = cos(pkin(11));
	t100 = t89 * t91;
	t88 = sin(pkin(11));
	t94 = cos(qJ(2));
	t92 = cos(pkin(6));
	t93 = sin(qJ(2));
	t97 = t92 * t93;
	t80 = t88 * t94 + t91 * t97;
	t86 = qJ(3) + qJ(4);
	t84 = sin(t86);
	t85 = cos(t86);
	t70 = t85 * t100 + t80 * t84;
	t99 = t89 * t93;
	t77 = t84 * t99 - t92 * t85;
	t69 = atan2(-t70, t77);
	t66 = sin(t69);
	t67 = cos(t69);
	t60 = -t66 * t70 + t67 * t77;
	t59 = 0.1e1 / t60 ^ 2;
	t101 = t88 * t89;
	t82 = -t88 * t97 + t91 * t94;
	t73 = -t85 * t101 + t82 * t84;
	t106 = t59 * t73;
	t96 = t92 * t94;
	t81 = t88 * t96 + t91 * t93;
	t87 = sin(pkin(12));
	t103 = t81 * t87;
	t74 = t84 * t101 + t82 * t85;
	t90 = cos(pkin(12));
	t65 = t74 * t90 + t103;
	t63 = 0.1e1 / t65 ^ 2;
	t102 = t81 * t90;
	t64 = t74 * t87 - t102;
	t105 = t63 * t64;
	t76 = 0.1e1 / t77 ^ 2;
	t104 = t70 * t76;
	t98 = t89 * t94;
	t95 = -t66 * t77 - t67 * t70;
	t79 = -t88 * t93 + t91 * t96;
	t78 = t92 * t84 + t85 * t99;
	t75 = 0.1e1 / t77;
	t72 = -t84 * t100 + t80 * t85;
	t68 = 0.1e1 / (t70 ^ 2 * t76 + 0.1e1);
	t62 = 0.1e1 / t65;
	t61 = 0.1e1 / (t64 ^ 2 * t63 + 0.1e1);
	t58 = 0.1e1 / t60;
	t57 = 0.1e1 / (t73 ^ 2 * t59 + 0.1e1);
	t56 = (t98 * t104 - t75 * t79) * t84 * t68;
	t55 = (t78 * t104 - t72 * t75) * t68;
	t54 = (t90 * t105 - t62 * t87) * t73 * t61;
	t53 = (t74 * t58 - (t95 * t55 - t66 * t72 + t67 * t78) * t106) * t57;
	t1 = [0, t56, t55, t55, 0, 0; 0, (-t81 * t84 * t58 - ((-t66 * t79 + t67 * t98) * t84 + t95 * t56) * t106) * t57, t53, t53, 0, 0; 0, ((-t85 * t103 - t82 * t90) * t62 - (-t85 * t102 + t82 * t87) * t105) * t61, t54, t54, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (1012->32), mult. (1387->80), div. (95->9), fcn. (1976->13), ass. (0->54)
	t102 = sin(pkin(6));
	t103 = cos(pkin(11));
	t113 = t102 * t103;
	t101 = sin(pkin(11));
	t106 = cos(qJ(2));
	t104 = cos(pkin(6));
	t105 = sin(qJ(2));
	t110 = t104 * t105;
	t91 = t101 * t106 + t103 * t110;
	t100 = qJ(3) + qJ(4);
	t97 = sin(t100);
	t98 = cos(t100);
	t81 = t98 * t113 + t91 * t97;
	t112 = t102 * t105;
	t88 = -t104 * t98 + t97 * t112;
	t80 = atan2(-t81, t88);
	t77 = sin(t80);
	t78 = cos(t80);
	t71 = -t77 * t81 + t78 * t88;
	t70 = 0.1e1 / t71 ^ 2;
	t114 = t101 * t102;
	t93 = -t101 * t110 + t103 * t106;
	t84 = -t98 * t114 + t93 * t97;
	t118 = t70 * t84;
	t85 = t97 * t114 + t93 * t98;
	t109 = t104 * t106;
	t92 = t101 * t109 + t103 * t105;
	t99 = pkin(12) + qJ(6);
	t95 = sin(t99);
	t96 = cos(t99);
	t76 = t85 * t96 + t92 * t95;
	t74 = 0.1e1 / t76 ^ 2;
	t75 = t85 * t95 - t92 * t96;
	t117 = t74 * t75;
	t87 = 0.1e1 / t88 ^ 2;
	t116 = t81 * t87;
	t115 = t92 * t98;
	t111 = t102 * t106;
	t108 = t75 ^ 2 * t74 + 0.1e1;
	t107 = -t77 * t88 - t78 * t81;
	t90 = -t101 * t105 + t103 * t109;
	t89 = t104 * t97 + t98 * t112;
	t86 = 0.1e1 / t88;
	t83 = -t97 * t113 + t91 * t98;
	t79 = 0.1e1 / (t81 ^ 2 * t87 + 0.1e1);
	t73 = 0.1e1 / t76;
	t72 = 0.1e1 / t108;
	t69 = 0.1e1 / t71;
	t68 = 0.1e1 / (t84 ^ 2 * t70 + 0.1e1);
	t67 = (t111 * t116 - t86 * t90) * t97 * t79;
	t66 = (t89 * t116 - t83 * t86) * t79;
	t65 = (t96 * t117 - t73 * t95) * t84 * t72;
	t64 = (t85 * t69 - (t107 * t66 - t77 * t83 + t78 * t89) * t118) * t68;
	t1 = [0, t67, t66, t66, 0, 0; 0, (-t92 * t97 * t69 - ((t78 * t111 - t77 * t90) * t97 + t107 * t67) * t118) * t68, t64, t64, 0, 0; 0, ((-t95 * t115 - t93 * t96) * t73 - (-t96 * t115 + t93 * t95) * t117) * t72, t65, t65, 0, t108 * t72;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end