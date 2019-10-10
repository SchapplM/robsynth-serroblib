% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S7RRRRRRR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% Ja_transl [3x7]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 17:10
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S7RRRRRRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(7,1),uint8(0),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobia_transl_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S7RRRRRRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobia_transl_sym_varpar: pkin has to be [4x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->4), mult. (20->10), div. (0->0), fcn. (20->4), ass. (0->7)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t6 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t5 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t7 = [t4 * r_i_i_C(3) - t2 * t6, t5 * t4, 0, 0, 0, 0, 0; t2 * r_i_i_C(3) + t4 * t6, t5 * t2, 0, 0, 0, 0, 0; 0, t6, 0, 0, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (26->15), mult. (70->31), div. (0->0), fcn. (78->6), ass. (0->17)
	t15 = pkin(2) + r_i_i_C(3);
	t6 = sin(qJ(2));
	t13 = t15 * t6;
	t7 = sin(qJ(1));
	t9 = cos(qJ(2));
	t16 = t7 * t9;
	t10 = cos(qJ(1));
	t14 = t10 * t9;
	t5 = sin(qJ(3));
	t8 = cos(qJ(3));
	t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5;
	t11 = -t12 * t6 - t15 * t9;
	t4 = t8 * t14 - t7 * t5;
	t3 = -t5 * t14 - t7 * t8;
	t2 = -t10 * t5 - t8 * t16;
	t1 = -t10 * t8 + t5 * t16;
	t17 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * t13, t11 * t10, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t10 * t13, t11 * t7, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0; 0, t12 * t9 - t13, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6, 0, 0, 0, 0;];
	Ja_transl = t17;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (56->31), mult. (157->63), div. (0->0), fcn. (189->8), ass. (0->29)
	t14 = cos(qJ(2));
	t9 = sin(qJ(3));
	t29 = t14 * t9;
	t15 = cos(qJ(1));
	t28 = t15 * t9;
	t10 = sin(qJ(2));
	t11 = sin(qJ(1));
	t27 = t10 * t11;
	t12 = cos(qJ(4));
	t26 = t10 * t12;
	t13 = cos(qJ(3));
	t25 = t10 * t13;
	t24 = t10 * t15;
	t23 = t13 * t14;
	t22 = t15 * t13;
	t8 = sin(qJ(4));
	t21 = t12 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t4 = t11 * t23 + t28;
	t20 = -t12 * t4 - t8 * t27;
	t19 = t11 * t26 - t4 * t8;
	t18 = t12 * t14 + t8 * t25;
	t17 = -t12 * t25 + t14 * t8;
	t16 = t10 * t9 * r_i_i_C(3) - t14 * pkin(2) + t17 * r_i_i_C(1) + t18 * r_i_i_C(2);
	t6 = -t11 * t9 + t14 * t22;
	t5 = -t11 * t13 - t14 * t28;
	t3 = t11 * t29 - t22;
	t2 = t6 * t12 + t8 * t24;
	t1 = t12 * t24 - t6 * t8;
	t7 = [pkin(2) * t27 + t20 * r_i_i_C(1) - t19 * r_i_i_C(2) + t3 * r_i_i_C(3), t16 * t15, -r_i_i_C(3) * t6 + t21 * t5, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; -pkin(2) * t24 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * r_i_i_C(3), t16 * t11, -r_i_i_C(3) * t4 - t21 * t3, t19 * r_i_i_C(1) + t20 * r_i_i_C(2), 0, 0, 0; 0, (t10 * t8 + t12 * t23) * r_i_i_C(1) + (-t8 * t23 + t26) * r_i_i_C(2) - r_i_i_C(3) * t29 - t10 * pkin(2), (-r_i_i_C(3) * t13 - t21 * t9) * t10, -t18 * r_i_i_C(1) + t17 * r_i_i_C(2), 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:04
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (141->54), mult. (389->104), div. (0->0), fcn. (489->10), ass. (0->47)
	t30 = cos(qJ(4));
	t24 = sin(qJ(5));
	t29 = cos(qJ(5));
	t38 = t29 * r_i_i_C(1) - t24 * r_i_i_C(2);
	t25 = sin(qJ(4));
	t56 = pkin(3) + r_i_i_C(3);
	t39 = t56 * t25;
	t60 = t38 * t30 + t39;
	t26 = sin(qJ(3));
	t31 = cos(qJ(3));
	t33 = cos(qJ(1));
	t42 = t33 * t31;
	t28 = sin(qJ(1));
	t32 = cos(qJ(2));
	t46 = t28 * t32;
	t15 = t26 * t46 - t42;
	t43 = t33 * t26;
	t16 = t31 * t46 + t43;
	t27 = sin(qJ(2));
	t50 = t27 * t28;
	t4 = t16 * t30 + t25 * t50;
	t59 = t15 * t29 + t4 * t24;
	t58 = t15 * t24 - t4 * t29;
	t53 = t26 * t27;
	t52 = t26 * t30;
	t51 = t26 * t32;
	t49 = t27 * t30;
	t48 = t27 * t31;
	t47 = t27 * t33;
	t45 = t32 * t25;
	t44 = t32 * t30;
	t41 = t24 * t53;
	t40 = t29 * t53;
	t37 = -t24 * r_i_i_C(1) - t29 * r_i_i_C(2);
	t36 = -t16 * t25 + t28 * t49;
	t14 = t30 * t48 - t45;
	t35 = t25 * t48 + t44;
	t34 = t56 * t35;
	t20 = -t28 * t26 + t32 * t42;
	t19 = -t28 * t31 - t32 * t43;
	t18 = t27 * t25 + t31 * t44;
	t10 = t14 * t28;
	t8 = t20 * t30 + t25 * t47;
	t7 = t20 * t25 - t30 * t47;
	t2 = t19 * t24 + t8 * t29;
	t1 = t19 * t29 - t8 * t24;
	t3 = [pkin(2) * t50 + t58 * r_i_i_C(1) + t59 * r_i_i_C(2) + t56 * t36, (-t32 * pkin(2) + t41 * r_i_i_C(1) + t40 * r_i_i_C(2) - t38 * t14 - t34) * t33, t60 * t19 + t37 * t20, -t38 * t7 + t56 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; -pkin(2) * t47 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t56 * t7, (-t10 * t29 + t28 * t41) * r_i_i_C(1) + (t10 * t24 + t28 * t40) * r_i_i_C(2) - pkin(2) * t46 - t28 * t34, -t60 * t15 + t37 * t16, t38 * t36 + t56 * t4, -t59 * r_i_i_C(1) + t58 * r_i_i_C(2), 0, 0; 0, (t18 * t29 - t24 * t51) * r_i_i_C(1) + (-t18 * t24 - t29 * t51) * r_i_i_C(2) - t27 * pkin(2) + t56 * (t31 * t45 - t49), ((-t24 * t31 - t29 * t52) * r_i_i_C(1) + (t24 * t52 - t29 * t31) * r_i_i_C(2) - t26 * t39) * t27, t56 * t14 - t38 * t35, (-t14 * t24 - t40) * r_i_i_C(1) + (-t14 * t29 + t41) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:05
	% DurationCPUTime: 0.73s
	% Computational Cost: add. (255->98), mult. (711->189), div. (0->0), fcn. (919->12), ass. (0->65)
	t48 = cos(qJ(3));
	t42 = sin(qJ(3));
	t50 = cos(qJ(1));
	t55 = t50 * t42;
	t44 = sin(qJ(1));
	t49 = cos(qJ(2));
	t60 = t44 * t49;
	t30 = t48 * t60 + t55;
	t41 = sin(qJ(4));
	t43 = sin(qJ(2));
	t47 = cos(qJ(4));
	t63 = t43 * t47;
	t14 = t30 * t41 - t44 * t63;
	t39 = sin(qJ(6));
	t64 = t43 * t44;
	t15 = t30 * t47 + t41 * t64;
	t54 = t50 * t48;
	t29 = t42 * t60 - t54;
	t40 = sin(qJ(5));
	t46 = cos(qJ(5));
	t4 = t15 * t46 - t29 * t40;
	t45 = cos(qJ(6));
	t79 = -t14 * t45 + t4 * t39;
	t78 = -t14 * t39 - t4 * t45;
	t75 = -t15 * t40 - t29 * t46;
	t74 = t40 * r_i_i_C(3);
	t73 = t41 * pkin(3);
	t70 = t39 * t41;
	t69 = t39 * t46;
	t68 = t40 * t47;
	t67 = t41 * t45;
	t66 = t42 * t43;
	t65 = t42 * t49;
	t62 = t43 * t48;
	t61 = t43 * t50;
	t59 = t45 * t46;
	t58 = t46 * t47;
	t57 = t49 * t41;
	t56 = t49 * t47;
	t53 = t40 * t66;
	t52 = t46 * t66;
	t51 = t45 * r_i_i_C(1) - t39 * r_i_i_C(2);
	t28 = t47 * t62 - t57;
	t27 = t41 * t62 + t56;
	t34 = -t44 * t42 + t49 * t54;
	t33 = -t44 * t48 - t49 * t55;
	t32 = t43 * t41 + t48 * t56;
	t31 = t48 * t57 - t63;
	t25 = t28 * t50;
	t24 = t27 * t50;
	t23 = t28 * t44;
	t22 = t27 * t44;
	t20 = t34 * t47 + t41 * t61;
	t19 = t34 * t41 - t47 * t61;
	t18 = t32 * t46 - t40 * t65;
	t13 = t28 * t46 - t53;
	t11 = -t25 * t46 + t50 * t53;
	t10 = -t23 * t46 + t44 * t53;
	t9 = t33 * t58 - t34 * t40;
	t8 = -t29 * t58 - t30 * t40;
	t7 = t20 * t46 + t33 * t40;
	t6 = t20 * t40 - t33 * t46;
	t2 = t19 * t39 + t7 * t45;
	t1 = t19 * t45 - t7 * t39;
	t3 = [pkin(2) * t64 - t14 * pkin(3) + t78 * r_i_i_C(1) + t79 * r_i_i_C(2) + t75 * r_i_i_C(3), (t11 * t45 - t24 * t39) * r_i_i_C(1) + (-t11 * t39 - t24 * t45) * r_i_i_C(2) + (-t25 * t40 - t50 * t52) * r_i_i_C(3) - t24 * pkin(3) - t50 * t49 * pkin(2), (t33 * t70 + t9 * t45) * r_i_i_C(1) + (t33 * t67 - t9 * t39) * r_i_i_C(2) + (t33 * t68 + t34 * t46) * r_i_i_C(3) + t33 * t73, (-t19 * t59 + t20 * t39) * r_i_i_C(1) + (t19 * t69 + t20 * t45) * r_i_i_C(2) - t19 * t74 + t20 * pkin(3), t7 * r_i_i_C(3) - t51 * t6, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; -pkin(2) * t61 + t19 * pkin(3) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t6 * r_i_i_C(3), (t10 * t45 - t22 * t39) * r_i_i_C(1) + (-t10 * t39 - t22 * t45) * r_i_i_C(2) + (-t23 * t40 - t44 * t52) * r_i_i_C(3) - t22 * pkin(3) - pkin(2) * t60, (-t29 * t70 + t8 * t45) * r_i_i_C(1) + (-t29 * t67 - t8 * t39) * r_i_i_C(2) + (-t29 * t68 + t30 * t46) * r_i_i_C(3) - t29 * t73, (-t14 * t59 + t15 * t39) * r_i_i_C(1) + (t14 * t69 + t15 * t45) * r_i_i_C(2) - t14 * t74 + t15 * pkin(3), t4 * r_i_i_C(3) + t51 * t75, -t79 * r_i_i_C(1) + t78 * r_i_i_C(2), 0; 0, (t18 * t45 + t31 * t39) * r_i_i_C(1) + (-t18 * t39 + t31 * t45) * r_i_i_C(2) + (t32 * t40 + t46 * t65) * r_i_i_C(3) + t31 * pkin(3) - t43 * pkin(2), (t51 * (-t40 * t48 - t42 * t58) + (-t42 * t68 + t46 * t48) * r_i_i_C(3) + (-t39 * r_i_i_C(1) - t45 * r_i_i_C(2) - pkin(3)) * t42 * t41) * t43, (-t27 * t59 + t28 * t39) * r_i_i_C(1) + (t27 * t69 + t28 * t45) * r_i_i_C(2) - t27 * t74 + t28 * pkin(3), t13 * r_i_i_C(3) + t51 * (-t28 * t40 - t52), (-t13 * t39 + t27 * t45) * r_i_i_C(1) + (-t13 * t45 - t27 * t39) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 17:10:04
	% EndTime: 2019-10-10 17:10:05
	% DurationCPUTime: 1.41s
	% Computational Cost: add. (562->129), mult. (1558->252), div. (0->0), fcn. (2050->14), ass. (0->95)
	t78 = sin(qJ(2));
	t79 = sin(qJ(1));
	t103 = t78 * t79;
	t84 = cos(qJ(3));
	t77 = sin(qJ(3));
	t86 = cos(qJ(1));
	t94 = t86 * t77;
	t85 = cos(qJ(2));
	t99 = t79 * t85;
	t64 = t84 * t99 + t94;
	t76 = sin(qJ(4));
	t83 = cos(qJ(4));
	t46 = t76 * t103 + t64 * t83;
	t93 = t86 * t84;
	t63 = t77 * t99 - t93;
	t75 = sin(qJ(5));
	t82 = cos(qJ(5));
	t23 = t46 * t75 + t63 * t82;
	t24 = t46 * t82 - t63 * t75;
	t102 = t78 * t83;
	t45 = -t79 * t102 + t64 * t76;
	t74 = sin(qJ(6));
	t81 = cos(qJ(6));
	t4 = t24 * t81 + t45 * t74;
	t73 = sin(qJ(7));
	t80 = cos(qJ(7));
	t123 = t23 * t80 + t4 * t73;
	t122 = t23 * t73 - t4 * t80;
	t119 = t24 * t74 - t45 * t81;
	t116 = pkin(4) + r_i_i_C(3);
	t115 = pkin(3) * t76;
	t112 = t73 * t75;
	t111 = t74 * t76;
	t110 = t74 * t82;
	t109 = t75 * t80;
	t108 = t75 * t83;
	t107 = t76 * t81;
	t106 = t77 * t78;
	t105 = t77 * t85;
	t104 = t78 * t76;
	t101 = t78 * t84;
	t100 = t78 * t86;
	t98 = t81 * t82;
	t97 = t82 * t83;
	t96 = t85 * t76;
	t95 = t85 * t83;
	t92 = t77 * t104;
	t91 = t75 * t106;
	t90 = t82 * t106;
	t89 = t80 * r_i_i_C(1) - t73 * r_i_i_C(2);
	t88 = -t73 * r_i_i_C(1) - t80 * r_i_i_C(2);
	t62 = t83 * t101 - t96;
	t61 = t76 * t101 + t95;
	t87 = t116 * t74 - t89 * t81;
	t68 = -t79 * t77 + t85 * t93;
	t67 = -t79 * t84 - t85 * t94;
	t66 = t84 * t95 + t104;
	t65 = t84 * t96 - t102;
	t58 = t62 * t86;
	t57 = t61 * t86;
	t56 = t62 * t79;
	t55 = t61 * t79;
	t54 = (-t75 * t84 - t77 * t97) * t78;
	t53 = (-t77 * t108 + t82 * t84) * t78;
	t52 = t76 * t100 + t68 * t83;
	t51 = -t83 * t100 + t68 * t76;
	t50 = -t75 * t105 + t66 * t82;
	t49 = t82 * t105 + t66 * t75;
	t44 = t62 * t82 - t91;
	t43 = t62 * t75 + t90;
	t42 = -t58 * t82 + t86 * t91;
	t41 = -t58 * t75 - t86 * t90;
	t40 = -t56 * t82 + t79 * t91;
	t39 = -t56 * t75 - t79 * t90;
	t38 = t54 * t81 - t74 * t92;
	t36 = t67 * t97 - t68 * t75;
	t35 = t67 * t108 + t68 * t82;
	t34 = -t63 * t97 - t64 * t75;
	t33 = -t63 * t108 + t64 * t82;
	t32 = -t61 * t98 + t62 * t74;
	t30 = t52 * t82 + t67 * t75;
	t29 = t52 * t75 - t67 * t82;
	t28 = t50 * t81 + t65 * t74;
	t22 = t44 * t81 + t61 * t74;
	t20 = t42 * t81 - t57 * t74;
	t18 = t40 * t81 - t55 * t74;
	t16 = t67 * t111 + t36 * t81;
	t14 = -t63 * t111 + t34 * t81;
	t12 = -t51 * t98 + t52 * t74;
	t10 = -t45 * t98 + t46 * t74;
	t8 = t30 * t81 + t51 * t74;
	t7 = -t30 * t74 + t51 * t81;
	t2 = -t29 * t73 + t8 * t80;
	t1 = -t29 * t80 - t8 * t73;
	t3 = [pkin(2) * t103 - t45 * pkin(3) + t122 * r_i_i_C(1) + t123 * r_i_i_C(2) + t116 * t119, (t20 * t80 - t41 * t73) * r_i_i_C(1) + (-t20 * t73 - t41 * t80) * r_i_i_C(2) - t57 * pkin(3) - t86 * t85 * pkin(2) + t116 * (-t42 * t74 - t57 * t81), (t16 * t80 - t35 * t73) * r_i_i_C(1) + (-t16 * t73 - t35 * t80) * r_i_i_C(2) + t67 * t115 + t116 * (t67 * t107 - t36 * t74), (t51 * t112 + t12 * t80) * r_i_i_C(1) + (t51 * t109 - t12 * t73) * r_i_i_C(2) + t52 * pkin(3) + t116 * (t51 * t110 + t52 * t81), t87 * t29 + t88 * t30, -t116 * t8 + t89 * t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -pkin(2) * t100 + t51 * pkin(3) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t116 * t7, (t18 * t80 - t39 * t73) * r_i_i_C(1) + (-t18 * t73 - t39 * t80) * r_i_i_C(2) - t55 * pkin(3) - pkin(2) * t99 + t116 * (-t40 * t74 - t55 * t81), (t14 * t80 - t33 * t73) * r_i_i_C(1) + (-t14 * t73 - t33 * t80) * r_i_i_C(2) - t63 * t115 + t116 * (-t63 * t107 - t34 * t74), (t10 * t80 + t45 * t112) * r_i_i_C(1) + (-t10 * t73 + t45 * t109) * r_i_i_C(2) + t46 * pkin(3) + t116 * (t45 * t110 + t46 * t81), t87 * t23 + t88 * t24, -t116 * t4 - t89 * t119, -t123 * r_i_i_C(1) + t122 * r_i_i_C(2); 0, (t28 * t80 - t49 * t73) * r_i_i_C(1) + (-t28 * t73 - t49 * t80) * r_i_i_C(2) + t65 * pkin(3) - t78 * pkin(2) + t116 * (-t50 * t74 + t65 * t81), (t38 * t80 - t53 * t73) * r_i_i_C(1) + (-t38 * t73 - t53 * t80) * r_i_i_C(2) - pkin(3) * t92 + t116 * (-t54 * t74 - t81 * t92), (t61 * t112 + t32 * t80) * r_i_i_C(1) + (t61 * t109 - t32 * t73) * r_i_i_C(2) + t62 * pkin(3) + t116 * (t61 * t110 + t62 * t81), t87 * t43 + t88 * t44, -t116 * t22 + t89 * (-t44 * t74 + t61 * t81), (-t22 * t73 - t43 * t80) * r_i_i_C(1) + (-t22 * t80 + t43 * t73) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,7);
end