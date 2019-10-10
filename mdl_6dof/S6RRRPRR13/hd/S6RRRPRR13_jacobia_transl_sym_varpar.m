% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR13
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:14
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR13_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR13_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR13_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (19->14), mult. (46->24), div. (0->0), fcn. (54->6), ass. (0->16)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t15 = t8 * t7;
	t9 = cos(qJ(2));
	t14 = t8 * t9;
	t10 = cos(qJ(1));
	t13 = t10 * t7;
	t12 = t10 * t9;
	t5 = sin(pkin(6));
	t11 = (pkin(9) + r_i_i_C(3)) * t5;
	t6 = cos(pkin(6));
	t4 = -t6 * t15 + t12;
	t3 = -t6 * t14 - t13;
	t2 = -t6 * t13 - t14;
	t1 = -t6 * t12 + t15;
	t16 = [-t8 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t10 * t11, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0, 0; t10 * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t11 * t8, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0; 0, (r_i_i_C(1) * t9 - r_i_i_C(2) * t7) * t5, 0, 0, 0, 0;];
	Ja_transl = t16;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (93->46), mult. (252->84), div. (0->0), fcn. (316->10), ass. (0->33)
	t39 = pkin(10) + r_i_i_C(3);
	t12 = cos(pkin(7));
	t13 = sin(qJ(3));
	t16 = cos(qJ(3));
	t10 = sin(pkin(7));
	t11 = sin(pkin(6));
	t18 = cos(qJ(1));
	t35 = t11 * t18;
	t27 = t10 * t35;
	t14 = sin(qJ(2));
	t15 = sin(qJ(1));
	t17 = cos(qJ(2));
	t28 = cos(pkin(6));
	t24 = t18 * t28;
	t3 = t15 * t14 - t17 * t24;
	t4 = t14 * t24 + t15 * t17;
	t38 = (t12 * t3 + t27) * t16 + t4 * t13;
	t36 = t11 * t15;
	t34 = t12 * t13;
	t33 = t12 * t16;
	t32 = t13 * t14;
	t31 = t13 * t17;
	t30 = t14 * t16;
	t29 = t16 * t17;
	t26 = t10 * t39;
	t25 = t15 * t28;
	t5 = -t18 * t14 - t17 * t25;
	t22 = t10 * t36 + t12 * t5;
	t19 = t13 * t27 - t4 * t16 + t3 * t34;
	t6 = -t14 * t25 + t18 * t17;
	t2 = t22 * t13 + t6 * t16;
	t1 = -t6 * t13 + t22 * t16;
	t7 = [t19 * r_i_i_C(1) + t38 * r_i_i_C(2) - t4 * pkin(2) - t15 * pkin(1) + pkin(9) * t35 + t39 * (-t3 * t10 + t12 * t35), (t5 * t16 - t6 * t34) * r_i_i_C(1) + (-t5 * t13 - t6 * t33) * r_i_i_C(2) + t5 * pkin(2) + t6 * t26, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; t18 * pkin(1) + t6 * pkin(2) + pkin(9) * t36 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t39 * (-t5 * t10 + t12 * t36), (-t3 * t16 - t4 * t34) * r_i_i_C(1) + (t3 * t13 - t4 * t33) * r_i_i_C(2) - t3 * pkin(2) + t4 * t26, -t38 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, ((-t12 * t32 + t29) * r_i_i_C(1) + (-t12 * t30 - t31) * r_i_i_C(2) + t17 * pkin(2) + t14 * t26) * t11, ((t12 * t29 - t32) * r_i_i_C(1) + (-t12 * t31 - t30) * r_i_i_C(2)) * t11 + (t16 * r_i_i_C(1) - t13 * r_i_i_C(2)) * t10 * t28, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (213->69), mult. (584->118), div. (0->0), fcn. (752->12), ass. (0->45)
	t30 = sin(pkin(7));
	t59 = t30 * pkin(10);
	t35 = sin(qJ(2));
	t36 = sin(qJ(1));
	t38 = cos(qJ(2));
	t39 = cos(qJ(1));
	t46 = cos(pkin(6));
	t42 = t39 * t46;
	t20 = t35 * t42 + t36 * t38;
	t34 = sin(qJ(3));
	t58 = t20 * t34;
	t29 = sin(pkin(13));
	t57 = t29 * t30;
	t32 = cos(pkin(13));
	t56 = t30 * t32;
	t31 = sin(pkin(6));
	t55 = t31 * t36;
	t54 = t31 * t39;
	t33 = cos(pkin(7));
	t53 = t33 * t34;
	t37 = cos(qJ(3));
	t52 = t33 * t37;
	t51 = t34 * t35;
	t50 = t34 * t38;
	t49 = t35 * t37;
	t48 = t37 * t38;
	t47 = r_i_i_C(3) + qJ(4);
	t45 = t30 * t55;
	t44 = t30 * t54;
	t43 = t36 * t46;
	t41 = t46 * t30;
	t40 = t32 * r_i_i_C(1) - t29 * r_i_i_C(2) + pkin(3);
	t19 = t36 * t35 - t38 * t42;
	t13 = -t19 * t30 + t33 * t54;
	t21 = -t39 * t35 - t38 * t43;
	t14 = -t21 * t30 + t33 * t55;
	t4 = t19 * t53 - t20 * t37 + t34 * t44;
	t22 = -t35 * t43 + t39 * t38;
	t11 = -t37 * t41 + (-t33 * t48 + t51) * t31;
	t10 = t21 * t37 - t22 * t53;
	t8 = -t19 * t37 - t20 * t53;
	t6 = t22 * t37 + (t21 * t33 + t45) * t34;
	t5 = -t21 * t52 + t22 * t34 - t37 * t45;
	t1 = t19 * t52 + t37 * t44 + t58;
	t2 = [(t13 * t29 + t4 * t32) * r_i_i_C(1) + (t13 * t32 - t4 * t29) * r_i_i_C(2) + t4 * pkin(3) - t20 * pkin(2) - t36 * pkin(1) + pkin(9) * t54 + t47 * (-t58 + (-t19 * t33 - t44) * t37) + t13 * pkin(10), (t10 * t32 + t22 * t57) * r_i_i_C(1) + (-t10 * t29 + t22 * t56) * r_i_i_C(2) + t10 * pkin(3) + t21 * pkin(2) + t22 * t59 + t47 * (t21 * t34 + t22 * t52), -t40 * t5 + t47 * t6, t5, 0, 0; (t14 * t29 + t6 * t32) * r_i_i_C(1) + (t14 * t32 - t6 * t29) * r_i_i_C(2) + t6 * pkin(3) + t22 * pkin(2) + t39 * pkin(1) + pkin(9) * t55 + t47 * t5 + t14 * pkin(10), (t20 * t57 + t8 * t32) * r_i_i_C(1) + (t20 * t56 - t8 * t29) * r_i_i_C(2) + t8 * pkin(3) - t19 * pkin(2) + t20 * t59 + t47 * (-t19 * t34 + t20 * t52), -t40 * t1 - t47 * t4, t1, 0, 0; 0, (t40 * (-t33 * t51 + t48) + t38 * pkin(2) + (r_i_i_C(1) * t29 + r_i_i_C(2) * t32 + pkin(10)) * t35 * t30 + t47 * (t33 * t49 + t50)) * t31, t47 * (t34 * t41 + (t33 * t50 + t49) * t31) - t40 * t11, t11, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:18
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (317->68), mult. (734->116), div. (0->0), fcn. (944->14), ass. (0->49)
	t68 = pkin(10) + pkin(4) * sin(pkin(13));
	t33 = cos(pkin(13)) * pkin(4) + pkin(3);
	t36 = pkin(13) + qJ(5);
	t34 = sin(t36);
	t35 = cos(t36);
	t49 = t35 * r_i_i_C(1) - t34 * r_i_i_C(2) + t33;
	t67 = t34 * r_i_i_C(1) + t35 * r_i_i_C(2) + t68;
	t65 = r_i_i_C(3) + pkin(11) + qJ(4);
	t43 = sin(qJ(2));
	t44 = sin(qJ(1));
	t46 = cos(qJ(2));
	t47 = cos(qJ(1));
	t55 = cos(pkin(6));
	t51 = t47 * t55;
	t24 = t43 * t51 + t44 * t46;
	t42 = sin(qJ(3));
	t64 = t24 * t42;
	t39 = sin(pkin(6));
	t63 = t39 * t44;
	t62 = t39 * t47;
	t40 = cos(pkin(7));
	t61 = t40 * t42;
	t45 = cos(qJ(3));
	t60 = t40 * t45;
	t59 = t42 * t43;
	t58 = t42 * t46;
	t57 = t43 * t45;
	t56 = t45 * t46;
	t38 = sin(pkin(7));
	t54 = t38 * t63;
	t53 = t38 * t62;
	t52 = t44 * t55;
	t50 = t55 * t38;
	t23 = t44 * t43 - t46 * t51;
	t15 = -t23 * t38 + t40 * t62;
	t25 = -t47 * t43 - t46 * t52;
	t17 = -t25 * t38 + t40 * t63;
	t6 = t23 * t61 - t24 * t45 + t42 * t53;
	t48 = t67 * t38;
	t26 = -t43 * t52 + t47 * t46;
	t22 = -t39 * t46 * t38 + t55 * t40;
	t14 = t42 * t50 + (t40 * t58 + t57) * t39;
	t13 = -t45 * t50 + (-t40 * t56 + t59) * t39;
	t8 = t26 * t45 + (t25 * t40 + t54) * t42;
	t7 = -t25 * t60 + t26 * t42 - t45 * t54;
	t3 = t23 * t60 + t45 * t53 + t64;
	t2 = t17 * t34 + t8 * t35;
	t1 = t17 * t35 - t8 * t34;
	t4 = [-t24 * pkin(2) - t44 * pkin(1) + pkin(9) * t62 + t65 * (-t64 + (-t23 * t40 - t53) * t45) + t49 * t6 + t67 * t15, t25 * pkin(2) + t49 * (t25 * t45 - t26 * t61) + t26 * t48 + t65 * (t25 * t42 + t26 * t60), -t49 * t7 + t65 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t47 * pkin(1) + t26 * pkin(2) + pkin(9) * t63 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t68 * t17 + t8 * t33 + t65 * t7, -t23 * pkin(2) + t65 * (-t23 * t42 + t24 * t60) + t49 * (-t23 * t45 - t24 * t61) + t24 * t48, -t49 * t3 - t6 * t65, t3, (-t15 * t35 + t6 * t34) * r_i_i_C(1) + (t15 * t34 + t6 * t35) * r_i_i_C(2), 0; 0, (t65 * (t40 * t57 + t58) + t49 * (-t40 * t59 + t56) + t46 * pkin(2) + t43 * t48) * t39, -t49 * t13 + t65 * t14, t13, (-t14 * t34 + t22 * t35) * r_i_i_C(1) + (-t14 * t35 - t22 * t34) * r_i_i_C(2), 0;];
	Ja_transl = t4;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:18
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (679->108), mult. (1554->181), div. (0->0), fcn. (2029->16), ass. (0->73)
	t56 = sin(qJ(2));
	t58 = cos(qJ(2));
	t59 = cos(qJ(1));
	t77 = cos(pkin(6));
	t70 = t59 * t77;
	t85 = sin(qJ(1));
	t37 = t85 * t56 - t58 * t70;
	t38 = t56 * t70 + t85 * t58;
	t55 = sin(qJ(3));
	t52 = cos(pkin(7));
	t86 = cos(qJ(3));
	t73 = t52 * t86;
	t50 = sin(pkin(7));
	t51 = sin(pkin(6));
	t81 = t51 * t59;
	t75 = t50 * t81;
	t15 = t37 * t73 + t38 * t55 + t86 * t75;
	t80 = t52 * t55;
	t16 = -t37 * t80 + t38 * t86 - t55 * t75;
	t30 = -t37 * t50 + t52 * t81;
	t48 = pkin(13) + qJ(5);
	t46 = sin(t48);
	t47 = cos(t48);
	t4 = t16 * t47 - t30 * t46;
	t54 = sin(qJ(6));
	t57 = cos(qJ(6));
	t94 = -t15 * t57 + t4 * t54;
	t93 = -t15 * t54 - t4 * t57;
	t90 = pkin(10) + pkin(4) * sin(pkin(13));
	t45 = cos(pkin(13)) * pkin(4) + pkin(3);
	t88 = r_i_i_C(3) + pkin(12);
	t89 = t88 * t46 + t45;
	t84 = t46 * t50;
	t83 = t47 * t50;
	t82 = t50 * t51;
	t79 = t55 * t56;
	t78 = t55 * t58;
	t76 = t56 * t82;
	t74 = t51 * t85;
	t72 = t86 * t56;
	t71 = t86 * t58;
	t69 = t77 * t50;
	t68 = t50 * t74;
	t67 = t90 * t50;
	t66 = t77 * t85;
	t65 = t57 * r_i_i_C(1) - t54 * r_i_i_C(2) + pkin(5);
	t53 = -pkin(11) - qJ(4);
	t64 = t54 * r_i_i_C(1) + t57 * r_i_i_C(2) - t53;
	t63 = t59 * t56 + t58 * t66;
	t62 = t63 * t86;
	t61 = -t65 * t47 - t89;
	t60 = t63 * t50 + t52 * t74;
	t39 = -t56 * t66 + t59 * t58;
	t36 = t77 * t52 - t58 * t82;
	t35 = (-t52 * t79 + t71) * t51;
	t34 = (t52 * t72 + t78) * t51;
	t29 = t55 * t69 + (t52 * t78 + t72) * t51;
	t28 = -t86 * t69 + (-t52 * t71 + t79) * t51;
	t26 = -t39 * t80 - t62;
	t25 = t39 * t73 - t63 * t55;
	t24 = -t37 * t86 - t38 * t80;
	t23 = -t37 * t55 + t38 * t73;
	t22 = t35 * t47 + t46 * t76;
	t20 = t39 * t86 + (-t63 * t52 + t68) * t55;
	t19 = t39 * t55 + t52 * t62 - t86 * t68;
	t14 = t29 * t47 + t36 * t46;
	t12 = t26 * t47 + t39 * t84;
	t10 = t24 * t47 + t38 * t84;
	t8 = t20 * t47 + t60 * t46;
	t7 = t20 * t46 - t60 * t47;
	t2 = t19 * t54 + t8 * t57;
	t1 = t19 * t57 - t8 * t54;
	t3 = [t93 * r_i_i_C(1) + t94 * r_i_i_C(2) - t4 * pkin(5) + t15 * t53 - t38 * pkin(2) - t85 * pkin(1) + pkin(9) * t81 - t89 * t16 + (-t88 * t47 + t90) * t30, (t12 * t57 + t25 * t54) * r_i_i_C(1) + (-t12 * t54 + t25 * t57) * r_i_i_C(2) + t12 * pkin(5) + t26 * t45 - t25 * t53 - t63 * pkin(2) + t39 * t67 + t88 * (t26 * t46 - t39 * t83), t61 * t19 + t64 * t20, t19, -t65 * t7 + t88 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t59 * pkin(1) + t39 * pkin(2) + t8 * pkin(5) + pkin(9) * t74 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t19 * t53 + t20 * t45 + t90 * t60 + t88 * t7, (t10 * t57 + t23 * t54) * r_i_i_C(1) + (-t10 * t54 + t23 * t57) * r_i_i_C(2) + t10 * pkin(5) + t24 * t45 - t23 * t53 - t37 * pkin(2) + t88 * (t24 * t46 - t38 * t83) + t38 * t67, t61 * t15 + t64 * t16, t15, t88 * t4 + t65 * (-t16 * t46 - t30 * t47), -t94 * r_i_i_C(1) + t93 * r_i_i_C(2); 0, (t22 * t57 + t34 * t54) * r_i_i_C(1) + (-t22 * t54 + t34 * t57) * r_i_i_C(2) + t22 * pkin(5) + t35 * t45 - t34 * t53 + t88 * (t35 * t46 - t47 * t76) + (t58 * pkin(2) + t56 * t67) * t51, t61 * t28 + t64 * t29, t28, t88 * t14 + t65 * (-t29 * t46 + t36 * t47), (-t14 * t54 + t28 * t57) * r_i_i_C(1) + (-t14 * t57 - t28 * t54) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end