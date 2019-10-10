% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:19
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(13));
	t1 = sin(pkin(13));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (55->24), mult. (156->50), div. (0->0), fcn. (196->10), ass. (0->24)
	t7 = sin(pkin(7));
	t8 = sin(pkin(6));
	t24 = t7 * t8;
	t10 = cos(pkin(7));
	t15 = cos(qJ(2));
	t23 = t10 * t15;
	t11 = cos(pkin(6));
	t13 = sin(qJ(2));
	t22 = t11 * t13;
	t21 = t11 * t15;
	t12 = sin(qJ(3));
	t14 = cos(qJ(3));
	t20 = r_i_i_C(1) * t14 - r_i_i_C(2) * t12;
	t6 = sin(pkin(13));
	t9 = cos(pkin(13));
	t1 = -t6 * t13 + t9 * t21;
	t19 = -t1 * t10 + t9 * t24;
	t3 = -t9 * t13 - t6 * t21;
	t18 = t10 * t3 + t6 * t24;
	t17 = pkin(2) + t20;
	t16 = (pkin(9) + r_i_i_C(3)) * t7 + (-t12 * r_i_i_C(1) - t14 * r_i_i_C(2)) * t10;
	t4 = t15 * t9 - t6 * t22;
	t2 = t15 * t6 + t9 * t22;
	t5 = [0, t16 * t4 + t17 * t3, (-t12 * t4 + t18 * t14) * r_i_i_C(1) + (-t18 * t12 - t14 * t4) * r_i_i_C(2), 0, 0, 0; 0, t17 * t1 + t16 * t2, (-t12 * t2 - t19 * t14) * r_i_i_C(1) + (t19 * t12 - t14 * t2) * r_i_i_C(2), 0, 0, 0; 1, (t16 * t13 + t17 * t15) * t8, t20 * t7 * t11 + ((-t12 * t13 + t14 * t23) * r_i_i_C(1) + (-t12 * t23 - t13 * t14) * r_i_i_C(2)) * t8, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (168->57), mult. (471->110), div. (0->0), fcn. (606->12), ass. (0->43)
	t50 = r_i_i_C(3) + pkin(10);
	t22 = sin(pkin(7));
	t49 = t22 * pkin(9);
	t23 = sin(pkin(6));
	t48 = t22 * t23;
	t26 = cos(pkin(6));
	t47 = t22 * t26;
	t27 = sin(qJ(4));
	t46 = t22 * t27;
	t30 = cos(qJ(4));
	t45 = t22 * t30;
	t25 = cos(pkin(7));
	t44 = t23 * t25;
	t28 = sin(qJ(3));
	t43 = t25 * t28;
	t31 = cos(qJ(3));
	t42 = t25 * t31;
	t29 = sin(qJ(2));
	t41 = t26 * t29;
	t32 = cos(qJ(2));
	t40 = t26 * t32;
	t39 = t28 * t29;
	t38 = t28 * t32;
	t37 = t29 * t31;
	t36 = t31 * t32;
	t35 = t30 * r_i_i_C(1) - t27 * r_i_i_C(2) + pkin(3);
	t21 = sin(pkin(13));
	t24 = cos(pkin(13));
	t16 = -t21 * t29 + t24 * t40;
	t34 = t16 * t25 - t24 * t48;
	t18 = -t21 * t40 - t24 * t29;
	t33 = t18 * t25 + t21 * t48;
	t19 = -t21 * t41 + t24 * t32;
	t17 = t21 * t32 + t24 * t41;
	t15 = t26 * t25 - t32 * t48;
	t12 = -t18 * t22 + t21 * t44;
	t11 = -t16 * t22 - t24 * t44;
	t10 = t28 * t47 + (t25 * t38 + t37) * t23;
	t8 = t18 * t31 - t19 * t43;
	t6 = t16 * t31 - t17 * t43;
	t4 = t19 * t31 + t33 * t28;
	t2 = t17 * t31 + t34 * t28;
	t1 = [0, (t19 * t46 + t8 * t30) * r_i_i_C(1) + (t19 * t45 - t8 * t27) * r_i_i_C(2) + t8 * pkin(3) + t18 * pkin(2) + t19 * t49 + t50 * (t18 * t28 + t19 * t42), t50 * t4 + t35 * (-t19 * t28 + t33 * t31), (t12 * t30 - t4 * t27) * r_i_i_C(1) + (-t12 * t27 - t4 * t30) * r_i_i_C(2), 0, 0; 0, (t17 * t46 + t6 * t30) * r_i_i_C(1) + (t17 * t45 - t6 * t27) * r_i_i_C(2) + t6 * pkin(3) + t16 * pkin(2) + t17 * t49 + t50 * (t16 * t28 + t17 * t42), t50 * t2 + t35 * (-t17 * t28 + t34 * t31), (t11 * t30 - t2 * t27) * r_i_i_C(1) + (-t11 * t27 - t2 * t30) * r_i_i_C(2), 0, 0; 1, (t35 * (-t25 * t39 + t36) + t32 * pkin(2) + (t27 * r_i_i_C(1) + t30 * r_i_i_C(2) + pkin(9)) * t29 * t22 + t50 * (t25 * t37 + t38)) * t23, t50 * t10 + t35 * (t31 * t47 + (t25 * t36 - t39) * t23), (-t10 * t27 + t15 * t30) * r_i_i_C(1) + (-t10 * t30 - t15 * t27) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (287->58), mult. (661->106), div. (0->0), fcn. (849->14), ass. (0->45)
	t31 = sin(pkin(13));
	t34 = cos(pkin(13));
	t39 = sin(qJ(2));
	t36 = cos(pkin(6));
	t42 = cos(qJ(2));
	t52 = t36 * t42;
	t22 = -t31 * t39 + t34 * t52;
	t32 = sin(pkin(7));
	t33 = sin(pkin(6));
	t35 = cos(pkin(7));
	t56 = t33 * t35;
	t17 = -t22 * t32 - t34 * t56;
	t30 = qJ(4) + qJ(5);
	t28 = sin(t30);
	t29 = cos(t30);
	t53 = t36 * t39;
	t23 = t31 * t42 + t34 * t53;
	t38 = sin(qJ(3));
	t41 = cos(qJ(3));
	t58 = t32 * t33;
	t46 = t22 * t35 - t34 * t58;
	t8 = t23 * t41 + t46 * t38;
	t62 = (t17 * t29 - t8 * t28) * r_i_i_C(1) + (-t17 * t28 - t8 * t29) * r_i_i_C(2);
	t25 = -t31 * t53 + t34 * t42;
	t24 = -t31 * t52 - t34 * t39;
	t45 = t24 * t35 + t31 * t58;
	t10 = t25 * t41 + t45 * t38;
	t18 = -t24 * t32 + t31 * t56;
	t61 = (-t10 * t28 + t18 * t29) * r_i_i_C(1) + (-t10 * t29 - t18 * t28) * r_i_i_C(2);
	t49 = t39 * t41;
	t50 = t38 * t42;
	t57 = t32 * t36;
	t16 = t38 * t57 + (t35 * t50 + t49) * t33;
	t21 = t36 * t35 - t42 * t58;
	t60 = (-t16 * t28 + t21 * t29) * r_i_i_C(1) + (-t16 * t29 - t21 * t28) * r_i_i_C(2);
	t59 = r_i_i_C(3) + pkin(11) + pkin(10);
	t55 = t35 * t38;
	t54 = t35 * t41;
	t51 = t38 * t39;
	t48 = t41 * t42;
	t40 = cos(qJ(4));
	t47 = t40 * pkin(4) + t29 * r_i_i_C(1) - t28 * r_i_i_C(2) + pkin(3);
	t37 = sin(qJ(4));
	t44 = (pkin(4) * t37 + r_i_i_C(1) * t28 + r_i_i_C(2) * t29 + pkin(9)) * t32;
	t1 = [0, t24 * pkin(2) + t47 * (t24 * t41 - t25 * t55) + t25 * t44 + t59 * (t24 * t38 + t25 * t54), t59 * t10 + t47 * (-t25 * t38 + t45 * t41), (-t10 * t37 + t18 * t40) * pkin(4) + t61, t61, 0; 0, t22 * pkin(2) + t47 * (t22 * t41 - t23 * t55) + t23 * t44 + t59 * (t22 * t38 + t23 * t54), t59 * t8 + t47 * (-t23 * t38 + t46 * t41), (t17 * t40 - t37 * t8) * pkin(4) + t62, t62, 0; 1, (t59 * (t35 * t49 + t50) + t47 * (-t35 * t51 + t48) + t42 * pkin(2) + t39 * t44) * t33, t59 * t16 + t47 * (t41 * t57 + (t35 * t48 - t51) * t33), (-t16 * t37 + t21 * t40) * pkin(4) + t60, t60, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (660->94), mult. (1496->168), div. (0->0), fcn. (1949->16), ass. (0->67)
	t90 = r_i_i_C(3) + pkin(12);
	t52 = sin(qJ(6));
	t56 = cos(qJ(6));
	t93 = t56 * r_i_i_C(1) - t52 * r_i_i_C(2) + pkin(5);
	t49 = sin(pkin(7));
	t55 = sin(qJ(2));
	t58 = cos(qJ(2));
	t51 = cos(pkin(13));
	t81 = cos(pkin(6));
	t75 = t51 * t81;
	t79 = sin(pkin(13));
	t63 = t79 * t55 - t58 * t75;
	t80 = cos(pkin(7));
	t50 = sin(pkin(6));
	t84 = t50 * t51;
	t92 = t49 * t84 + t63 * t80;
	t69 = t81 * t79;
	t64 = t51 * t55 + t58 * t69;
	t76 = t50 * t79;
	t91 = -t49 * t76 + t64 * t80;
	t87 = cos(qJ(3));
	t48 = qJ(4) + qJ(5);
	t46 = sin(t48);
	t86 = t46 * t49;
	t47 = cos(t48);
	t85 = t47 * t49;
	t83 = t50 * t55;
	t82 = t50 * t58;
	t78 = t49 * t83;
	t54 = sin(qJ(3));
	t74 = t54 * t80;
	t73 = t81 * t49;
	t53 = sin(qJ(4));
	t71 = (pkin(4) * t53 + pkin(9)) * t49;
	t70 = t80 * t87;
	t59 = -pkin(11) - pkin(10);
	t68 = t52 * r_i_i_C(1) + t56 * r_i_i_C(2) - t59;
	t39 = t55 * t75 + t79 * t58;
	t21 = t39 * t87 - t92 * t54;
	t32 = t63 * t49 - t80 * t84;
	t8 = t21 * t47 + t32 * t46;
	t67 = t90 * t8 + t93 * (-t21 * t46 + t32 * t47);
	t40 = t51 * t58 - t55 * t69;
	t23 = t40 * t87 - t91 * t54;
	t33 = t64 * t49 + t80 * t76;
	t10 = t23 * t47 + t33 * t46;
	t66 = t90 * t10 + t93 * (-t23 * t46 + t33 * t47);
	t31 = t54 * t73 + (t87 * t55 + t58 * t74) * t50;
	t38 = -t49 * t82 + t81 * t80;
	t19 = t31 * t47 + t38 * t46;
	t65 = t90 * t19 + t93 * (-t31 * t46 + t38 * t47);
	t57 = cos(qJ(4));
	t45 = t57 * pkin(4) + pkin(3);
	t62 = -t90 * t46 - t93 * t47 - t45;
	t37 = (-t55 * t74 + t87 * t58) * t50;
	t36 = (t54 * t58 + t55 * t70) * t50;
	t30 = t54 * t83 - t70 * t82 - t87 * t73;
	t29 = t37 * t47 + t46 * t78;
	t27 = -t40 * t74 - t64 * t87;
	t26 = t40 * t70 - t64 * t54;
	t25 = -t39 * t74 - t63 * t87;
	t24 = t39 * t70 - t63 * t54;
	t22 = t40 * t54 + t91 * t87;
	t20 = t39 * t54 + t92 * t87;
	t14 = t27 * t47 + t40 * t86;
	t12 = t25 * t47 + t39 * t86;
	t1 = [0, (t14 * t56 + t26 * t52) * r_i_i_C(1) + (-t14 * t52 + t26 * t56) * r_i_i_C(2) + t14 * pkin(5) + t27 * t45 - t26 * t59 - t64 * pkin(2) + t40 * t71 + t90 * (t27 * t46 - t40 * t85), t62 * t22 + t68 * t23, (-t23 * t53 + t33 * t57) * pkin(4) + t66, t66, (-t10 * t52 + t22 * t56) * r_i_i_C(1) + (-t10 * t56 - t22 * t52) * r_i_i_C(2); 0, (t12 * t56 + t24 * t52) * r_i_i_C(1) + (-t12 * t52 + t24 * t56) * r_i_i_C(2) + t12 * pkin(5) + t25 * t45 - t24 * t59 - t63 * pkin(2) + t39 * t71 + t90 * (t25 * t46 - t39 * t85), t62 * t20 + t68 * t21, (-t21 * t53 + t32 * t57) * pkin(4) + t67, t67, (t20 * t56 - t8 * t52) * r_i_i_C(1) + (-t20 * t52 - t8 * t56) * r_i_i_C(2); 1, (t29 * t56 + t36 * t52) * r_i_i_C(1) + (-t29 * t52 + t36 * t56) * r_i_i_C(2) + t29 * pkin(5) + t37 * t45 - t36 * t59 + t90 * (t37 * t46 - t47 * t78) + (t58 * pkin(2) + t55 * t71) * t50, t62 * t30 + t68 * t31, (-t31 * t53 + t38 * t57) * pkin(4) + t65, t65, (-t19 * t52 + t30 * t56) * r_i_i_C(1) + (-t19 * t56 - t30 * t52) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end