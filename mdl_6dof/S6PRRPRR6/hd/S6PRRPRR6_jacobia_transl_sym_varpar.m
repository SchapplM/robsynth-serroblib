% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(12));
	t1 = sin(pkin(12));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
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
	t6 = sin(pkin(12));
	t9 = cos(pkin(12));
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
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (139->47), mult. (392->87), div. (0->0), fcn. (504->12), ass. (0->39)
	t24 = sin(pkin(7));
	t51 = t24 * pkin(9);
	t23 = sin(pkin(13));
	t50 = t23 * t24;
	t25 = sin(pkin(6));
	t49 = t24 * t25;
	t26 = cos(pkin(13));
	t48 = t24 * t26;
	t28 = cos(pkin(7));
	t29 = sin(qJ(3));
	t47 = t28 * t29;
	t31 = cos(qJ(3));
	t46 = t28 * t31;
	t30 = sin(qJ(2));
	t45 = t29 * t30;
	t32 = cos(qJ(2));
	t44 = t29 * t32;
	t43 = t30 * t31;
	t42 = t31 * t32;
	t41 = r_i_i_C(3) + qJ(4);
	t40 = cos(pkin(6));
	t39 = sin(pkin(12));
	t27 = cos(pkin(12));
	t38 = t27 * t49;
	t37 = t27 * t40;
	t36 = t40 * t24;
	t35 = t39 * t49;
	t34 = t40 * t39;
	t33 = t26 * r_i_i_C(1) - t23 * r_i_i_C(2) + pkin(3);
	t18 = t27 * t32 - t30 * t34;
	t17 = -t27 * t30 - t32 * t34;
	t16 = t30 * t37 + t39 * t32;
	t15 = -t39 * t30 + t32 * t37;
	t9 = -t31 * t36 + (-t28 * t42 + t45) * t25;
	t8 = t17 * t31 - t18 * t47;
	t6 = t15 * t31 - t16 * t47;
	t3 = -t17 * t46 + t18 * t29 - t31 * t35;
	t1 = -t15 * t46 + t16 * t29 + t31 * t38;
	t2 = [0, (t18 * t50 + t8 * t26) * r_i_i_C(1) + (t18 * t48 - t8 * t23) * r_i_i_C(2) + t8 * pkin(3) + t17 * pkin(2) + t18 * t51 + t41 * (t17 * t29 + t18 * t46), t41 * (t18 * t31 + (t17 * t28 + t35) * t29) - t33 * t3, t3, 0, 0; 0, (t16 * t50 + t6 * t26) * r_i_i_C(1) + (t16 * t48 - t6 * t23) * r_i_i_C(2) + t6 * pkin(3) + t15 * pkin(2) + t16 * t51 + t41 * (t15 * t29 + t16 * t46), t41 * (t16 * t31 + (t15 * t28 - t38) * t29) - t33 * t1, t1, 0, 0; 1, (t33 * (-t28 * t45 + t42) + t32 * pkin(2) + (r_i_i_C(1) * t23 + r_i_i_C(2) * t26 + pkin(9)) * t30 * t24 + t41 * (t28 * t43 + t44)) * t25, t41 * (t29 * t36 + (t28 * t44 + t43) * t25) - t33 * t9, t9, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (225->54), mult. (524->100), div. (0->0), fcn. (674->14), ass. (0->44)
	t58 = r_i_i_C(3) + pkin(10) + qJ(4);
	t31 = sin(pkin(7));
	t32 = sin(pkin(6));
	t57 = t31 * t32;
	t34 = cos(pkin(7));
	t56 = t32 * t34;
	t36 = sin(qJ(3));
	t55 = t34 * t36;
	t38 = cos(qJ(3));
	t54 = t34 * t38;
	t37 = sin(qJ(2));
	t53 = t36 * t37;
	t39 = cos(qJ(2));
	t52 = t36 * t39;
	t51 = t37 * t38;
	t50 = t38 * t39;
	t49 = cos(pkin(6));
	t48 = sin(pkin(12));
	t33 = cos(pkin(12));
	t47 = t33 * t57;
	t46 = t32 * t48;
	t45 = t33 * t49;
	t44 = t49 * t31;
	t43 = t31 * t46;
	t42 = t49 * t48;
	t29 = pkin(13) + qJ(5);
	t27 = sin(t29);
	t28 = cos(t29);
	t41 = t28 * r_i_i_C(1) - t27 * r_i_i_C(2) + cos(pkin(13)) * pkin(4) + pkin(3);
	t40 = (pkin(4) * sin(pkin(13)) + t27 * r_i_i_C(1) + t28 * r_i_i_C(2) + pkin(9)) * t31;
	t21 = t33 * t39 - t37 * t42;
	t20 = -t33 * t37 - t39 * t42;
	t19 = t37 * t45 + t48 * t39;
	t18 = -t48 * t37 + t39 * t45;
	t17 = t49 * t34 - t39 * t57;
	t12 = -t20 * t31 + t34 * t46;
	t11 = -t18 * t31 - t33 * t56;
	t10 = t36 * t44 + (t34 * t52 + t51) * t32;
	t9 = t32 * t53 - t38 * t44 - t50 * t56;
	t4 = t21 * t38 + (t20 * t34 + t43) * t36;
	t3 = -t20 * t54 + t21 * t36 - t38 * t43;
	t2 = t19 * t38 + (t18 * t34 - t47) * t36;
	t1 = -t18 * t54 + t19 * t36 + t38 * t47;
	t5 = [0, t20 * pkin(2) + t58 * (t20 * t36 + t21 * t54) + t41 * (t20 * t38 - t21 * t55) + t21 * t40, -t41 * t3 + t58 * t4, t3, (t12 * t28 - t4 * t27) * r_i_i_C(1) + (-t12 * t27 - t4 * t28) * r_i_i_C(2), 0; 0, t18 * pkin(2) + t58 * (t18 * t36 + t19 * t54) + t41 * (t18 * t38 - t19 * t55) + t19 * t40, -t41 * t1 + t58 * t2, t1, (t11 * t28 - t2 * t27) * r_i_i_C(1) + (-t11 * t27 - t2 * t28) * r_i_i_C(2), 0; 1, (t58 * (t34 * t51 + t52) + t41 * (-t34 * t53 + t50) + t39 * pkin(2) + t37 * t40) * t32, t58 * t10 - t41 * t9, t9, (-t10 * t27 + t17 * t28) * r_i_i_C(1) + (-t10 * t28 - t17 * t27) * r_i_i_C(2), 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (517->88), mult. (1188->159), div. (0->0), fcn. (1549->16), ass. (0->62)
	t41 = sin(pkin(7));
	t47 = sin(qJ(2));
	t49 = cos(qJ(2));
	t43 = cos(pkin(12));
	t69 = cos(pkin(6));
	t63 = t43 * t69;
	t67 = sin(pkin(12));
	t53 = t67 * t47 - t49 * t63;
	t68 = cos(pkin(7));
	t42 = sin(pkin(6));
	t72 = t42 * t43;
	t78 = t41 * t72 + t53 * t68;
	t57 = t69 * t67;
	t54 = t43 * t47 + t49 * t57;
	t64 = t42 * t67;
	t77 = -t41 * t64 + t54 * t68;
	t76 = r_i_i_C(3) + pkin(11);
	t75 = cos(qJ(3));
	t39 = pkin(13) + qJ(5);
	t37 = sin(t39);
	t74 = t37 * t41;
	t38 = cos(t39);
	t73 = t38 * t41;
	t71 = t42 * t47;
	t70 = t42 * t49;
	t66 = t41 * t71;
	t46 = sin(qJ(3));
	t62 = t46 * t68;
	t61 = t69 * t41;
	t59 = (pkin(4) * sin(pkin(13)) + pkin(9)) * t41;
	t58 = t68 * t75;
	t45 = sin(qJ(6));
	t48 = cos(qJ(6));
	t56 = t48 * r_i_i_C(1) - t45 * r_i_i_C(2) + pkin(5);
	t44 = -pkin(10) - qJ(4);
	t55 = t45 * r_i_i_C(1) + t48 * r_i_i_C(2) - t44;
	t36 = cos(pkin(13)) * pkin(4) + pkin(3);
	t52 = -t76 * t37 - t56 * t38 - t36;
	t31 = t43 * t49 - t47 * t57;
	t30 = t47 * t63 + t67 * t49;
	t29 = -t41 * t70 + t69 * t68;
	t28 = (-t47 * t62 + t75 * t49) * t42;
	t27 = (t46 * t49 + t47 * t58) * t42;
	t24 = t54 * t41 + t68 * t64;
	t23 = t53 * t41 - t68 * t72;
	t22 = t46 * t61 + (t75 * t47 + t49 * t62) * t42;
	t21 = t46 * t71 - t58 * t70 - t75 * t61;
	t20 = t28 * t38 + t37 * t66;
	t18 = -t31 * t62 - t54 * t75;
	t17 = t31 * t58 - t54 * t46;
	t16 = -t30 * t62 - t53 * t75;
	t15 = t30 * t58 - t53 * t46;
	t14 = t31 * t75 - t77 * t46;
	t13 = t31 * t46 + t77 * t75;
	t12 = t30 * t75 - t78 * t46;
	t11 = t30 * t46 + t78 * t75;
	t10 = t22 * t38 + t29 * t37;
	t8 = t18 * t38 + t31 * t74;
	t6 = t16 * t38 + t30 * t74;
	t4 = t14 * t38 + t24 * t37;
	t2 = t12 * t38 + t23 * t37;
	t1 = [0, (t17 * t45 + t8 * t48) * r_i_i_C(1) + (t17 * t48 - t8 * t45) * r_i_i_C(2) + t8 * pkin(5) + t18 * t36 - t17 * t44 - t54 * pkin(2) + t76 * (t18 * t37 - t31 * t73) + t31 * t59, t52 * t13 + t55 * t14, t13, t76 * t4 + t56 * (-t14 * t37 + t24 * t38), (t13 * t48 - t4 * t45) * r_i_i_C(1) + (-t13 * t45 - t4 * t48) * r_i_i_C(2); 0, (t15 * t45 + t6 * t48) * r_i_i_C(1) + (t15 * t48 - t6 * t45) * r_i_i_C(2) + t6 * pkin(5) + t16 * t36 - t15 * t44 - t53 * pkin(2) + t76 * (t16 * t37 - t30 * t73) + t30 * t59, t52 * t11 + t55 * t12, t11, t76 * t2 + t56 * (-t12 * t37 + t23 * t38), (t11 * t48 - t2 * t45) * r_i_i_C(1) + (-t11 * t45 - t2 * t48) * r_i_i_C(2); 1, (t20 * t48 + t27 * t45) * r_i_i_C(1) + (-t20 * t45 + t27 * t48) * r_i_i_C(2) + t20 * pkin(5) + t28 * t36 - t27 * t44 + t76 * (t28 * t37 - t38 * t66) + (t49 * pkin(2) + t47 * t59) * t42, t52 * t21 + t55 * t22, t21, t76 * t10 + t56 * (-t22 * t37 + t29 * t38), (-t10 * t45 + t21 * t48) * r_i_i_C(1) + (-t10 * t48 - t21 * t45) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end