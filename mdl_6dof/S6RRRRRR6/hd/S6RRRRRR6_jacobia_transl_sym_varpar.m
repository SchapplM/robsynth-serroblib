% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRR6
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:24:46
	% EndTime: 2019-10-10 13:24:46
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:24:46
	% EndTime: 2019-10-10 13:24:46
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
	% StartTime: 2019-10-10 13:24:46
	% EndTime: 2019-10-10 13:24:46
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
	t11 = (pkin(8) + r_i_i_C(3)) * t5;
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
	% StartTime: 2019-10-10 13:24:46
	% EndTime: 2019-10-10 13:24:46
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (61->29), mult. (150->51), div. (0->0), fcn. (184->8), ass. (0->26)
	t27 = -r_i_i_C(3) - pkin(9);
	t13 = sin(qJ(1));
	t9 = sin(pkin(6));
	t26 = t13 * t9;
	t14 = cos(qJ(3));
	t25 = t14 * t9;
	t16 = cos(qJ(1));
	t24 = t16 * t9;
	t12 = sin(qJ(2));
	t23 = t12 * t13;
	t22 = t12 * t16;
	t15 = cos(qJ(2));
	t21 = t13 * t15;
	t20 = t15 * t16;
	t11 = sin(qJ(3));
	t10 = cos(pkin(6));
	t4 = t10 * t22 + t21;
	t19 = t11 * t24 - t14 * t4;
	t18 = t14 * r_i_i_C(1) - t11 * r_i_i_C(2) + pkin(2);
	t17 = t4 * t11 + t14 * t24;
	t6 = -t10 * t23 + t20;
	t5 = t10 * t21 + t22;
	t3 = -t10 * t20 + t23;
	t2 = t11 * t26 + t14 * t6;
	t1 = -t11 * t6 + t13 * t25;
	t7 = [-t13 * pkin(1) - t4 * pkin(2) + pkin(8) * t24 + t19 * r_i_i_C(1) + t17 * r_i_i_C(2) + t27 * t3, -t18 * t5 - t27 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; pkin(1) * t16 + t6 * pkin(2) + pkin(8) * t26 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t27 * t5, -t18 * t3 - t27 * t4, -t17 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, (-t27 * t12 + t18 * t15) * t9, (-t11 * t12 * t9 + t10 * t14) * r_i_i_C(1) + (-t10 * t11 - t12 * t25) * r_i_i_C(2), 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:24:46
	% EndTime: 2019-10-10 13:24:46
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (132->38), mult. (218->61), div. (0->0), fcn. (264->10), ass. (0->33)
	t20 = cos(pkin(6));
	t22 = sin(qJ(2));
	t26 = cos(qJ(1));
	t31 = t26 * t22;
	t23 = sin(qJ(1));
	t25 = cos(qJ(2));
	t32 = t23 * t25;
	t10 = t20 * t31 + t32;
	t18 = qJ(3) + qJ(4);
	t16 = sin(t18);
	t19 = sin(pkin(6));
	t34 = t19 * t26;
	t13 = t16 * t34;
	t17 = cos(t18);
	t40 = (-t10 * t16 - t17 * t34) * r_i_i_C(1) + (-t10 * t17 + t13) * r_i_i_C(2);
	t30 = t26 * t25;
	t33 = t23 * t22;
	t12 = -t20 * t33 + t30;
	t35 = t19 * t23;
	t5 = -t12 * t16 + t17 * t35;
	t6 = t12 * t17 + t16 * t35;
	t39 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t36 = t19 * t22;
	t38 = (-t16 * t36 + t20 * t17) * r_i_i_C(1) + (-t20 * t16 - t17 * t36) * r_i_i_C(2);
	t37 = r_i_i_C(3) + pkin(10) + pkin(9);
	t21 = sin(qJ(3));
	t29 = pkin(3) * t21 + pkin(8);
	t24 = cos(qJ(3));
	t15 = t24 * pkin(3) + pkin(2);
	t28 = t17 * r_i_i_C(1) - t16 * r_i_i_C(2) + t15;
	t11 = t20 * t32 + t31;
	t9 = -t20 * t30 + t33;
	t1 = [-t23 * pkin(1) + t13 * r_i_i_C(1) - t37 * t9 - t28 * t10 + (r_i_i_C(2) * t17 + t29) * t34, -t11 * t28 + t12 * t37, (-t12 * t21 + t24 * t35) * pkin(3) + t39, t39, 0, 0; t26 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t11 * t37 + t12 * t15 + t29 * t35, t10 * t37 - t28 * t9, (-t10 * t21 - t24 * t34) * pkin(3) + t40, t40, 0, 0; 0, (t22 * t37 + t25 * t28) * t19, (t20 * t24 - t21 * t36) * pkin(3) + t38, t38, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:24:46
	% EndTime: 2019-10-10 13:24:47
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (314->56), mult. (509->91), div. (0->0), fcn. (637->12), ass. (0->43)
	t59 = pkin(11) + r_i_i_C(3);
	t32 = sin(qJ(5));
	t36 = cos(qJ(5));
	t63 = t36 * r_i_i_C(1) - t32 * r_i_i_C(2) + pkin(4);
	t34 = sin(qJ(2));
	t35 = sin(qJ(1));
	t38 = cos(qJ(2));
	t39 = cos(qJ(1));
	t50 = cos(pkin(6));
	t48 = t39 * t50;
	t21 = t34 * t48 + t35 * t38;
	t30 = qJ(3) + qJ(4);
	t28 = sin(t30);
	t29 = cos(t30);
	t31 = sin(pkin(6));
	t51 = t31 * t39;
	t10 = t21 * t29 - t28 * t51;
	t20 = t35 * t34 - t38 * t48;
	t62 = t10 * t32 - t20 * t36;
	t61 = -t10 * t36 - t20 * t32;
	t37 = cos(qJ(3));
	t27 = t37 * pkin(3) + pkin(2);
	t60 = t59 * t28 + t63 * t29 + t27;
	t54 = t31 * t34;
	t53 = t31 * t35;
	t52 = t31 * t38;
	t49 = t35 * t50;
	t33 = sin(qJ(3));
	t47 = t31 * (pkin(3) * t33 + pkin(8));
	t40 = -pkin(10) - pkin(9);
	t45 = t32 * r_i_i_C(1) + t36 * r_i_i_C(2) - t40;
	t9 = -t21 * t28 - t29 * t51;
	t44 = t59 * t10 + t63 * t9;
	t23 = -t34 * t49 + t39 * t38;
	t13 = t23 * t28 - t29 * t53;
	t14 = t23 * t29 + t28 * t53;
	t43 = -t63 * t13 + t59 * t14;
	t19 = t50 * t28 + t29 * t54;
	t42 = t59 * t19 + t63 * (-t28 * t54 + t50 * t29);
	t22 = t39 * t34 + t38 * t49;
	t2 = t14 * t36 + t22 * t32;
	t1 = -t14 * t32 + t22 * t36;
	t3 = [-t35 * pkin(1) - t10 * pkin(4) + t61 * r_i_i_C(1) + t62 * r_i_i_C(2) + t20 * t40 - t21 * t27 + t39 * t47 + t59 * t9, -t22 * t60 + t45 * t23, (-t23 * t33 + t37 * t53) * pkin(3) + t43, t43, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t39 * pkin(1) + t14 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t59 * t13 - t22 * t40 + t23 * t27 + t35 * t47, -t20 * t60 + t45 * t21, (-t21 * t33 - t37 * t51) * pkin(3) + t44, t44, -t62 * r_i_i_C(1) + t61 * r_i_i_C(2), 0; 0, (t45 * t34 + t60 * t38) * t31, (-t33 * t54 + t50 * t37) * pkin(3) + t42, t42, (-t19 * t32 - t36 * t52) * r_i_i_C(1) + (-t19 * t36 + t32 * t52) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:24:46
	% EndTime: 2019-10-10 13:24:47
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (455->63), mult. (639->100), div. (0->0), fcn. (798->14), ass. (0->48)
	t71 = r_i_i_C(3) + pkin(12) + pkin(11);
	t49 = cos(qJ(5));
	t36 = t49 * pkin(5) + pkin(4);
	t42 = qJ(5) + qJ(6);
	t38 = sin(t42);
	t40 = cos(t42);
	t60 = r_i_i_C(1) * t40 - r_i_i_C(2) * t38 + t36;
	t50 = cos(qJ(3));
	t37 = t50 * pkin(3) + pkin(2);
	t43 = qJ(3) + qJ(4);
	t39 = sin(t43);
	t41 = cos(t43);
	t76 = t71 * t39 + t60 * t41 + t37;
	t47 = sin(qJ(2));
	t48 = sin(qJ(1));
	t51 = cos(qJ(2));
	t52 = cos(qJ(1));
	t65 = cos(pkin(6));
	t62 = t52 * t65;
	t30 = t47 * t62 + t48 * t51;
	t44 = sin(pkin(6));
	t67 = t44 * t52;
	t18 = t30 * t41 - t39 * t67;
	t29 = t48 * t47 - t51 * t62;
	t75 = (-t18 * t38 + t29 * t40) * r_i_i_C(1) + (-t18 * t40 - t29 * t38) * r_i_i_C(2);
	t63 = t48 * t65;
	t32 = -t47 * t63 + t52 * t51;
	t69 = t44 * t48;
	t22 = t32 * t41 + t39 * t69;
	t31 = t52 * t47 + t51 * t63;
	t5 = -t22 * t38 + t31 * t40;
	t6 = t22 * t40 + t31 * t38;
	t74 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t70 = t44 * t47;
	t68 = t44 * t51;
	t28 = t65 * t39 + t41 * t70;
	t66 = (-t28 * t38 - t40 * t68) * r_i_i_C(1) + (-t28 * t40 + t38 * t68) * r_i_i_C(2);
	t45 = sin(qJ(5));
	t64 = t45 * pkin(5) + pkin(9) + pkin(10);
	t46 = sin(qJ(3));
	t61 = t44 * (pkin(3) * t46 + pkin(8));
	t17 = -t30 * t39 - t41 * t67;
	t59 = t60 * t17 + t71 * t18;
	t21 = t32 * t39 - t41 * t69;
	t58 = -t60 * t21 + t71 * t22;
	t57 = t71 * t28 + t60 * (-t39 * t70 + t65 * t41);
	t56 = t38 * r_i_i_C(1) + t40 * r_i_i_C(2) + t64;
	t1 = [-t48 * pkin(1) + t71 * t17 - t60 * t18 - t56 * t29 - t30 * t37 + t52 * t61, -t31 * t76 + t56 * t32, (-t32 * t46 + t50 * t69) * pkin(3) + t58, t58, (-t22 * t45 + t31 * t49) * pkin(5) + t74, t74; t52 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t71 * t21 + t22 * t36 + t64 * t31 + t32 * t37 + t48 * t61, -t29 * t76 + t56 * t30, (-t30 * t46 - t50 * t67) * pkin(3) + t59, t59, (-t18 * t45 + t29 * t49) * pkin(5) + t75, t75; 0, (t56 * t47 + t76 * t51) * t44, (-t46 * t70 + t65 * t50) * pkin(3) + t57, t57, (-t28 * t45 - t49 * t68) * pkin(5) + t66, t66;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end