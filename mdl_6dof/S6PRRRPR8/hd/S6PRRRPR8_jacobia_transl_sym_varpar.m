% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
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
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
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
	t21 = sin(pkin(12));
	t24 = cos(pkin(12));
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
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (284->61), mult. (791->116), div. (0->0), fcn. (1027->12), ass. (0->51)
	t66 = pkin(4) - r_i_i_C(2);
	t65 = pkin(10) + r_i_i_C(1);
	t34 = sin(pkin(7));
	t64 = pkin(9) * t34;
	t35 = sin(pkin(6));
	t63 = t34 * t35;
	t38 = cos(pkin(6));
	t62 = t34 * t38;
	t39 = sin(qJ(4));
	t61 = t34 * t39;
	t41 = sin(qJ(2));
	t60 = t34 * t41;
	t42 = cos(qJ(4));
	t59 = t34 * t42;
	t37 = cos(pkin(7));
	t58 = t35 * t37;
	t40 = sin(qJ(3));
	t57 = t37 * t40;
	t43 = cos(qJ(3));
	t56 = t37 * t43;
	t55 = t38 * t41;
	t44 = cos(qJ(2));
	t54 = t38 * t44;
	t53 = t40 * t41;
	t52 = t40 * t44;
	t51 = t41 * t43;
	t50 = t43 * t44;
	t49 = r_i_i_C(3) + qJ(5);
	t48 = t35 * t60;
	t33 = sin(pkin(12));
	t36 = cos(pkin(12));
	t28 = -t33 * t41 + t36 * t54;
	t47 = t28 * t37 - t36 * t63;
	t30 = -t33 * t54 - t36 * t41;
	t46 = t30 * t37 + t33 * t63;
	t45 = t49 * t39 + t66 * t42 + pkin(3);
	t31 = -t33 * t55 + t36 * t44;
	t29 = t33 * t44 + t36 * t55;
	t27 = t38 * t37 - t44 * t63;
	t26 = (-t37 * t53 + t50) * t35;
	t24 = -t30 * t34 + t33 * t58;
	t23 = -t28 * t34 - t36 * t58;
	t22 = t40 * t62 + (t37 * t52 + t51) * t35;
	t18 = t30 * t43 - t31 * t57;
	t16 = t28 * t43 - t29 * t57;
	t13 = t22 * t39 - t27 * t42;
	t12 = t31 * t43 + t46 * t40;
	t10 = t29 * t43 + t47 * t40;
	t3 = t12 * t39 - t24 * t42;
	t1 = t10 * t39 - t23 * t42;
	t2 = [0, t31 * t64 + t30 * pkin(2) + t18 * pkin(3) + t66 * (t18 * t42 + t31 * t61) + t49 * (t18 * t39 - t31 * t59) + t65 * (t30 * t40 + t31 * t56), t65 * t12 + t45 * (-t31 * t40 + t46 * t43), t49 * (t12 * t42 + t24 * t39) - t66 * t3, t3, 0; 0, t29 * t64 + t28 * pkin(2) + t16 * pkin(3) + t66 * (t16 * t42 + t29 * t61) + t49 * (t16 * t39 - t29 * t59) + t65 * (t28 * t40 + t29 * t56), t65 * t10 + t45 * (-t29 * t40 + t47 * t43), t49 * (t10 * t42 + t23 * t39) - t66 * t1, t1, 0; 1, t26 * pkin(3) + t66 * (t26 * t42 + t39 * t48) + t49 * (t26 * t39 - t42 * t48) + (pkin(2) * t44 + pkin(9) * t60 + t65 * (t37 * t51 + t52)) * t35, t65 * t22 + t45 * (t43 * t62 + (t37 * t50 - t53) * t35), t49 * (t22 * t42 + t27 * t39) - t66 * t13, t13, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.43s
	% Computational Cost: add. (493->73), mult. (1381->138), div. (0->0), fcn. (1805->14), ass. (0->53)
	t36 = sin(pkin(7));
	t42 = sin(qJ(2));
	t45 = cos(qJ(2));
	t38 = cos(pkin(12));
	t68 = cos(pkin(6));
	t61 = t38 * t68;
	t66 = sin(pkin(12));
	t52 = t66 * t42 - t45 * t61;
	t67 = cos(pkin(7));
	t37 = sin(pkin(6));
	t71 = t37 * t38;
	t77 = t36 * t71 + t52 * t67;
	t56 = t68 * t66;
	t53 = t38 * t42 + t45 * t56;
	t62 = t37 * t66;
	t76 = -t36 * t62 + t53 * t67;
	t75 = pkin(9) * t36;
	t74 = cos(qJ(3));
	t40 = sin(qJ(4));
	t73 = t36 * t40;
	t44 = cos(qJ(4));
	t72 = t36 * t44;
	t70 = t37 * t42;
	t69 = t37 * t45;
	t65 = r_i_i_C(3) + pkin(11) + pkin(4);
	t64 = t36 * t70;
	t41 = sin(qJ(3));
	t60 = t41 * t67;
	t59 = t68 * t36;
	t57 = t67 * t74;
	t39 = sin(qJ(6));
	t43 = cos(qJ(6));
	t55 = t39 * r_i_i_C(1) + t43 * r_i_i_C(2) + qJ(5);
	t54 = t43 * r_i_i_C(1) - t39 * r_i_i_C(2) + pkin(5) + pkin(10);
	t51 = -t36 * t69 + t68 * t67;
	t48 = -t55 * t40 - t65 * t44 - pkin(3);
	t47 = t52 * t36 - t67 * t71;
	t46 = t53 * t36 + t67 * t62;
	t31 = t38 * t45 - t42 * t56;
	t30 = t42 * t61 + t66 * t45;
	t28 = (-t42 * t60 + t74 * t45) * t37;
	t24 = t41 * t59 + (t74 * t42 + t45 * t60) * t37;
	t23 = t41 * t70 - t57 * t69 - t74 * t59;
	t18 = -t31 * t60 - t53 * t74;
	t16 = -t30 * t60 - t52 * t74;
	t13 = t24 * t40 - t51 * t44;
	t12 = t31 * t74 - t76 * t41;
	t11 = t31 * t41 + t76 * t74;
	t10 = t30 * t74 - t77 * t41;
	t9 = t30 * t41 + t77 * t74;
	t3 = t12 * t40 - t46 * t44;
	t1 = t10 * t40 - t47 * t44;
	t2 = [0, t18 * pkin(3) - t53 * pkin(2) + t31 * t75 + t55 * (t18 * t40 - t31 * t72) + t54 * (t31 * t57 - t53 * t41) + t65 * (t18 * t44 + t31 * t73), t48 * t11 + t54 * t12, t55 * (t12 * t44 + t46 * t40) - t65 * t3, t3, (-t11 * t39 + t3 * t43) * r_i_i_C(1) + (-t11 * t43 - t3 * t39) * r_i_i_C(2); 0, t16 * pkin(3) - t52 * pkin(2) + t30 * t75 + t55 * (t16 * t40 - t30 * t72) + t54 * (t30 * t57 - t52 * t41) + t65 * (t16 * t44 + t30 * t73), t54 * t10 + t48 * t9, t55 * (t10 * t44 + t47 * t40) - t65 * t1, t1, (t1 * t43 - t9 * t39) * r_i_i_C(1) + (-t1 * t39 - t9 * t43) * r_i_i_C(2); 1, t28 * pkin(3) + t55 * (t28 * t40 - t44 * t64) + t65 * (t28 * t44 + t40 * t64) + (t45 * pkin(2) + t42 * t75 + t54 * (t41 * t45 + t42 * t57)) * t37, t48 * t23 + t54 * t24, t55 * (t24 * t44 + t51 * t40) - t65 * t13, t13, (t13 * t43 - t23 * t39) * r_i_i_C(1) + (-t13 * t39 - t23 * t43) * r_i_i_C(2);];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end