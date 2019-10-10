% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:43
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR12_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR12_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR12_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
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
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (13->12), mult. (28->19), div. (0->0), fcn. (35->6), ass. (0->12)
	t1 = sin(pkin(12));
	t5 = sin(qJ(1));
	t11 = t5 * t1;
	t3 = cos(pkin(12));
	t10 = t5 * t3;
	t6 = cos(qJ(1));
	t9 = t6 * t1;
	t8 = t6 * t3;
	t2 = sin(pkin(6));
	t7 = t2 * (r_i_i_C(3) + qJ(2));
	t4 = cos(pkin(6));
	t12 = [(-t4 * t9 - t10) * r_i_i_C(1) + (-t4 * t8 + t11) * r_i_i_C(2) - t5 * pkin(1) + t6 * t7, t5 * t2, 0, 0, 0, 0; (-t4 * t11 + t8) * r_i_i_C(1) + (-t4 * t10 - t9) * r_i_i_C(2) + t6 * pkin(1) + t5 * t7, -t6 * t2, 0, 0, 0, 0; 0, t4, 0, 0, 0, 0;];
	Ja_transl = t12;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (62->32), mult. (166->58), div. (0->0), fcn. (213->10), ass. (0->29)
	t34 = r_i_i_C(3) + pkin(9);
	t13 = cos(pkin(7));
	t15 = sin(qJ(3));
	t17 = cos(qJ(3));
	t10 = sin(pkin(7));
	t11 = sin(pkin(6));
	t18 = cos(qJ(1));
	t26 = t18 * t11;
	t24 = t10 * t26;
	t12 = cos(pkin(12));
	t14 = cos(pkin(6));
	t29 = t14 * t18;
	t16 = sin(qJ(1));
	t9 = sin(pkin(12));
	t32 = t16 * t9;
	t3 = -t12 * t29 + t32;
	t27 = t16 * t12;
	t4 = t9 * t29 + t27;
	t33 = (t13 * t3 + t24) * t17 + t4 * t15;
	t30 = t13 * t15;
	t28 = t16 * t11;
	t25 = t11 * qJ(2);
	t5 = -t14 * t27 - t18 * t9;
	t22 = t10 * t28 + t13 * t5;
	t19 = t15 * t24 - t17 * t4 + t3 * t30;
	t6 = t12 * t18 - t14 * t32;
	t2 = t22 * t15 + t17 * t6;
	t1 = -t15 * t6 + t22 * t17;
	t7 = [t19 * r_i_i_C(1) + t33 * r_i_i_C(2) - t4 * pkin(2) - t16 * pkin(1) + t18 * t25 + t34 * (-t3 * t10 + t13 * t26), t28, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; t18 * pkin(1) + t6 * pkin(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t25 + t34 * (-t10 * t5 + t13 * t28), -t26, -t33 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, t14, (t17 * r_i_i_C(1) - t15 * r_i_i_C(2)) * t14 * t10 + ((t12 * t13 * t17 - t15 * t9) * r_i_i_C(1) + (-t12 * t30 - t17 * t9) * r_i_i_C(2)) * t11, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (176->48), mult. (480->84), div. (0->0), fcn. (625->12), ass. (0->42)
	t27 = cos(pkin(6));
	t25 = cos(pkin(12));
	t33 = cos(qJ(1));
	t41 = t33 * t25;
	t22 = sin(pkin(12));
	t30 = sin(qJ(1));
	t46 = t30 * t22;
	t16 = -t27 * t41 + t46;
	t23 = sin(pkin(7));
	t26 = cos(pkin(7));
	t24 = sin(pkin(6));
	t42 = t33 * t24;
	t11 = -t16 * t23 + t26 * t42;
	t28 = sin(qJ(4));
	t31 = cos(qJ(4));
	t43 = t33 * t22;
	t44 = t30 * t25;
	t17 = t27 * t43 + t44;
	t29 = sin(qJ(3));
	t32 = cos(qJ(3));
	t39 = t23 * t42;
	t47 = t26 * t29;
	t6 = t16 * t47 - t17 * t32 + t29 * t39;
	t52 = t11 * t31 - t6 * t28;
	t51 = t11 * t28 + t6 * t31;
	t36 = t27 * t44 + t43;
	t45 = t30 * t24;
	t50 = -t23 * t45 + t36 * t26;
	t49 = r_i_i_C(3) + pkin(10);
	t48 = t23 * t27;
	t40 = t24 * qJ(2);
	t37 = t31 * r_i_i_C(1) - t28 * r_i_i_C(2) + pkin(3);
	t34 = -t17 * t29 + (-t16 * t26 - t39) * t32;
	t13 = t23 * t36 + t26 * t45;
	t18 = -t27 * t46 + t41;
	t15 = -t24 * t25 * t23 + t27 * t26;
	t10 = t29 * t48 + (t22 * t32 + t25 * t47) * t24;
	t8 = t18 * t32 - t50 * t29;
	t7 = t18 * t29 + t50 * t32;
	t2 = t13 * t28 + t8 * t31;
	t1 = t13 * t31 - t8 * t28;
	t3 = [-t30 * pkin(1) - t17 * pkin(2) + t6 * pkin(3) + t11 * pkin(9) + t51 * r_i_i_C(1) + t52 * r_i_i_C(2) + t33 * t40 + t49 * t34, t45, -t37 * t7 + t49 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t33 * pkin(1) + t18 * pkin(2) + t8 * pkin(3) + t13 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t30 * t40 + t49 * t7, -t42, t37 * t34 - t49 * t6, -t52 * r_i_i_C(1) + t51 * r_i_i_C(2), 0, 0; 0, t27, t49 * t10 + t37 * (t32 * t48 + (t25 * t26 * t32 - t22 * t29) * t24), (-t10 * t28 + t15 * t31) * r_i_i_C(1) + (-t10 * t31 - t15 * t28) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (294->51), mult. (800->84), div. (0->0), fcn. (1050->12), ass. (0->45)
	t33 = cos(pkin(6));
	t31 = cos(pkin(12));
	t39 = cos(qJ(1));
	t49 = t39 * t31;
	t28 = sin(pkin(12));
	t36 = sin(qJ(1));
	t54 = t36 * t28;
	t22 = -t33 * t49 + t54;
	t51 = t39 * t28;
	t52 = t36 * t31;
	t23 = t33 * t51 + t52;
	t35 = sin(qJ(3));
	t38 = cos(qJ(3));
	t29 = sin(pkin(7));
	t30 = sin(pkin(6));
	t50 = t39 * t30;
	t46 = t29 * t50;
	t32 = cos(pkin(7));
	t55 = t32 * t35;
	t12 = t22 * t55 - t23 * t38 + t35 * t46;
	t18 = -t22 * t29 + t32 * t50;
	t34 = sin(qJ(4));
	t37 = cos(qJ(4));
	t61 = t12 * t37 + t18 * t34;
	t60 = t12 * t34 - t18 * t37;
	t44 = t33 * t52 + t51;
	t53 = t36 * t30;
	t59 = -t29 * t53 + t44 * t32;
	t48 = r_i_i_C(3) + qJ(5);
	t57 = -r_i_i_C(2) + pkin(4);
	t42 = t48 * t34 + t57 * t37 + pkin(3);
	t58 = r_i_i_C(1) + pkin(10);
	t56 = t29 * t33;
	t47 = t30 * qJ(2);
	t41 = -t23 * t35 + (-t22 * t32 - t46) * t38;
	t40 = t44 * t29 + t32 * t53;
	t24 = -t33 * t54 + t49;
	t21 = -t30 * t31 * t29 + t33 * t32;
	t17 = t35 * t56 + (t28 * t38 + t31 * t55) * t30;
	t14 = t24 * t38 - t59 * t35;
	t13 = t24 * t35 + t38 * t59;
	t7 = t17 * t34 - t21 * t37;
	t6 = t14 * t37 + t40 * t34;
	t5 = t14 * t34 - t40 * t37;
	t1 = [-t36 * pkin(1) - t23 * pkin(2) + t12 * pkin(3) + t18 * pkin(9) + t39 * t47 + t58 * t41 + t48 * t60 + t57 * t61, t53, -t13 * t42 + t58 * t14, t48 * t6 - t57 * t5, t5, 0; t39 * pkin(1) + t24 * pkin(2) + t14 * pkin(3) + t40 * pkin(9) + t58 * t13 + t36 * t47 + t48 * t5 + t57 * t6, -t50, -t12 * t58 + t42 * t41, -t48 * t61 + t57 * t60, -t60, 0; 0, t33, t58 * t17 + t42 * (t38 * t56 + (t31 * t32 * t38 - t28 * t35) * t30), t48 * (t17 * t37 + t21 * t34) - t57 * t7, t7, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (505->64), mult. (1388->108), div. (0->0), fcn. (1830->14), ass. (0->50)
	t38 = cos(qJ(1));
	t65 = cos(pkin(12));
	t67 = cos(pkin(6));
	t55 = t67 * t65;
	t63 = sin(pkin(12));
	t69 = sin(qJ(1));
	t48 = -t38 * t55 + t69 * t63;
	t32 = sin(pkin(6));
	t64 = sin(pkin(7));
	t59 = t32 * t64;
	t66 = cos(pkin(7));
	t78 = t38 * t59 + t48 * t66;
	t53 = t67 * t63;
	t25 = t38 * t53 + t69 * t65;
	t35 = sin(qJ(3));
	t70 = cos(qJ(3));
	t14 = -t25 * t70 + t78 * t35;
	t34 = sin(qJ(4));
	t37 = cos(qJ(4));
	t60 = t32 * t66;
	t40 = -t38 * t60 + t48 * t64;
	t77 = t14 * t34 + t40 * t37;
	t76 = t14 * t37 - t40 * t34;
	t11 = t25 * t35 + t78 * t70;
	t44 = t38 * t63 + t69 * t55;
	t72 = t44 * t66 - t69 * t59;
	t71 = pkin(5) + pkin(10);
	t68 = t38 * t32;
	t62 = r_i_i_C(3) + pkin(11) + pkin(4);
	t61 = t69 * t32;
	t54 = t67 * t64;
	t52 = t66 * t65;
	t33 = sin(qJ(6));
	t36 = cos(qJ(6));
	t51 = t33 * r_i_i_C(1) + t36 * r_i_i_C(2) + qJ(5);
	t49 = t36 * r_i_i_C(1) - t33 * r_i_i_C(2) + t71;
	t43 = -t65 * t59 + t67 * t66;
	t42 = -t51 * t34 - t62 * t37 - pkin(3);
	t39 = t44 * t64 + t69 * t60;
	t26 = t38 * t65 - t69 * t53;
	t20 = t35 * t54 + (t35 * t52 + t63 * t70) * t32;
	t19 = -t70 * t54 + (t35 * t63 - t52 * t70) * t32;
	t16 = t26 * t70 - t72 * t35;
	t15 = t26 * t35 + t72 * t70;
	t9 = t20 * t34 - t43 * t37;
	t8 = t16 * t37 + t39 * t34;
	t7 = t16 * t34 - t39 * t37;
	t2 = t15 * t36 + t7 * t33;
	t1 = -t15 * t33 + t7 * t36;
	t3 = [-t69 * pkin(1) - t25 * pkin(2) + t14 * pkin(3) - t40 * pkin(9) + qJ(2) * t68 - t49 * t11 + t51 * t77 + t62 * t76, t61, t42 * t15 + t49 * t16, t51 * t8 - t62 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t38 * pkin(1) + t26 * pkin(2) + t16 * pkin(3) + t39 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t61 + t7 * qJ(5) + t71 * t15 + t62 * t8, -t68, t42 * t11 - t14 * t49, -t51 * t76 + t62 * t77, -t77, (-t11 * t33 - t36 * t77) * r_i_i_C(1) + (-t11 * t36 + t33 * t77) * r_i_i_C(2); 0, t67, t42 * t19 + t49 * t20, -t62 * t9 + t51 * (t20 * t37 + t43 * t34), t9, (-t19 * t33 + t9 * t36) * r_i_i_C(1) + (-t19 * t36 - t9 * t33) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end