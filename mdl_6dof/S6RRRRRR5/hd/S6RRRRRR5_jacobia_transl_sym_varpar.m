% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRRR5
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
% Datum: 2019-10-10 13:22
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
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
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
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
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
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
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
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
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (252->47), mult. (288->69), div. (0->0), fcn. (346->12), ass. (0->36)
	t28 = cos(pkin(6));
	t29 = sin(qJ(2));
	t32 = cos(qJ(1));
	t35 = t32 * t29;
	t30 = sin(qJ(1));
	t31 = cos(qJ(2));
	t36 = t30 * t31;
	t10 = t28 * t35 + t36;
	t26 = qJ(3) + qJ(4);
	t23 = qJ(5) + t26;
	t19 = sin(t23);
	t27 = sin(pkin(6));
	t38 = t27 * t32;
	t13 = t19 * t38;
	t20 = cos(t23);
	t45 = (-t10 * t19 - t20 * t38) * r_i_i_C(1) + (-t10 * t20 + t13) * r_i_i_C(2);
	t34 = t32 * t31;
	t37 = t30 * t29;
	t12 = -t28 * t37 + t34;
	t39 = t27 * t30;
	t5 = -t12 * t19 + t20 * t39;
	t6 = t12 * t20 + t19 * t39;
	t44 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t40 = t27 * t29;
	t43 = (-t19 * t40 + t28 * t20) * r_i_i_C(1) + (-t28 * t19 - t20 * t40) * r_i_i_C(2);
	t21 = sin(t26);
	t15 = pkin(4) * t21 + sin(qJ(3)) * pkin(3);
	t42 = pkin(8) + t15;
	t41 = r_i_i_C(3) + pkin(11) + pkin(10) + pkin(9);
	t22 = cos(t26);
	t16 = pkin(4) * t22 + cos(qJ(3)) * pkin(3);
	t14 = pkin(2) + t16;
	t33 = t20 * r_i_i_C(1) - t19 * r_i_i_C(2) + t14;
	t11 = t28 * t36 + t35;
	t9 = -t28 * t34 + t37;
	t1 = [-t30 * pkin(1) + t13 * r_i_i_C(1) - t41 * t9 - t33 * t10 + (r_i_i_C(2) * t20 + t42) * t38, -t33 * t11 + t41 * t12, -t12 * t15 + t16 * t39 + t44, (-t12 * t21 + t22 * t39) * pkin(4) + t44, t44, 0; t32 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t41 * t11 + t12 * t14 + t42 * t39, t41 * t10 - t33 * t9, -t10 * t15 - t16 * t38 + t45, (-t10 * t21 - t22 * t38) * pkin(4) + t45, t45, 0; 0, (t41 * t29 + t33 * t31) * t27, -t15 * t40 + t28 * t16 + t43, (-t21 * t40 + t22 * t28) * pkin(4) + t43, t43, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (563->65), mult. (639->99), div. (0->0), fcn. (794->14), ass. (0->46)
	t64 = pkin(12) + r_i_i_C(3);
	t40 = sin(qJ(6));
	t43 = cos(qJ(6));
	t68 = t43 * r_i_i_C(1) - t40 * r_i_i_C(2) + pkin(5);
	t41 = sin(qJ(2));
	t42 = sin(qJ(1));
	t44 = cos(qJ(2));
	t45 = cos(qJ(1));
	t55 = cos(pkin(6));
	t52 = t45 * t55;
	t21 = t41 * t52 + t42 * t44;
	t38 = qJ(3) + qJ(4);
	t35 = qJ(5) + t38;
	t31 = sin(t35);
	t32 = cos(t35);
	t39 = sin(pkin(6));
	t56 = t39 * t45;
	t10 = t21 * t32 - t31 * t56;
	t20 = t42 * t41 - t44 * t52;
	t67 = t10 * t40 - t20 * t43;
	t66 = -t10 * t43 - t20 * t40;
	t34 = cos(t38);
	t28 = pkin(4) * t34 + cos(qJ(3)) * pkin(3);
	t26 = pkin(2) + t28;
	t65 = t64 * t31 + t68 * t32 + t26;
	t59 = t39 * t41;
	t58 = t39 * t42;
	t57 = t39 * t44;
	t33 = sin(t38);
	t27 = pkin(4) * t33 + sin(qJ(3)) * pkin(3);
	t54 = t39 * (pkin(8) + t27);
	t53 = t42 * t55;
	t37 = -pkin(11) - pkin(10) - pkin(9);
	t50 = t40 * r_i_i_C(1) + t43 * r_i_i_C(2) - t37;
	t9 = -t21 * t31 - t32 * t56;
	t49 = t64 * t10 + t68 * t9;
	t23 = -t41 * t53 + t45 * t44;
	t13 = t23 * t31 - t32 * t58;
	t14 = t23 * t32 + t31 * t58;
	t48 = -t68 * t13 + t64 * t14;
	t19 = t55 * t31 + t32 * t59;
	t47 = t64 * t19 + t68 * (-t31 * t59 + t55 * t32);
	t22 = t45 * t41 + t44 * t53;
	t2 = t14 * t43 + t22 * t40;
	t1 = -t14 * t40 + t22 * t43;
	t3 = [-t42 * pkin(1) - t10 * pkin(5) + t66 * r_i_i_C(1) + t67 * r_i_i_C(2) + t20 * t37 - t21 * t26 + t45 * t54 + t64 * t9, -t22 * t65 + t50 * t23, -t23 * t27 + t28 * t58 + t48, (-t23 * t33 + t34 * t58) * pkin(4) + t48, t48, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t45 * pkin(1) + t14 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t64 * t13 - t22 * t37 + t23 * t26 + t42 * t54, -t20 * t65 + t50 * t21, -t21 * t27 - t28 * t56 + t49, (-t21 * t33 - t34 * t56) * pkin(4) + t49, t49, -t67 * r_i_i_C(1) + t66 * r_i_i_C(2); 0, (t50 * t41 + t65 * t44) * t39, -t27 * t59 + t55 * t28 + t47, (-t33 * t59 + t55 * t34) * pkin(4) + t47, t47, (-t19 * t40 - t43 * t57) * r_i_i_C(1) + (-t19 * t43 + t40 * t57) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end