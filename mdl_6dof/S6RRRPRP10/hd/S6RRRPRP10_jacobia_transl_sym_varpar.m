% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRP10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:49
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP10_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP10_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP10_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
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
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
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
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
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
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (132->41), mult. (334->69), div. (0->0), fcn. (423->10), ass. (0->29)
	t20 = sin(qJ(3));
	t23 = cos(qJ(3));
	t17 = sin(pkin(11));
	t19 = cos(pkin(11));
	t28 = r_i_i_C(1) * t19 - r_i_i_C(2) * t17 + pkin(3);
	t33 = r_i_i_C(3) + qJ(4);
	t37 = t33 * t20 + t28 * t23 + pkin(2);
	t18 = sin(pkin(6));
	t22 = sin(qJ(1));
	t36 = t18 * t22;
	t35 = t18 * t23;
	t25 = cos(qJ(1));
	t34 = t18 * t25;
	t32 = cos(pkin(6));
	t21 = sin(qJ(2));
	t24 = cos(qJ(2));
	t29 = t25 * t32;
	t10 = t21 * t29 + t22 * t24;
	t31 = t10 * t23 - t20 * t34;
	t30 = t22 * t32;
	t27 = t17 * r_i_i_C(1) + t19 * r_i_i_C(2) + pkin(9);
	t1 = t10 * t20 + t23 * t34;
	t12 = -t21 * t30 + t25 * t24;
	t11 = t25 * t21 + t24 * t30;
	t9 = t22 * t21 - t24 * t29;
	t7 = t18 * t21 * t20 - t32 * t23;
	t6 = t12 * t23 + t20 * t36;
	t5 = t12 * t20 - t22 * t35;
	t2 = [(-t9 * t17 - t19 * t31) * r_i_i_C(1) + (t17 * t31 - t9 * t19) * r_i_i_C(2) - t31 * pkin(3) - t10 * pkin(2) - t9 * pkin(9) - t22 * pkin(1) + pkin(8) * t34 - t33 * t1, -t11 * t37 + t27 * t12, -t28 * t5 + t33 * t6, t5, 0, 0; (t11 * t17 + t6 * t19) * r_i_i_C(1) + (t11 * t19 - t6 * t17) * r_i_i_C(2) + t6 * pkin(3) + t12 * pkin(2) + t11 * pkin(9) + t25 * pkin(1) + pkin(8) * t36 + t33 * t5, t27 * t10 - t37 * t9, -t28 * t1 + t33 * t31, t1, 0, 0; 0, (t27 * t21 + t37 * t24) * t18, t33 * (t32 * t20 + t21 * t35) - t28 * t7, t7, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (218->47), mult. (427->80), div. (0->0), fcn. (540->12), ass. (0->36)
	t26 = sin(qJ(3));
	t29 = cos(qJ(3));
	t19 = cos(pkin(11)) * pkin(4) + pkin(3);
	t22 = pkin(11) + qJ(5);
	t20 = sin(t22);
	t21 = cos(t22);
	t33 = t21 * r_i_i_C(1) - t20 * r_i_i_C(2) + t19;
	t43 = r_i_i_C(3) + pkin(10) + qJ(4);
	t44 = t43 * t26 + t33 * t29 + pkin(2);
	t42 = cos(qJ(1));
	t24 = sin(pkin(6));
	t28 = sin(qJ(1));
	t41 = t24 * t28;
	t40 = t24 * t29;
	t30 = cos(qJ(2));
	t39 = t24 * t30;
	t38 = cos(pkin(6));
	t37 = sin(pkin(11)) * pkin(4) + pkin(9);
	t36 = t24 * t42;
	t27 = sin(qJ(2));
	t34 = t38 * t42;
	t12 = t27 * t34 + t28 * t30;
	t4 = t12 * t29 - t26 * t36;
	t35 = t28 * t38;
	t3 = t12 * t26 + t29 * t36;
	t32 = t20 * r_i_i_C(1) + t21 * r_i_i_C(2) + t37;
	t14 = -t27 * t35 + t42 * t30;
	t13 = t42 * t27 + t30 * t35;
	t11 = t28 * t27 - t30 * t34;
	t10 = t38 * t26 + t27 * t40;
	t9 = t24 * t27 * t26 - t38 * t29;
	t8 = t14 * t29 + t26 * t41;
	t7 = t14 * t26 - t28 * t40;
	t2 = t13 * t20 + t8 * t21;
	t1 = t13 * t21 - t8 * t20;
	t5 = [-t28 * pkin(1) - t12 * pkin(2) + pkin(8) * t36 - t32 * t11 - t43 * t3 - t33 * t4, -t13 * t44 + t32 * t14, -t33 * t7 + t43 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t42 * pkin(1) + t14 * pkin(2) + pkin(8) * t41 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t37 * t13 + t8 * t19 + t43 * t7, -t11 * t44 + t32 * t12, -t33 * t3 + t43 * t4, t3, (t11 * t21 - t4 * t20) * r_i_i_C(1) + (-t11 * t20 - t4 * t21) * r_i_i_C(2), 0; 0, (t32 * t27 + t44 * t30) * t24, t43 * t10 - t33 * t9, t9, (-t10 * t20 - t21 * t39) * r_i_i_C(1) + (-t10 * t21 + t20 * t39) * r_i_i_C(2), 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:49:37
	% EndTime: 2019-10-10 11:49:37
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (354->61), mult. (659->101), div. (0->0), fcn. (845->12), ass. (0->43)
	t41 = sin(qJ(2));
	t42 = sin(qJ(1));
	t44 = cos(qJ(2));
	t52 = cos(pkin(6));
	t62 = cos(qJ(1));
	t47 = t52 * t62;
	t26 = t41 * t47 + t42 * t44;
	t40 = sin(qJ(3));
	t43 = cos(qJ(3));
	t38 = sin(pkin(6));
	t50 = t38 * t62;
	t14 = t26 * t43 - t40 * t50;
	t25 = t42 * t41 - t44 * t47;
	t36 = pkin(11) + qJ(5);
	t34 = sin(t36);
	t35 = cos(t36);
	t1 = t14 * t34 - t25 * t35;
	t66 = t14 * t35 + t25 * t34;
	t33 = cos(pkin(11)) * pkin(4) + pkin(3);
	t63 = r_i_i_C(2) + pkin(10) + qJ(4);
	t65 = t33 * t43 + t63 * t40 + pkin(2);
	t64 = pkin(5) + r_i_i_C(1);
	t59 = t34 * t43;
	t58 = t35 * t43;
	t57 = t38 * t42;
	t56 = t38 * t43;
	t55 = t38 * t44;
	t54 = t43 * t44;
	t53 = r_i_i_C(3) + qJ(6);
	t51 = pkin(4) * sin(pkin(11)) + pkin(9);
	t48 = t42 * t52;
	t13 = t26 * t40 + t43 * t50;
	t45 = -t53 * t34 - t64 * t35 - t33;
	t28 = -t41 * t48 + t62 * t44;
	t27 = t62 * t41 + t44 * t48;
	t24 = t52 * t40 + t41 * t56;
	t23 = t38 * t41 * t40 - t52 * t43;
	t18 = t28 * t43 + t40 * t57;
	t17 = t28 * t40 - t42 * t56;
	t11 = t24 * t34 + t35 * t55;
	t6 = t18 * t35 + t27 * t34;
	t5 = t18 * t34 - t27 * t35;
	t2 = [-t42 * pkin(1) - t26 * pkin(2) + pkin(8) * t50 - t53 * t1 - t63 * t13 - t14 * t33 - t51 * t25 - t64 * t66, t53 * (-t27 * t59 - t28 * t35) + t51 * t28 + t64 * (-t27 * t58 + t28 * t34) - t65 * t27, t45 * t17 + t63 * t18, t17, -t64 * t5 + t53 * t6, t5; t62 * pkin(1) + t28 * pkin(2) + pkin(8) * t57 + t63 * t17 + t18 * t33 + t51 * t27 + t53 * t5 + t64 * t6, t64 * (-t25 * t58 + t26 * t34) + t53 * (-t25 * t59 - t26 * t35) + t51 * t26 - t65 * t25, t45 * t13 + t63 * t14, t13, -t64 * t1 + t53 * t66, t1; 0, (t64 * (t34 * t41 + t35 * t54) + t53 * (t34 * t54 - t35 * t41) + t51 * t41 + t65 * t44) * t38, t45 * t23 + t63 * t24, t23, t53 * (t24 * t35 - t34 * t55) - t64 * t11, t11;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end