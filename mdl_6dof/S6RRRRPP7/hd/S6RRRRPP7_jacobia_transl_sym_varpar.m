% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:20
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
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:20
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
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:20
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
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:20
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (150->46), mult. (381->80), div. (0->0), fcn. (482->10), ass. (0->34)
	t21 = sin(qJ(2));
	t22 = sin(qJ(1));
	t25 = cos(qJ(2));
	t26 = cos(qJ(1));
	t33 = cos(pkin(6));
	t31 = t26 * t33;
	t11 = t22 * t21 - t25 * t31;
	t19 = sin(qJ(4));
	t23 = cos(qJ(4));
	t12 = t21 * t31 + t22 * t25;
	t20 = sin(qJ(3));
	t24 = cos(qJ(3));
	t18 = sin(pkin(6));
	t34 = t18 * t26;
	t4 = t12 * t24 - t20 * t34;
	t43 = -t11 * t23 + t4 * t19;
	t42 = -t11 * t19 - t4 * t23;
	t30 = t23 * r_i_i_C(1) - t19 * r_i_i_C(2) + pkin(3);
	t40 = pkin(10) + r_i_i_C(3);
	t41 = t40 * t20 + t30 * t24 + pkin(2);
	t37 = t18 * t22;
	t36 = t18 * t24;
	t35 = t18 * t25;
	t32 = t22 * t33;
	t29 = t19 * r_i_i_C(1) + t23 * r_i_i_C(2) + pkin(9);
	t28 = -t12 * t20 - t24 * t34;
	t14 = -t21 * t32 + t26 * t25;
	t13 = t26 * t21 + t25 * t32;
	t10 = t33 * t20 + t21 * t36;
	t8 = t14 * t24 + t20 * t37;
	t7 = t14 * t20 - t22 * t36;
	t2 = t13 * t19 + t8 * t23;
	t1 = t13 * t23 - t8 * t19;
	t3 = [-t22 * pkin(1) - t12 * pkin(2) - t4 * pkin(3) + pkin(8) * t34 - t11 * pkin(9) + t42 * r_i_i_C(1) + t43 * r_i_i_C(2) + t40 * t28, -t13 * t41 + t29 * t14, -t30 * t7 + t40 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t26 * pkin(1) + t14 * pkin(2) + t8 * pkin(3) + pkin(8) * t37 + t13 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t40 * t7, -t11 * t41 + t29 * t12, t30 * t28 + t40 * t4, -t43 * r_i_i_C(1) + t42 * r_i_i_C(2), 0, 0; 0, (t29 * t21 + t41 * t25) * t18, t40 * t10 + t30 * (-t18 * t21 * t20 + t33 * t24), (-t10 * t19 - t23 * t35) * r_i_i_C(1) + (-t10 * t23 + t19 * t35) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:20
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (231->53), mult. (461->89), div. (0->0), fcn. (581->12), ass. (0->38)
	t26 = sin(qJ(3));
	t30 = cos(qJ(3));
	t29 = cos(qJ(4));
	t19 = t29 * pkin(4) + pkin(3);
	t22 = qJ(4) + pkin(11);
	t20 = sin(t22);
	t21 = cos(t22);
	t34 = t21 * r_i_i_C(1) - t20 * r_i_i_C(2) + t19;
	t44 = r_i_i_C(3) + qJ(5) + pkin(10);
	t45 = t44 * t26 + t34 * t30 + pkin(2);
	t43 = cos(qJ(1));
	t23 = sin(pkin(6));
	t28 = sin(qJ(1));
	t42 = t23 * t28;
	t41 = t23 * t30;
	t31 = cos(qJ(2));
	t40 = t23 * t31;
	t39 = cos(pkin(6));
	t25 = sin(qJ(4));
	t38 = t25 * pkin(4) + pkin(9);
	t37 = t23 * t43;
	t27 = sin(qJ(2));
	t35 = t39 * t43;
	t12 = t27 * t35 + t28 * t31;
	t4 = t12 * t30 - t26 * t37;
	t36 = t28 * t39;
	t3 = t12 * t26 + t30 * t37;
	t33 = t20 * r_i_i_C(1) + t21 * r_i_i_C(2) + t38;
	t14 = -t27 * t36 + t43 * t31;
	t13 = t43 * t27 + t31 * t36;
	t11 = t28 * t27 - t31 * t35;
	t10 = t39 * t26 + t27 * t41;
	t9 = t23 * t27 * t26 - t39 * t30;
	t8 = t14 * t30 + t26 * t42;
	t7 = t14 * t26 - t28 * t41;
	t2 = t13 * t20 + t8 * t21;
	t1 = t13 * t21 - t8 * t20;
	t5 = [-t28 * pkin(1) - t12 * pkin(2) + pkin(8) * t37 - t33 * t11 - t44 * t3 - t34 * t4, -t13 * t45 + t33 * t14, -t34 * t7 + t44 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + (t13 * t29 - t25 * t8) * pkin(4), t7, 0; t43 * pkin(1) + t14 * pkin(2) + pkin(8) * t42 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t38 * t13 + t8 * t19 + t44 * t7, -t11 * t45 + t33 * t12, -t34 * t3 + t44 * t4, (t11 * t21 - t4 * t20) * r_i_i_C(1) + (-t11 * t20 - t4 * t21) * r_i_i_C(2) + (t11 * t29 - t4 * t25) * pkin(4), t3, 0; 0, (t33 * t27 + t31 * t45) * t23, t44 * t10 - t34 * t9, (-t10 * t20 - t21 * t40) * r_i_i_C(1) + (-t10 * t21 + t20 * t40) * r_i_i_C(2) + (-t10 * t25 - t29 * t40) * pkin(4), t9, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:29:20
	% EndTime: 2019-10-10 12:29:20
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (367->67), mult. (693->110), div. (0->0), fcn. (886->12), ass. (0->45)
	t41 = sin(qJ(2));
	t42 = sin(qJ(1));
	t45 = cos(qJ(2));
	t53 = cos(pkin(6));
	t63 = cos(qJ(1));
	t48 = t53 * t63;
	t26 = t41 * t48 + t42 * t45;
	t40 = sin(qJ(3));
	t44 = cos(qJ(3));
	t37 = sin(pkin(6));
	t51 = t37 * t63;
	t14 = t26 * t44 - t40 * t51;
	t25 = t42 * t41 - t45 * t48;
	t36 = qJ(4) + pkin(11);
	t34 = sin(t36);
	t35 = cos(t36);
	t1 = t14 * t34 - t25 * t35;
	t67 = t14 * t35 + t25 * t34;
	t43 = cos(qJ(4));
	t33 = t43 * pkin(4) + pkin(3);
	t64 = r_i_i_C(2) + qJ(5) + pkin(10);
	t66 = t33 * t44 + t64 * t40 + pkin(2);
	t65 = pkin(5) + r_i_i_C(1);
	t60 = t34 * t44;
	t59 = t35 * t44;
	t58 = t37 * t42;
	t57 = t37 * t44;
	t56 = t37 * t45;
	t55 = t44 * t45;
	t54 = r_i_i_C(3) + qJ(6);
	t39 = sin(qJ(4));
	t52 = pkin(4) * t39 + pkin(9);
	t49 = t42 * t53;
	t13 = t26 * t40 + t44 * t51;
	t46 = -t54 * t34 - t65 * t35 - t33;
	t28 = -t41 * t49 + t63 * t45;
	t27 = t63 * t41 + t45 * t49;
	t24 = t53 * t40 + t41 * t57;
	t23 = t37 * t41 * t40 - t53 * t44;
	t18 = t28 * t44 + t40 * t58;
	t17 = t28 * t40 - t42 * t57;
	t11 = t24 * t34 + t35 * t56;
	t6 = t18 * t35 + t27 * t34;
	t5 = t18 * t34 - t27 * t35;
	t2 = [-t42 * pkin(1) - t26 * pkin(2) + pkin(8) * t51 - t54 * t1 - t64 * t13 - t14 * t33 - t52 * t25 - t65 * t67, t54 * (-t27 * t60 - t28 * t35) + t52 * t28 + t65 * (-t27 * t59 + t28 * t34) - t66 * t27, t46 * t17 + t64 * t18, t54 * t6 - t65 * t5 + (-t18 * t39 + t27 * t43) * pkin(4), t17, t5; t63 * pkin(1) + t28 * pkin(2) + pkin(8) * t58 + t64 * t17 + t18 * t33 + t52 * t27 + t54 * t5 + t65 * t6, t65 * (-t25 * t59 + t26 * t34) + t54 * (-t25 * t60 - t26 * t35) + t52 * t26 - t66 * t25, t46 * t13 + t64 * t14, t54 * t67 - t65 * t1 + (-t14 * t39 + t25 * t43) * pkin(4), t13, t1; 0, (t65 * (t34 * t41 + t35 * t55) + t54 * (t34 * t55 - t35 * t41) + t52 * t41 + t66 * t45) * t37, t46 * t23 + t64 * t24, t54 * (t24 * t35 - t34 * t56) - t65 * t11 + (-t24 * t39 - t43 * t56) * pkin(4), t23, t11;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end