% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR12_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR12_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR12_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
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
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
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
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (40->17), mult. (89->25), div. (0->0), fcn. (110->6), ass. (0->18)
	t19 = pkin(2) - r_i_i_C(2);
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t18 = t10 * t9;
	t12 = cos(qJ(1));
	t17 = t12 * t9;
	t11 = cos(qJ(2));
	t16 = t10 * t11;
	t15 = t11 * t12;
	t14 = r_i_i_C(3) + qJ(3);
	t7 = sin(pkin(6));
	t13 = (pkin(8) + r_i_i_C(1)) * t7;
	t8 = cos(pkin(6));
	t4 = -t8 * t18 + t15;
	t3 = t8 * t16 + t17;
	t2 = t8 * t17 + t16;
	t1 = -t8 * t15 + t18;
	t5 = [-t10 * pkin(1) - t14 * t1 + t12 * t13 - t19 * t2, t14 * t4 - t19 * t3, t3, 0, 0, 0; t12 * pkin(1) + t10 * t13 + t14 * t3 + t19 * t4, -t19 * t1 + t14 * t2, t1, 0, 0, 0; 0, (t19 * t11 + t14 * t9) * t7, -t7 * t11, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (75->31), mult. (179->50), div. (0->0), fcn. (222->8), ass. (0->26)
	t27 = pkin(3) + pkin(8);
	t10 = sin(pkin(6));
	t14 = sin(qJ(1));
	t26 = t10 * t14;
	t16 = cos(qJ(2));
	t25 = t10 * t16;
	t17 = cos(qJ(1));
	t24 = t10 * t17;
	t13 = sin(qJ(2));
	t23 = t13 * t14;
	t22 = t13 * t17;
	t21 = t14 * t16;
	t20 = t16 * t17;
	t19 = -r_i_i_C(3) - pkin(9) - pkin(2);
	t12 = sin(qJ(4));
	t15 = cos(qJ(4));
	t18 = t12 * r_i_i_C(1) + t15 * r_i_i_C(2) + qJ(3);
	t11 = cos(pkin(6));
	t8 = t15 * t24;
	t6 = -t11 * t23 + t20;
	t5 = t11 * t21 + t22;
	t4 = t11 * t22 + t21;
	t3 = -t11 * t20 + t23;
	t2 = t12 * t5 + t15 * t26;
	t1 = -t12 * t26 + t15 * t5;
	t7 = [-t14 * pkin(1) + t8 * r_i_i_C(1) - t18 * t3 + (-t12 * r_i_i_C(2) + t27) * t24 + t19 * t4, t18 * t6 + t19 * t5, t5, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0; pkin(1) * t17 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * qJ(3) - t19 * t6 + t27 * t26, t18 * t4 + t19 * t3, t3, (t12 * t24 + t3 * t15) * r_i_i_C(1) + (-t12 * t3 + t8) * r_i_i_C(2), 0, 0; 0, (t18 * t13 - t19 * t16) * t10, -t25, (-t11 * t12 - t15 * t25) * r_i_i_C(1) + (-t11 * t15 + t12 * t25) * r_i_i_C(2), 0, 0;];
	Ja_transl = t7;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (127->40), mult. (232->62), div. (0->0), fcn. (286->10), ass. (0->30)
	t20 = cos(qJ(4));
	t33 = t20 * pkin(4) + pkin(3) + pkin(8);
	t14 = sin(pkin(6));
	t19 = sin(qJ(1));
	t32 = t14 * t19;
	t21 = cos(qJ(2));
	t31 = t14 * t21;
	t22 = cos(qJ(1));
	t30 = t14 * t22;
	t18 = sin(qJ(2));
	t29 = t19 * t18;
	t28 = t19 * t21;
	t27 = t22 * t18;
	t26 = t22 * t21;
	t25 = -r_i_i_C(3) - qJ(5) - pkin(9) - pkin(2);
	t17 = sin(qJ(4));
	t24 = t17 * pkin(4) + qJ(3);
	t13 = qJ(4) + pkin(11);
	t11 = sin(t13);
	t12 = cos(t13);
	t23 = t11 * r_i_i_C(1) + t12 * r_i_i_C(2) + t24;
	t15 = cos(pkin(6));
	t7 = t12 * t30;
	t6 = -t15 * t29 + t26;
	t5 = t15 * t28 + t27;
	t4 = t15 * t27 + t28;
	t3 = -t15 * t26 + t29;
	t2 = t5 * t11 + t12 * t32;
	t1 = -t11 * t32 + t5 * t12;
	t8 = [-t19 * pkin(1) + t7 * r_i_i_C(1) - t23 * t3 + (-t11 * r_i_i_C(2) + t33) * t30 + t25 * t4, t23 * t6 + t25 * t5, t5, t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + (-t17 * t32 + t20 * t5) * pkin(4), t6, 0; t22 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t24 * t5 - t25 * t6 + t33 * t32, t23 * t4 + t25 * t3, t3, (t11 * t30 + t3 * t12) * r_i_i_C(1) + (-t3 * t11 + t7) * r_i_i_C(2) + (t17 * t30 + t3 * t20) * pkin(4), t4, 0; 0, (t23 * t18 - t25 * t21) * t14, -t31, (-t15 * t11 - t12 * t31) * r_i_i_C(1) + (t11 * t31 - t15 * t12) * r_i_i_C(2) + (-t15 * t17 - t20 * t31) * pkin(4), t14 * t18, 0;];
	Ja_transl = t8;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (267->56), mult. (463->89), div. (0->0), fcn. (584->12), ass. (0->38)
	t47 = pkin(10) + r_i_i_C(3);
	t46 = pkin(2) + qJ(5) + pkin(9);
	t23 = sin(pkin(6));
	t27 = sin(qJ(2));
	t45 = t23 * t27;
	t28 = sin(qJ(1));
	t44 = t23 * t28;
	t31 = cos(qJ(2));
	t43 = t23 * t31;
	t32 = cos(qJ(1));
	t42 = t23 * t32;
	t41 = cos(pkin(6));
	t30 = cos(qJ(4));
	t40 = t23 * (t30 * pkin(4) + pkin(3) + pkin(8));
	t26 = sin(qJ(4));
	t39 = t26 * pkin(4) + qJ(3);
	t38 = t28 * t41;
	t37 = t32 * t41;
	t25 = sin(qJ(6));
	t29 = cos(qJ(6));
	t36 = t29 * r_i_i_C(1) - t25 * r_i_i_C(2) + pkin(5);
	t13 = t28 * t27 - t31 * t37;
	t22 = qJ(4) + pkin(11);
	t20 = sin(t22);
	t21 = cos(t22);
	t7 = -t13 * t20 + t21 * t42;
	t35 = t13 * t21 + t20 * t42;
	t34 = t25 * r_i_i_C(1) + t29 * r_i_i_C(2) + t46;
	t33 = t36 * t20 - t47 * t21 + t39;
	t16 = -t27 * t38 + t32 * t31;
	t15 = t32 * t27 + t31 * t38;
	t14 = t27 * t37 + t28 * t31;
	t12 = -t20 * t43 + t41 * t21;
	t4 = t15 * t20 + t21 * t44;
	t3 = -t15 * t21 + t20 * t44;
	t2 = t16 * t25 + t4 * t29;
	t1 = t16 * t29 - t4 * t25;
	t5 = [-t28 * pkin(1) - t39 * t13 - t34 * t14 + t32 * t40 + t47 * t35 + t36 * t7, -t34 * t15 + t33 * t16, t15, t47 * t4 + (t15 * t30 - t26 * t44) * pkin(4) - t36 * t3, t16, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t32 * pkin(1) + t4 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t39 * t15 + t46 * t16 + t28 * t40 + t47 * t3, -t34 * t13 + t33 * t14, t13, -t47 * t7 + (t13 * t30 + t26 * t42) * pkin(4) + t36 * t35, t14, (t14 * t29 + t7 * t25) * r_i_i_C(1) + (-t14 * t25 + t7 * t29) * r_i_i_C(2); 0, (t33 * t27 + t34 * t31) * t23, -t43, t47 * t12 + (-t41 * t26 - t30 * t43) * pkin(4) + t36 * (-t41 * t20 - t21 * t43), t45, (-t12 * t25 + t29 * t45) * r_i_i_C(1) + (-t12 * t29 - t25 * t45) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end