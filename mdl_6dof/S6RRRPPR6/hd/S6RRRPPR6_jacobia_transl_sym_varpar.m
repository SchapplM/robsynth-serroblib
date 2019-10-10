% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:26
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:55
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:55
	% EndTime: 2019-10-10 11:25:56
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
	% StartTime: 2019-10-10 11:25:55
	% EndTime: 2019-10-10 11:25:56
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
	% StartTime: 2019-10-10 11:25:55
	% EndTime: 2019-10-10 11:25:56
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
	% StartTime: 2019-10-10 11:25:55
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (110->39), mult. (189->62), div. (0->0), fcn. (230->10), ass. (0->30)
	t31 = r_i_i_C(3) + qJ(4) + pkin(9);
	t13 = sin(pkin(6));
	t17 = sin(qJ(2));
	t30 = t13 * t17;
	t18 = sin(qJ(1));
	t29 = t13 * t18;
	t21 = cos(qJ(1));
	t28 = t13 * t21;
	t27 = t17 * t18;
	t26 = t17 * t21;
	t20 = cos(qJ(2));
	t25 = t18 * t20;
	t24 = t20 * t21;
	t16 = sin(qJ(3));
	t23 = t16 * pkin(3) + pkin(8);
	t12 = qJ(3) + pkin(11);
	t10 = sin(t12);
	t11 = cos(t12);
	t19 = cos(qJ(3));
	t9 = t19 * pkin(3) + pkin(2);
	t22 = t11 * r_i_i_C(1) - t10 * r_i_i_C(2) + t9;
	t14 = cos(pkin(6));
	t7 = t10 * t28;
	t6 = -t14 * t27 + t24;
	t5 = t14 * t25 + t26;
	t4 = t14 * t26 + t25;
	t3 = -t14 * t24 + t27;
	t2 = t10 * t29 + t11 * t6;
	t1 = -t10 * t6 + t11 * t29;
	t8 = [-t18 * pkin(1) + t7 * r_i_i_C(1) - t22 * t4 - t31 * t3 + (t11 * r_i_i_C(2) + t23) * t28, -t22 * t5 + t31 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2 + (-t16 * t6 + t19 * t29) * pkin(3), t5, 0, 0; pkin(1) * t21 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t23 * t29 + t31 * t5 + t6 * t9, -t22 * t3 + t31 * t4, (-t4 * t10 - t11 * t28) * r_i_i_C(1) + (-t11 * t4 + t7) * r_i_i_C(2) + (-t4 * t16 - t19 * t28) * pkin(3), t3, 0, 0; 0, (t31 * t17 + t22 * t20) * t13, (-t10 * t30 + t11 * t14) * r_i_i_C(1) + (-t10 * t14 - t11 * t30) * r_i_i_C(2) + (t14 * t19 - t16 * t30) * pkin(3), -t13 * t20, 0, 0;];
	Ja_transl = t8;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:55
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (185->41), mult. (294->63), div. (0->0), fcn. (365->10), ass. (0->34)
	t25 = cos(qJ(3));
	t15 = t25 * pkin(3) + pkin(2);
	t18 = qJ(3) + pkin(11);
	t16 = sin(t18);
	t17 = cos(t18);
	t31 = r_i_i_C(3) + qJ(5);
	t40 = pkin(4) - r_i_i_C(2);
	t41 = t31 * t16 + t40 * t17 + t15;
	t39 = r_i_i_C(1) + qJ(4) + pkin(9);
	t19 = sin(pkin(6));
	t23 = sin(qJ(2));
	t38 = t19 * t23;
	t24 = sin(qJ(1));
	t37 = t19 * t24;
	t27 = cos(qJ(1));
	t36 = t19 * t27;
	t35 = t24 * t23;
	t26 = cos(qJ(2));
	t34 = t24 * t26;
	t33 = t27 * t23;
	t32 = t27 * t26;
	t22 = sin(qJ(3));
	t30 = t19 * (pkin(3) * t22 + pkin(8));
	t20 = cos(pkin(6));
	t10 = t20 * t33 + t34;
	t1 = t10 * t16 + t17 * t36;
	t29 = -t10 * t17 + t16 * t36;
	t12 = -t20 * t35 + t32;
	t11 = t20 * t34 + t33;
	t9 = -t20 * t32 + t35;
	t7 = t16 * t38 - t20 * t17;
	t6 = t12 * t17 + t16 * t37;
	t5 = t12 * t16 - t17 * t37;
	t2 = [-t24 * pkin(1) - t31 * t1 - t10 * t15 + t27 * t30 + t40 * t29 - t39 * t9, -t11 * t41 + t39 * t12, t31 * t6 - t40 * t5 + (-t12 * t22 + t25 * t37) * pkin(3), t11, t5, 0; t27 * pkin(1) + t39 * t11 + t12 * t15 + t24 * t30 + t31 * t5 + t40 * t6, t39 * t10 - t41 * t9, -t31 * t29 - t40 * t1 + (-t10 * t22 - t25 * t36) * pkin(3), t9, t1, 0; 0, (t39 * t23 + t41 * t26) * t19, t31 * (t20 * t16 + t17 * t38) - t40 * t7 + (t20 * t25 - t22 * t38) * pkin(3), -t19 * t26, t7, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:55
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (302->55), mult. (498->89), div. (0->0), fcn. (629->12), ass. (0->38)
	t29 = cos(qJ(3));
	t18 = t29 * pkin(3) + pkin(2);
	t21 = qJ(3) + pkin(11);
	t19 = sin(t21);
	t20 = cos(t21);
	t24 = sin(qJ(6));
	t28 = cos(qJ(6));
	t35 = t24 * r_i_i_C(1) + t28 * r_i_i_C(2) + qJ(5);
	t39 = pkin(4) + pkin(10) + r_i_i_C(3);
	t46 = t35 * t19 + t39 * t20 + t18;
	t45 = pkin(5) + qJ(4) + pkin(9);
	t22 = sin(pkin(6));
	t26 = sin(qJ(2));
	t44 = t22 * t26;
	t27 = sin(qJ(1));
	t43 = t22 * t27;
	t30 = cos(qJ(2));
	t42 = t22 * t30;
	t31 = cos(qJ(1));
	t41 = t22 * t31;
	t40 = cos(pkin(6));
	t38 = t27 * t40;
	t37 = t31 * t40;
	t25 = sin(qJ(3));
	t36 = t22 * (pkin(3) * t25 + pkin(8));
	t12 = t26 * t37 + t27 * t30;
	t3 = t12 * t19 + t20 * t41;
	t34 = -t12 * t20 + t19 * t41;
	t33 = t28 * r_i_i_C(1) - t24 * r_i_i_C(2) + t45;
	t14 = -t26 * t38 + t31 * t30;
	t13 = t31 * t26 + t30 * t38;
	t11 = t27 * t26 - t30 * t37;
	t9 = t19 * t44 - t40 * t20;
	t8 = t14 * t20 + t19 * t43;
	t7 = t14 * t19 - t20 * t43;
	t2 = t13 * t28 + t7 * t24;
	t1 = -t13 * t24 + t7 * t28;
	t4 = [-t27 * pkin(1) - t33 * t11 - t12 * t18 - t35 * t3 + t31 * t36 + t39 * t34, -t13 * t46 + t33 * t14, (-t14 * t25 + t29 * t43) * pkin(3) + t35 * t8 - t39 * t7, t13, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t31 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(5) + t45 * t13 + t14 * t18 + t27 * t36 + t39 * t8, -t11 * t46 + t33 * t12, (-t12 * t25 - t29 * t41) * pkin(3) - t35 * t34 - t39 * t3, t11, t3, (-t11 * t24 + t3 * t28) * r_i_i_C(1) + (-t11 * t28 - t3 * t24) * r_i_i_C(2); 0, (t33 * t26 + t46 * t30) * t22, (-t25 * t44 + t40 * t29) * pkin(3) - t39 * t9 + t35 * (t40 * t19 + t20 * t44), -t42, t9, (t24 * t42 + t9 * t28) * r_i_i_C(1) + (-t9 * t24 + t28 * t42) * r_i_i_C(2);];
	Ja_transl = t4;
else
	Ja_transl=NaN(3,6);
end