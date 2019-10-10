% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPPR1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
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
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t2 * t5 + t4 * t8, t6 * t4, 0, 0, 0, 0; t2 * t8 + t4 * t5, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (29->10), mult. (31->12), div. (0->0), fcn. (33->6), ass. (0->10)
	t5 = qJ(2) + pkin(9);
	t2 = sin(t5);
	t3 = cos(t5);
	t14 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2) + cos(qJ(2)) * pkin(2);
	t13 = r_i_i_C(3) + qJ(3) + pkin(7);
	t11 = pkin(1) + t14;
	t10 = -sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t11 * t8 + t13 * t9, t10 * t9, t8, 0, 0, 0; t11 * t9 + t13 * t8, t10 * t8, -t9, 0, 0, 0; 0, t14, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (65->16), mult. (72->18), div. (0->0), fcn. (81->8), ass. (0->14)
	t6 = sin(pkin(10));
	t7 = cos(pkin(10));
	t15 = r_i_i_C(1) * t7 - r_i_i_C(2) * t6 + pkin(3);
	t16 = r_i_i_C(3) + qJ(4);
	t5 = qJ(2) + pkin(9);
	t2 = sin(t5);
	t3 = cos(t5);
	t18 = cos(qJ(2)) * pkin(2) + t15 * t3 + t16 * t2;
	t17 = pkin(1) + t18;
	t14 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7 + pkin(7) + qJ(3);
	t12 = -sin(qJ(2)) * pkin(2) - t15 * t2 + t16 * t3;
	t11 = cos(qJ(1));
	t10 = sin(qJ(1));
	t1 = [-t17 * t10 + t14 * t11, t12 * t11, t10, t11 * t2, 0, 0; t14 * t10 + t17 * t11, t12 * t10, -t11, t10 * t2, 0, 0; 0, t18, 0, -t3, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (94->25), mult. (111->30), div. (0->0), fcn. (129->8), ass. (0->22)
	t20 = r_i_i_C(2) + qJ(4);
	t9 = qJ(2) + pkin(9);
	t6 = sin(t9);
	t28 = t20 * t6 + cos(qJ(2)) * pkin(2);
	t7 = cos(t9);
	t27 = t7 * pkin(3) + pkin(1) + t28;
	t10 = sin(pkin(10));
	t11 = cos(pkin(10));
	t19 = r_i_i_C(3) + qJ(5);
	t25 = pkin(4) + r_i_i_C(1);
	t26 = t19 * t10 + t25 * t11 + pkin(3);
	t15 = cos(qJ(1));
	t24 = t10 * t15;
	t23 = t11 * t15;
	t14 = sin(qJ(1));
	t22 = t14 * t10;
	t21 = t14 * t11;
	t16 = -sin(qJ(2)) * pkin(2) + t20 * t7 - t26 * t6;
	t12 = -qJ(3) - pkin(7);
	t3 = t7 * t24 - t21;
	t1 = t7 * t22 + t23;
	t2 = [-t12 * t15 + t25 * (-t7 * t21 + t24) - t19 * t1 - t27 * t14, t16 * t15, t14, t15 * t6, t3, 0; -t14 * t12 + t25 * (t7 * t23 + t22) + t19 * t3 + t27 * t15, t16 * t14, -t15, t14 * t6, t1, 0; 0, t26 * t7 + t28, 0, -t7, t6 * t10, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:18:08
	% EndTime: 2019-10-10 09:18:08
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (169->39), mult. (224->55), div. (0->0), fcn. (272->10), ass. (0->29)
	t26 = -r_i_i_C(3) - pkin(8) + qJ(4);
	t12 = qJ(2) + pkin(9);
	t9 = sin(t12);
	t34 = cos(qJ(2)) * pkin(2) + t26 * t9;
	t10 = cos(t12);
	t33 = t10 * pkin(3) + pkin(1) + t34;
	t13 = sin(pkin(10));
	t14 = cos(pkin(10));
	t16 = sin(qJ(6));
	t19 = cos(qJ(6));
	t31 = pkin(4) + pkin(5);
	t22 = t19 * r_i_i_C(1) - t16 * r_i_i_C(2) + t31;
	t23 = t16 * r_i_i_C(1) + t19 * r_i_i_C(2) + qJ(5);
	t32 = t23 * t13 + t22 * t14 + pkin(3);
	t18 = sin(qJ(1));
	t30 = t18 * t13;
	t29 = t18 * t14;
	t20 = cos(qJ(1));
	t28 = t20 * t13;
	t27 = t20 * t14;
	t21 = -sin(qJ(2)) * pkin(2) + t26 * t10 - t32 * t9;
	t15 = -qJ(3) - pkin(7);
	t6 = t10 * t27 + t30;
	t5 = t10 * t28 - t29;
	t4 = t10 * t29 - t28;
	t3 = t10 * t30 + t27;
	t2 = t5 * t16 + t6 * t19;
	t1 = -t6 * t16 + t5 * t19;
	t7 = [-t20 * t15 - t33 * t18 - t22 * t4 - t23 * t3, t21 * t20, t18, t20 * t9, t5, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * qJ(5) - t18 * t15 + t33 * t20 + t31 * t6, t21 * t18, -t20, t18 * t9, t3, (-t4 * t16 + t3 * t19) * r_i_i_C(1) + (-t3 * t16 - t4 * t19) * r_i_i_C(2); 0, t32 * t10 + t34, 0, -t10, t9 * t13, ((t13 * t19 - t14 * t16) * r_i_i_C(1) + (-t13 * t16 - t14 * t19) * r_i_i_C(2)) * t9;];
	Ja_transl = t7;
else
	Ja_transl=NaN(3,6);
end