% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR14
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPRR14_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR14_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR14_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:24:51
	% EndTime: 2019-12-29 19:24:51
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:24:58
	% EndTime: 2019-12-29 19:24:58
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:24:58
	% EndTime: 2019-12-29 19:24:58
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (19->14), mult. (46->24), div. (0->0), fcn. (54->6), ass. (0->16)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t15 = t8 * t7;
	t9 = cos(qJ(2));
	t14 = t8 * t9;
	t10 = cos(qJ(1));
	t13 = t10 * t7;
	t12 = t10 * t9;
	t5 = sin(pkin(5));
	t11 = (pkin(7) + r_i_i_C(3)) * t5;
	t6 = cos(pkin(5));
	t4 = -t15 * t6 + t12;
	t3 = -t14 * t6 - t13;
	t2 = -t13 * t6 - t14;
	t1 = -t12 * t6 + t15;
	t16 = [-t8 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t10 * t11, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0; t10 * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t11 * t8, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0; 0, (r_i_i_C(1) * t9 - r_i_i_C(2) * t7) * t5, 0, 0, 0;];
	Ja_transl = t16;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:24:53
	% EndTime: 2019-12-29 19:24:54
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (51->25), mult. (121->42), div. (0->0), fcn. (150->8), ass. (0->19)
	t10 = sin(qJ(1));
	t7 = sin(pkin(5));
	t19 = t10 * t7;
	t12 = cos(qJ(1));
	t18 = t12 * t7;
	t17 = r_i_i_C(3) + qJ(3);
	t16 = cos(pkin(5));
	t15 = t10 * t16;
	t14 = t12 * t16;
	t6 = sin(pkin(10));
	t8 = cos(pkin(10));
	t13 = t8 * r_i_i_C(1) - t6 * r_i_i_C(2) + pkin(2);
	t11 = cos(qJ(2));
	t9 = sin(qJ(2));
	t4 = t12 * t11 - t9 * t15;
	t3 = t11 * t15 + t12 * t9;
	t2 = t10 * t11 + t9 * t14;
	t1 = t10 * t9 - t11 * t14;
	t5 = [(t6 * t18 - t2 * t8) * r_i_i_C(1) + (t8 * t18 + t2 * t6) * r_i_i_C(2) - t2 * pkin(2) - t10 * pkin(1) + pkin(7) * t18 - t17 * t1, -t13 * t3 + t17 * t4, t3, 0, 0; (t6 * t19 + t4 * t8) * r_i_i_C(1) + (t8 * t19 - t4 * t6) * r_i_i_C(2) + t4 * pkin(2) + t12 * pkin(1) + pkin(7) * t19 + t17 * t3, -t13 * t1 + t17 * t2, t1, 0, 0; 0, (t13 * t11 + t17 * t9) * t7, -t7 * t11, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:24:51
	% EndTime: 2019-12-29 19:24:51
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (102->33), mult. (168->53), div. (0->0), fcn. (207->10), ass. (0->28)
	t30 = r_i_i_C(3) + pkin(8) + qJ(3);
	t14 = sin(pkin(5));
	t17 = sin(qJ(2));
	t29 = t14 * t17;
	t18 = sin(qJ(1));
	t28 = t14 * t18;
	t20 = cos(qJ(1));
	t27 = t14 * t20;
	t26 = t17 * t18;
	t25 = t17 * t20;
	t19 = cos(qJ(2));
	t24 = t18 * t19;
	t23 = t19 * t20;
	t22 = pkin(3) * sin(pkin(10)) + pkin(7);
	t12 = pkin(10) + qJ(4);
	t10 = sin(t12);
	t11 = cos(t12);
	t9 = cos(pkin(10)) * pkin(3) + pkin(2);
	t21 = t11 * r_i_i_C(1) - t10 * r_i_i_C(2) + t9;
	t15 = cos(pkin(5));
	t7 = t10 * t27;
	t6 = -t15 * t26 + t23;
	t5 = t15 * t24 + t25;
	t4 = t15 * t25 + t24;
	t3 = -t15 * t23 + t26;
	t2 = t10 * t28 + t11 * t6;
	t1 = -t10 * t6 + t11 * t28;
	t8 = [-t18 * pkin(1) + t7 * r_i_i_C(1) - t21 * t4 - t30 * t3 + (t11 * r_i_i_C(2) + t22) * t27, -t21 * t5 + t30 * t6, t5, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; pkin(1) * t20 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t22 * t28 + t30 * t5 + t6 * t9, -t21 * t3 + t30 * t4, t3, (-t4 * t10 - t11 * t27) * r_i_i_C(1) + (-t11 * t4 + t7) * r_i_i_C(2), 0; 0, (t30 * t17 + t21 * t19) * t14, -t14 * t19, (-t10 * t29 + t11 * t15) * r_i_i_C(1) + (-t15 * t10 - t11 * t29) * r_i_i_C(2), 0;];
	Ja_transl = t8;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:24:58
	% EndTime: 2019-12-29 19:24:59
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (242->51), mult. (399->82), div. (0->0), fcn. (505->12), ass. (0->38)
	t26 = sin(qJ(2));
	t27 = sin(qJ(1));
	t29 = cos(qJ(2));
	t30 = cos(qJ(1));
	t38 = cos(pkin(5));
	t36 = t30 * t38;
	t11 = t27 * t26 - t29 * t36;
	t25 = sin(qJ(5));
	t28 = cos(qJ(5));
	t12 = t26 * t36 + t27 * t29;
	t21 = pkin(10) + qJ(4);
	t19 = sin(t21);
	t20 = cos(t21);
	t23 = sin(pkin(5));
	t39 = t23 * t30;
	t4 = t12 * t20 - t19 * t39;
	t48 = -t11 * t28 + t4 * t25;
	t47 = -t11 * t25 - t4 * t28;
	t18 = cos(pkin(10)) * pkin(3) + pkin(2);
	t34 = t28 * r_i_i_C(1) - t25 * r_i_i_C(2) + pkin(4);
	t45 = pkin(9) + r_i_i_C(3);
	t46 = t45 * t19 + t34 * t20 + t18;
	t42 = t23 * t26;
	t41 = t23 * t27;
	t40 = t23 * t29;
	t37 = t27 * t38;
	t35 = t23 * (pkin(3) * sin(pkin(10)) + pkin(7));
	t24 = -pkin(8) - qJ(3);
	t33 = t25 * r_i_i_C(1) + t28 * r_i_i_C(2) - t24;
	t32 = -t12 * t19 - t20 * t39;
	t14 = -t26 * t37 + t30 * t29;
	t13 = t30 * t26 + t29 * t37;
	t10 = t38 * t19 + t20 * t42;
	t8 = t14 * t20 + t19 * t41;
	t7 = t14 * t19 - t20 * t41;
	t2 = t13 * t25 + t8 * t28;
	t1 = t13 * t28 - t8 * t25;
	t3 = [-t27 * pkin(1) - t4 * pkin(4) + t47 * r_i_i_C(1) + t48 * r_i_i_C(2) + t11 * t24 - t12 * t18 + t30 * t35 + t45 * t32, -t13 * t46 + t33 * t14, t13, -t34 * t7 + t45 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t30 * pkin(1) + t8 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t13 * t24 + t14 * t18 + t27 * t35 + t45 * t7, -t11 * t46 + t33 * t12, t11, t34 * t32 + t45 * t4, -t48 * r_i_i_C(1) + t47 * r_i_i_C(2); 0, (t33 * t26 + t29 * t46) * t23, -t40, t45 * t10 + t34 * (-t19 * t42 + t38 * t20), (-t10 * t25 - t28 * t40) * r_i_i_C(1) + (-t10 * t28 + t25 * t40) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,5);
end