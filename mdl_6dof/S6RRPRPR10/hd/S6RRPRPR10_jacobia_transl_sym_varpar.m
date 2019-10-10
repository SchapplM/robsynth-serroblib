% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR10_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR10_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR10_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
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
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (51->25), mult. (121->42), div. (0->0), fcn. (150->8), ass. (0->19)
	t10 = sin(qJ(1));
	t7 = sin(pkin(6));
	t19 = t10 * t7;
	t12 = cos(qJ(1));
	t18 = t12 * t7;
	t17 = r_i_i_C(3) + qJ(3);
	t16 = cos(pkin(6));
	t15 = t10 * t16;
	t14 = t12 * t16;
	t6 = sin(pkin(11));
	t8 = cos(pkin(11));
	t13 = t8 * r_i_i_C(1) - t6 * r_i_i_C(2) + pkin(2);
	t11 = cos(qJ(2));
	t9 = sin(qJ(2));
	t4 = t12 * t11 - t9 * t15;
	t3 = t11 * t15 + t12 * t9;
	t2 = t10 * t11 + t9 * t14;
	t1 = t10 * t9 - t11 * t14;
	t5 = [(t6 * t18 - t2 * t8) * r_i_i_C(1) + (t8 * t18 + t2 * t6) * r_i_i_C(2) - t2 * pkin(2) - t10 * pkin(1) + pkin(8) * t18 - t17 * t1, -t13 * t3 + t17 * t4, t3, 0, 0, 0; (t6 * t19 + t4 * t8) * r_i_i_C(1) + (t8 * t19 - t4 * t6) * r_i_i_C(2) + t4 * pkin(2) + t12 * pkin(1) + pkin(8) * t19 + t17 * t3, -t13 * t1 + t17 * t2, t1, 0, 0, 0; 0, (t13 * t11 + t17 * t9) * t7, -t7 * t11, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (102->33), mult. (168->53), div. (0->0), fcn. (207->10), ass. (0->28)
	t30 = r_i_i_C(3) + pkin(9) + qJ(3);
	t14 = sin(pkin(6));
	t17 = sin(qJ(2));
	t29 = t14 * t17;
	t18 = sin(qJ(1));
	t28 = t14 * t18;
	t20 = cos(qJ(1));
	t27 = t14 * t20;
	t26 = t18 * t17;
	t19 = cos(qJ(2));
	t25 = t18 * t19;
	t24 = t20 * t17;
	t23 = t20 * t19;
	t22 = pkin(3) * sin(pkin(11)) + pkin(8);
	t12 = pkin(11) + qJ(4);
	t10 = sin(t12);
	t11 = cos(t12);
	t9 = cos(pkin(11)) * pkin(3) + pkin(2);
	t21 = t11 * r_i_i_C(1) - t10 * r_i_i_C(2) + t9;
	t15 = cos(pkin(6));
	t7 = t10 * t27;
	t6 = -t15 * t26 + t23;
	t5 = t15 * t25 + t24;
	t4 = t15 * t24 + t25;
	t3 = -t15 * t23 + t26;
	t2 = t10 * t28 + t6 * t11;
	t1 = -t6 * t10 + t11 * t28;
	t8 = [-t18 * pkin(1) + t7 * r_i_i_C(1) - t21 * t4 - t30 * t3 + (t11 * r_i_i_C(2) + t22) * t27, -t21 * t5 + t30 * t6, t5, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t20 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t22 * t28 + t30 * t5 + t6 * t9, -t21 * t3 + t30 * t4, t3, (-t4 * t10 - t11 * t27) * r_i_i_C(1) + (-t4 * t11 + t7) * r_i_i_C(2), 0, 0; 0, (t30 * t17 + t21 * t19) * t14, -t14 * t19, (-t10 * t29 + t15 * t11) * r_i_i_C(1) + (-t15 * t10 - t11 * t29) * r_i_i_C(2), 0, 0;];
	Ja_transl = t8;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (177->35), mult. (273->56), div. (0->0), fcn. (342->10), ass. (0->30)
	t15 = cos(pkin(11)) * pkin(3) + pkin(2);
	t18 = pkin(11) + qJ(4);
	t16 = sin(t18);
	t17 = cos(t18);
	t32 = r_i_i_C(3) + qJ(5);
	t37 = pkin(4) - r_i_i_C(2);
	t38 = t32 * t16 + t37 * t17 + t15;
	t36 = r_i_i_C(1) + pkin(9) + qJ(3);
	t20 = sin(pkin(6));
	t22 = sin(qJ(2));
	t35 = t20 * t22;
	t23 = sin(qJ(1));
	t34 = t20 * t23;
	t25 = cos(qJ(1));
	t33 = t20 * t25;
	t31 = cos(pkin(6));
	t30 = t23 * t31;
	t29 = t25 * t31;
	t28 = t20 * (pkin(3) * sin(pkin(11)) + pkin(8));
	t24 = cos(qJ(2));
	t10 = t22 * t29 + t23 * t24;
	t1 = t10 * t16 + t17 * t33;
	t27 = -t10 * t17 + t16 * t33;
	t12 = -t22 * t30 + t25 * t24;
	t11 = t25 * t22 + t24 * t30;
	t9 = t23 * t22 - t24 * t29;
	t7 = t16 * t35 - t31 * t17;
	t6 = t12 * t17 + t16 * t34;
	t5 = t12 * t16 - t17 * t34;
	t2 = [-t23 * pkin(1) - t32 * t1 - t10 * t15 + t25 * t28 + t37 * t27 - t36 * t9, -t11 * t38 + t36 * t12, t11, t32 * t6 - t37 * t5, t5, 0; t25 * pkin(1) + t36 * t11 + t12 * t15 + t23 * t28 + t32 * t5 + t37 * t6, t36 * t10 - t38 * t9, t9, -t37 * t1 - t32 * t27, t1, 0; 0, (t36 * t22 + t38 * t24) * t20, -t20 * t24, t32 * (t31 * t16 + t17 * t35) - t37 * t7, t7, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (294->49), mult. (477->80), div. (0->0), fcn. (606->12), ass. (0->36)
	t18 = cos(pkin(11)) * pkin(3) + pkin(2);
	t21 = pkin(11) + qJ(4);
	t19 = sin(t21);
	t20 = cos(t21);
	t25 = sin(qJ(6));
	t28 = cos(qJ(6));
	t34 = t25 * r_i_i_C(1) + t28 * r_i_i_C(2) + qJ(5);
	t38 = pkin(4) + pkin(10) + r_i_i_C(3);
	t45 = t34 * t19 + t38 * t20 + t18;
	t44 = pkin(5) + pkin(9) + qJ(3);
	t23 = sin(pkin(6));
	t26 = sin(qJ(2));
	t43 = t23 * t26;
	t27 = sin(qJ(1));
	t42 = t23 * t27;
	t29 = cos(qJ(2));
	t41 = t23 * t29;
	t30 = cos(qJ(1));
	t40 = t23 * t30;
	t39 = cos(pkin(6));
	t37 = t27 * t39;
	t36 = t30 * t39;
	t35 = t23 * (pkin(3) * sin(pkin(11)) + pkin(8));
	t12 = t26 * t36 + t27 * t29;
	t3 = t12 * t19 + t20 * t40;
	t33 = -t12 * t20 + t19 * t40;
	t32 = t28 * r_i_i_C(1) - t25 * r_i_i_C(2) + t44;
	t14 = -t26 * t37 + t30 * t29;
	t13 = t30 * t26 + t29 * t37;
	t11 = t27 * t26 - t29 * t36;
	t9 = t19 * t43 - t39 * t20;
	t8 = t14 * t20 + t19 * t42;
	t7 = t14 * t19 - t20 * t42;
	t2 = t13 * t28 + t7 * t25;
	t1 = -t13 * t25 + t7 * t28;
	t4 = [-t27 * pkin(1) - t32 * t11 - t12 * t18 - t34 * t3 + t30 * t35 + t38 * t33, -t13 * t45 + t32 * t14, t13, t34 * t8 - t38 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t30 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(5) + t44 * t13 + t14 * t18 + t27 * t35 + t38 * t8, -t11 * t45 + t32 * t12, t11, -t38 * t3 - t34 * t33, t3, (-t11 * t25 + t3 * t28) * r_i_i_C(1) + (-t11 * t28 - t3 * t25) * r_i_i_C(2); 0, (t32 * t26 + t45 * t29) * t23, -t41, -t38 * t9 + t34 * (t39 * t19 + t20 * t43), t9, (t25 * t41 + t9 * t28) * r_i_i_C(1) + (-t9 * t25 + t28 * t41) * r_i_i_C(2);];
	Ja_transl = t4;
else
	Ja_transl=NaN(3,6);
end