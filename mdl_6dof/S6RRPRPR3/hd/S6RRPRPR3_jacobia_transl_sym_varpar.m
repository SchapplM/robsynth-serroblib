% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:52
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:52
	% EndTime: 2019-10-10 10:07:53
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
	% StartTime: 2019-10-10 10:07:52
	% EndTime: 2019-10-10 10:07:53
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
	% StartTime: 2019-10-10 10:07:52
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (29->10), mult. (31->12), div. (0->0), fcn. (33->6), ass. (0->10)
	t5 = qJ(2) + pkin(10);
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
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (74->25), mult. (90->34), div. (0->0), fcn. (100->8), ass. (0->22)
	t24 = pkin(8) + r_i_i_C(3);
	t9 = qJ(2) + pkin(10);
	t6 = sin(t9);
	t26 = t24 * t6 + cos(qJ(2)) * pkin(2);
	t7 = cos(t9);
	t25 = t7 * pkin(3) + pkin(1) + t26;
	t11 = sin(qJ(4));
	t15 = cos(qJ(1));
	t23 = t11 * t15;
	t13 = sin(qJ(1));
	t22 = t13 * t11;
	t14 = cos(qJ(4));
	t21 = t13 * t14;
	t20 = t14 * t15;
	t17 = r_i_i_C(1) * t14 - r_i_i_C(2) * t11 + pkin(3);
	t16 = -sin(qJ(2)) * pkin(2) - t17 * t6 + t24 * t7;
	t10 = -qJ(3) - pkin(7);
	t4 = t7 * t20 + t22;
	t3 = -t7 * t23 + t21;
	t2 = -t7 * t21 + t23;
	t1 = t7 * t22 + t20;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t15 * t10 - t25 * t13, t16 * t15, t13, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t13 * t10 + t25 * t15, t16 * t13, -t15, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t17 * t7 + t26, 0, (-r_i_i_C(1) * t11 - r_i_i_C(2) * t14) * t6, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:52
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (122->33), mult. (113->46), div. (0->0), fcn. (126->10), ass. (0->25)
	t28 = r_i_i_C(3) + qJ(5) + pkin(8);
	t13 = qJ(2) + pkin(10);
	t8 = sin(t13);
	t31 = cos(qJ(2)) * pkin(2) + t28 * t8;
	t10 = cos(t13);
	t19 = cos(qJ(4));
	t5 = pkin(4) * t19 + pkin(3);
	t30 = t10 * t5 + pkin(1) + t31;
	t16 = sin(qJ(4));
	t29 = pkin(4) * t16;
	t18 = sin(qJ(1));
	t27 = t10 * t18;
	t20 = cos(qJ(1));
	t26 = t10 * t20;
	t23 = qJ(3) + pkin(7) + t29;
	t12 = qJ(4) + pkin(11);
	t7 = sin(t12);
	t9 = cos(t12);
	t22 = r_i_i_C(1) * t9 - r_i_i_C(2) * t7 + t5;
	t21 = -sin(qJ(2)) * pkin(2) + t28 * t10 - t22 * t8;
	t4 = t18 * t7 + t9 * t26;
	t3 = t18 * t9 - t7 * t26;
	t2 = t20 * t7 - t9 * t27;
	t1 = t20 * t9 + t7 * t27;
	t6 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t30 * t18 + t23 * t20, t21 * t20, t18, t3 * r_i_i_C(1) - t4 * r_i_i_C(2) + (-t16 * t26 + t18 * t19) * pkin(4), t20 * t8, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t23 * t18 + t30 * t20, t21 * t18, -t20, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + (-t16 * t27 - t19 * t20) * pkin(4), t18 * t8, 0; 0, t22 * t10 + t31, 0, (-r_i_i_C(1) * t7 - r_i_i_C(2) * t9 - t29) * t8, -t10, 0;];
	Ja_transl = t6;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:07:53
	% EndTime: 2019-10-10 10:07:53
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (208->36), mult. (145->47), div. (0->0), fcn. (162->12), ass. (0->29)
	t23 = qJ(2) + pkin(10);
	t16 = sin(t23);
	t37 = r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
	t41 = cos(qJ(2)) * pkin(2) + t37 * t16;
	t17 = cos(t23);
	t22 = qJ(4) + pkin(11);
	t11 = pkin(5) * cos(t22) + cos(qJ(4)) * pkin(4);
	t9 = pkin(3) + t11;
	t40 = t17 * t9 + pkin(1) + t41;
	t18 = qJ(6) + t22;
	t14 = cos(t18);
	t27 = cos(qJ(1));
	t13 = sin(t18);
	t26 = sin(qJ(1));
	t35 = t26 * t13;
	t5 = t14 * t27 + t17 * t35;
	t34 = t26 * t14;
	t6 = t13 * t27 - t17 * t34;
	t39 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t36 = t17 * t27;
	t7 = -t13 * t36 + t34;
	t8 = t14 * t36 + t35;
	t38 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t10 = pkin(5) * sin(t22) + sin(qJ(4)) * pkin(4);
	t33 = t10 + qJ(3) + pkin(7);
	t30 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
	t29 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
	t28 = -sin(qJ(2)) * pkin(2) - t16 * t29 + t17 * t37;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t40 * t26 + t33 * t27, t28 * t27, t26, -t10 * t36 + t26 * t11 + t38, t27 * t16, t38; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t33 * t26 + t40 * t27, t28 * t26, -t27, -t26 * t17 * t10 - t11 * t27 + t39, t26 * t16, t39; 0, t17 * t29 + t41, 0, (-t10 + t30) * t16, -t17, t30 * t16;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end