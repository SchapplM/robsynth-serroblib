% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:15
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:23
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:23
	% EndTime: 2019-10-10 10:15:24
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
	% StartTime: 2019-10-10 10:15:23
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-10-10 10:15:23
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (22->9), mult. (44->12), div. (0->0), fcn. (47->4), ass. (0->11)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(3) + qJ(3);
	t9 = pkin(2) + r_i_i_C(1);
	t6 = t7 * t1 + t3 * t9;
	t10 = pkin(1) + t6;
	t8 = pkin(7) + r_i_i_C(2);
	t5 = -t1 * t9 + t7 * t3;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t11 = [-t10 * t2 + t8 * t4, t5 * t4, t4 * t1, 0, 0, 0; t10 * t4 + t8 * t2, t5 * t2, t2 * t1, 0, 0, 0; 0, t6, -t3, 0, 0, 0;];
	Ja_transl = t11;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:23
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (48->21), mult. (109->31), div. (0->0), fcn. (128->6), ass. (0->23)
	t12 = cos(qJ(2));
	t24 = pkin(2) + pkin(3);
	t9 = sin(qJ(2));
	t14 = t9 * qJ(3) + t24 * t12;
	t25 = pkin(1) + t14;
	t8 = sin(qJ(4));
	t23 = t12 * t8;
	t13 = cos(qJ(1));
	t22 = t13 * t9;
	t20 = pkin(7) - pkin(8) - r_i_i_C(3);
	t10 = sin(qJ(1));
	t11 = cos(qJ(4));
	t6 = t11 * t9 - t23;
	t1 = t6 * t10;
	t16 = t11 * t12 + t8 * t9;
	t2 = t16 * t10;
	t19 = t1 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t3 = -t11 * t22 + t13 * t23;
	t4 = t16 * t13;
	t18 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t17 = -r_i_i_C(1) * t16 - t6 * r_i_i_C(2);
	t15 = qJ(3) * t12 - t24 * t9;
	t5 = [-t2 * r_i_i_C(1) - t1 * r_i_i_C(2) - t25 * t10 + t20 * t13, t15 * t13 - t18, t22, t18, 0, 0; t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t20 * t10 + t25 * t13, t15 * t10 - t19, t10 * t9, t19, 0, 0; 0, t14 - t17, -t12, t17, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:23
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (99->30), mult. (144->41), div. (0->0), fcn. (165->8), ass. (0->28)
	t14 = sin(qJ(2));
	t17 = cos(qJ(2));
	t13 = sin(qJ(4));
	t26 = pkin(4) * t13 + qJ(3);
	t16 = cos(qJ(4));
	t30 = pkin(4) * t16 + pkin(2) + pkin(3);
	t20 = t26 * t14 + t30 * t17;
	t31 = pkin(1) + t20;
	t11 = qJ(4) + pkin(10);
	t9 = sin(t11);
	t29 = t17 * t9;
	t18 = cos(qJ(1));
	t28 = t18 * t14;
	t27 = pkin(7) - r_i_i_C(3) - qJ(5) - pkin(8);
	t15 = sin(qJ(1));
	t10 = cos(t11);
	t6 = t10 * t14 - t29;
	t1 = t6 * t15;
	t22 = t10 * t17 + t14 * t9;
	t2 = t22 * t15;
	t25 = r_i_i_C(1) * t1 - r_i_i_C(2) * t2;
	t3 = -t10 * t28 + t18 * t29;
	t4 = t22 * t18;
	t24 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t23 = -r_i_i_C(1) * t22 - t6 * r_i_i_C(2);
	t21 = pkin(4) * (-t13 * t17 + t14 * t16);
	t19 = -t30 * t14 + t26 * t17;
	t5 = [-t2 * r_i_i_C(1) - t1 * r_i_i_C(2) - t31 * t15 + t27 * t18, t19 * t18 - t24, t28, t18 * t21 + t24, -t15, 0; t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t27 * t15 + t31 * t18, t19 * t15 - t25, t15 * t14, t15 * t21 + t25, t18, 0; 0, t20 - t23, -t17, (-t13 * t14 - t16 * t17) * pkin(4) + t23, 0, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:15:23
	% EndTime: 2019-10-10 10:15:24
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (222->44), mult. (306->63), div. (0->0), fcn. (365->10), ass. (0->35)
	t18 = sin(qJ(2));
	t22 = cos(qJ(2));
	t17 = sin(qJ(4));
	t35 = pkin(4) * t17 + qJ(3);
	t21 = cos(qJ(4));
	t38 = t21 * pkin(4) + pkin(2) + pkin(3);
	t25 = t35 * t18 + t38 * t22;
	t43 = pkin(1) + t25;
	t36 = qJ(4) + pkin(10);
	t33 = sin(t36);
	t34 = cos(t36);
	t8 = t18 * t34 - t22 * t33;
	t16 = sin(qJ(6));
	t20 = cos(qJ(6));
	t27 = t20 * r_i_i_C(1) - t16 * r_i_i_C(2) + pkin(5);
	t19 = sin(qJ(1));
	t3 = t8 * t19;
	t39 = -r_i_i_C(3) - pkin(9);
	t7 = t18 * t33 + t22 * t34;
	t4 = t7 * t19;
	t42 = -t27 * t3 + t39 * t4;
	t23 = cos(qJ(1));
	t28 = t23 * t33;
	t29 = t23 * t34;
	t5 = -t18 * t28 - t22 * t29;
	t6 = -t18 * t29 + t22 * t28;
	t41 = t27 * t6 - t39 * t5;
	t40 = -t27 * t7 - t39 * t8;
	t37 = pkin(7) - qJ(5) - pkin(8);
	t32 = -t16 * r_i_i_C(1) - t20 * r_i_i_C(2);
	t26 = pkin(4) * (-t17 * t22 + t18 * t21);
	t24 = -t38 * t18 + t35 * t22;
	t2 = -t19 * t16 - t5 * t20;
	t1 = t5 * t16 - t19 * t20;
	t9 = [-t39 * t3 - t27 * t4 + (t32 + t37) * t23 - t43 * t19, t24 * t23 + t41, t23 * t18, t23 * t26 - t41, -t19, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -t5 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t37 * t19 + t43 * t23 - t39 * t6, t24 * t19 + t42, t19 * t18, t19 * t26 - t42, t23, (-t4 * t16 + t23 * t20) * r_i_i_C(1) + (-t23 * t16 - t4 * t20) * r_i_i_C(2); 0, t25 - t40, -t22, (-t17 * t18 - t21 * t22) * pkin(4) + t40, 0, t32 * t8;];
	Ja_transl = t9;
else
	Ja_transl=NaN(3,6);
end