% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
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
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
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
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (40->17), mult. (89->25), div. (0->0), fcn. (110->6), ass. (0->18)
	t8 = sin(qJ(2));
	t9 = sin(qJ(1));
	t18 = t8 * t9;
	t17 = pkin(2) + r_i_i_C(1);
	t11 = cos(qJ(1));
	t16 = t11 * t8;
	t10 = cos(qJ(2));
	t15 = t9 * t10;
	t14 = t10 * t11;
	t13 = r_i_i_C(3) + qJ(3);
	t6 = sin(pkin(6));
	t12 = (pkin(8) + r_i_i_C(2)) * t6;
	t7 = cos(pkin(6));
	t4 = -t7 * t18 + t14;
	t3 = t7 * t15 + t16;
	t2 = t7 * t16 + t15;
	t1 = -t7 * t14 + t18;
	t5 = [-t9 * pkin(1) - t13 * t1 + t11 * t12 - t17 * t2, t13 * t4 - t17 * t3, t3, 0, 0, 0; t11 * pkin(1) + t9 * t12 + t13 * t3 + t17 * t4, -t17 * t1 + t13 * t2, t1, 0, 0, 0; 0, (t17 * t10 + t13 * t8) * t6, -t6 * t10, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (53->20), mult. (113->27), div. (0->0), fcn. (141->6), ass. (0->18)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t19 = t10 * t9;
	t12 = cos(qJ(1));
	t18 = t12 * t9;
	t11 = cos(qJ(2));
	t17 = t10 * t11;
	t16 = t11 * t12;
	t15 = r_i_i_C(1) + qJ(3);
	t14 = pkin(2) + pkin(3) - r_i_i_C(2);
	t7 = sin(pkin(6));
	t13 = (pkin(8) - r_i_i_C(3) - qJ(4)) * t7;
	t8 = cos(pkin(6));
	t4 = -t8 * t19 + t16;
	t3 = t8 * t17 + t18;
	t2 = t8 * t18 + t17;
	t1 = -t8 * t16 + t19;
	t5 = [-t10 * pkin(1) - t15 * t1 + t12 * t13 - t14 * t2, -t14 * t3 + t15 * t4, t3, -t10 * t7, 0, 0; pkin(1) * t12 + t10 * t13 + t14 * t4 + t15 * t3, -t14 * t1 + t15 * t2, t1, t12 * t7, 0, 0; 0, (t14 * t11 + t15 * t9) * t7, -t7 * t11, -t8, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (95->34), mult. (217->52), div. (0->0), fcn. (271->8), ass. (0->24)
	t12 = sin(qJ(1));
	t9 = sin(pkin(6));
	t25 = t12 * t9;
	t15 = cos(qJ(1));
	t24 = t15 * t9;
	t14 = cos(qJ(2));
	t23 = t9 * t14;
	t22 = pkin(4) + qJ(3);
	t21 = pkin(8) - qJ(4);
	t20 = cos(pkin(6));
	t19 = -r_i_i_C(3) - pkin(9) - pkin(3) - pkin(2);
	t18 = t12 * t20;
	t17 = t15 * t20;
	t10 = sin(qJ(5));
	t13 = cos(qJ(5));
	t16 = t13 * r_i_i_C(1) - t10 * r_i_i_C(2) + t22;
	t11 = sin(qJ(2));
	t6 = -t11 * t18 + t15 * t14;
	t5 = t15 * t11 + t14 * t18;
	t4 = t11 * t17 + t12 * t14;
	t3 = t12 * t11 - t14 * t17;
	t2 = -t10 * t25 + t5 * t13;
	t1 = -t5 * t10 - t13 * t25;
	t7 = [-t12 * pkin(1) - t16 * t3 + (-t10 * r_i_i_C(1) - t13 * r_i_i_C(2) + t21) * t24 + t19 * t4, t16 * t6 + t19 * t5, t5, -t25, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t15 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t19 * t6 + t21 * t25 + t22 * t5, t16 * t4 + t19 * t3, t3, t24, (-t3 * t10 + t13 * t24) * r_i_i_C(1) + (-t10 * t24 - t3 * t13) * r_i_i_C(2), 0; 0, (t16 * t11 - t19 * t14) * t9, -t23, -t20, (t10 * t23 - t20 * t13) * r_i_i_C(1) + (t20 * t10 + t13 * t23) * r_i_i_C(2), 0;];
	Ja_transl = t7;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (184->50), mult. (448->78), div. (0->0), fcn. (569->10), ass. (0->35)
	t40 = pkin(10) + r_i_i_C(3);
	t17 = sin(pkin(6));
	t20 = sin(qJ(2));
	t39 = t17 * t20;
	t24 = cos(qJ(2));
	t38 = t17 * t24;
	t21 = sin(qJ(1));
	t37 = t21 * t17;
	t25 = cos(qJ(1));
	t36 = t25 * t17;
	t35 = pkin(4) + qJ(3);
	t34 = cos(pkin(6));
	t33 = pkin(2) + pkin(3) + pkin(9);
	t32 = t17 * (pkin(8) - qJ(4));
	t31 = t21 * t34;
	t30 = t25 * t34;
	t18 = sin(qJ(6));
	t22 = cos(qJ(6));
	t29 = t22 * r_i_i_C(1) - t18 * r_i_i_C(2) + pkin(5);
	t11 = t21 * t20 - t24 * t30;
	t19 = sin(qJ(5));
	t23 = cos(qJ(5));
	t28 = -t11 * t19 + t23 * t36;
	t4 = t11 * t23 + t19 * t36;
	t27 = t18 * r_i_i_C(1) + t22 * r_i_i_C(2) + t33;
	t26 = t40 * t19 + t29 * t23 + t35;
	t14 = -t20 * t31 + t25 * t24;
	t13 = t25 * t20 + t24 * t31;
	t12 = t20 * t30 + t21 * t24;
	t10 = t34 * t19 + t23 * t38;
	t8 = t13 * t23 - t19 * t37;
	t7 = t13 * t19 + t23 * t37;
	t2 = t14 * t18 + t8 * t22;
	t1 = t14 * t22 - t8 * t18;
	t3 = [-t21 * pkin(1) - t35 * t11 - t27 * t12 + t25 * t32 + t40 * t28 - t29 * t4, -t27 * t13 + t26 * t14, t13, -t37, -t29 * t7 + t40 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t25 * pkin(1) + t8 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t35 * t13 + t33 * t14 + t21 * t32 + t40 * t7, -t27 * t11 + t26 * t12, t11, t36, t29 * t28 + t40 * t4, (t12 * t22 - t4 * t18) * r_i_i_C(1) + (-t12 * t18 - t4 * t22) * r_i_i_C(2); 0, (t26 * t20 + t27 * t24) * t17, -t38, -t34, -t40 * t10 + t29 * (t19 * t38 - t34 * t23), (t10 * t18 + t22 * t39) * r_i_i_C(1) + (t10 * t22 - t18 * t39) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end