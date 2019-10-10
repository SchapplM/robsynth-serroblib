% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR9
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
% Datum: 2019-10-10 09:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR9_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR9_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR9_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
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
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
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
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
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
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (53->18), mult. (118->28), div. (0->0), fcn. (148->6), ass. (0->16)
	t17 = r_i_i_C(2) + qJ(3);
	t16 = cos(pkin(6));
	t15 = pkin(2) + r_i_i_C(3) + qJ(4);
	t9 = sin(qJ(1));
	t14 = t9 * t16;
	t11 = cos(qJ(1));
	t13 = t11 * t16;
	t7 = sin(pkin(6));
	t12 = (pkin(3) + pkin(8) + r_i_i_C(1)) * t7;
	t10 = cos(qJ(2));
	t8 = sin(qJ(2));
	t4 = t10 * t11 - t8 * t14;
	t3 = t10 * t14 + t11 * t8;
	t2 = t9 * t10 + t8 * t13;
	t1 = -t10 * t13 + t8 * t9;
	t5 = [-t9 * pkin(1) - t17 * t1 + t11 * t12 - t15 * t2, -t15 * t3 + t17 * t4, t3, t4, 0, 0; pkin(1) * t11 + t9 * t12 + t15 * t4 + t17 * t3, -t15 * t1 + t17 * t2, t1, t2, 0, 0; 0, (t15 * t10 + t17 * t8) * t7, -t7 * t10, t7 * t8, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (88->32), mult. (208->54), div. (0->0), fcn. (260->8), ass. (0->26)
	t10 = sin(pkin(6));
	t12 = sin(qJ(2));
	t27 = t10 * t12;
	t13 = sin(qJ(1));
	t26 = t10 * t13;
	t14 = cos(qJ(5));
	t25 = t10 * t14;
	t16 = cos(qJ(1));
	t24 = t10 * t16;
	t23 = pkin(2) + qJ(4);
	t22 = cos(pkin(6));
	t21 = pkin(3) + pkin(4) + pkin(8);
	t20 = r_i_i_C(3) + pkin(9) - qJ(3);
	t19 = t13 * t22;
	t18 = t16 * t22;
	t11 = sin(qJ(5));
	t17 = t11 * r_i_i_C(1) + t14 * r_i_i_C(2) + t23;
	t15 = cos(qJ(2));
	t8 = t14 * t24;
	t6 = -t12 * t19 + t16 * t15;
	t5 = t16 * t12 + t15 * t19;
	t4 = t12 * t18 + t13 * t15;
	t3 = t13 * t12 - t15 * t18;
	t2 = t6 * t11 + t13 * t25;
	t1 = -t11 * t26 + t6 * t14;
	t7 = [-t13 * pkin(1) + t8 * r_i_i_C(1) - t17 * t4 + t20 * t3 + (-t11 * r_i_i_C(2) + t21) * t24, -t17 * t5 - t20 * t6, t5, t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t16 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t20 * t5 + t21 * t26 + t23 * t6, -t17 * t3 - t20 * t4, t3, t4, (t11 * t24 + t4 * t14) * r_i_i_C(1) + (-t4 * t11 + t8) * r_i_i_C(2), 0; 0, (-t20 * t12 + t17 * t15) * t10, -t10 * t15, t27, (-t22 * t11 + t12 * t25) * r_i_i_C(1) + (-t11 * t27 - t22 * t14) * r_i_i_C(2), 0;];
	Ja_transl = t7;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (177->49), mult. (439->79), div. (0->0), fcn. (558->10), ass. (0->35)
	t21 = sin(qJ(5));
	t25 = cos(qJ(5));
	t20 = sin(qJ(6));
	t24 = cos(qJ(6));
	t31 = t24 * r_i_i_C(1) - t20 * r_i_i_C(2) + pkin(5);
	t37 = pkin(2) + qJ(4);
	t42 = pkin(10) + r_i_i_C(3);
	t43 = t31 * t21 - t42 * t25 + t37;
	t19 = sin(pkin(6));
	t22 = sin(qJ(2));
	t41 = t19 * t22;
	t40 = t19 * t25;
	t26 = cos(qJ(2));
	t39 = t19 * t26;
	t27 = cos(qJ(1));
	t38 = t19 * t27;
	t36 = -pkin(9) + qJ(3);
	t35 = cos(pkin(6));
	t23 = sin(qJ(1));
	t34 = t23 * t35;
	t33 = t27 * t35;
	t32 = t19 * (pkin(3) + pkin(4) + pkin(8));
	t14 = t22 * t33 + t23 * t26;
	t7 = -t14 * t21 + t25 * t38;
	t30 = t14 * t25 + t21 * t38;
	t29 = t20 * r_i_i_C(1) + t24 * r_i_i_C(2) - t36;
	t16 = -t22 * t34 + t27 * t26;
	t15 = t27 * t22 + t26 * t34;
	t13 = t23 * t22 - t26 * t33;
	t12 = t21 * t41 + t35 * t25;
	t4 = t16 * t21 + t23 * t40;
	t3 = t23 * t19 * t21 - t16 * t25;
	t2 = -t15 * t20 + t4 * t24;
	t1 = -t15 * t24 - t4 * t20;
	t5 = [-t23 * pkin(1) + t29 * t13 - t37 * t14 + t27 * t32 + t42 * t30 + t31 * t7, -t15 * t43 - t29 * t16, t15, t16, -t31 * t3 + t42 * t4, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t27 * pkin(1) + t4 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t36 * t15 + t37 * t16 + t23 * t32 + t42 * t3, -t13 * t43 - t29 * t14, t13, t14, t31 * t30 - t42 * t7, (-t13 * t24 + t7 * t20) * r_i_i_C(1) + (t13 * t20 + t7 * t24) * r_i_i_C(2); 0, (-t29 * t22 + t43 * t26) * t19, -t39, t41, t42 * t12 + t31 * (-t35 * t21 + t22 * t40), (-t12 * t20 + t24 * t39) * r_i_i_C(1) + (-t12 * t24 - t20 * t39) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end