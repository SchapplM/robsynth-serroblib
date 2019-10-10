% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:53
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR11_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR11_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR11_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
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
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
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
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
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
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (64->20), mult. (150->32), div. (0->0), fcn. (188->8), ass. (0->18)
	t19 = cos(pkin(6));
	t18 = r_i_i_C(3) + qJ(4) + pkin(2);
	t11 = sin(qJ(1));
	t17 = t11 * t19;
	t13 = cos(qJ(1));
	t16 = t13 * t19;
	t7 = sin(pkin(11));
	t9 = cos(pkin(11));
	t15 = t7 * r_i_i_C(1) + t9 * r_i_i_C(2) + qJ(3);
	t8 = sin(pkin(6));
	t14 = (r_i_i_C(1) * t9 - r_i_i_C(2) * t7 + pkin(3) + pkin(8)) * t8;
	t12 = cos(qJ(2));
	t10 = sin(qJ(2));
	t4 = -t10 * t17 + t12 * t13;
	t3 = t13 * t10 + t12 * t17;
	t2 = t10 * t16 + t11 * t12;
	t1 = t10 * t11 - t12 * t16;
	t5 = [-t11 * pkin(1) - t15 * t1 + t13 * t14 - t18 * t2, t15 * t4 - t18 * t3, t3, t4, 0, 0; t13 * pkin(1) + t11 * t14 + t15 * t3 + t18 * t4, -t18 * t1 + t15 * t2, t1, t2, 0, 0; 0, (t15 * t10 + t18 * t12) * t8, -t8 * t12, t8 * t10, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (119->34), mult. (211->55), div. (0->0), fcn. (263->10), ass. (0->26)
	t30 = pkin(8) + cos(pkin(11)) * pkin(4) + pkin(3);
	t15 = sin(pkin(6));
	t18 = sin(qJ(1));
	t29 = t15 * t18;
	t19 = cos(qJ(2));
	t28 = t15 * t19;
	t20 = cos(qJ(1));
	t27 = t15 * t20;
	t26 = cos(pkin(6));
	t25 = -r_i_i_C(3) - pkin(9) - qJ(4) - pkin(2);
	t24 = sin(pkin(11)) * pkin(4) + qJ(3);
	t23 = t18 * t26;
	t22 = t20 * t26;
	t13 = pkin(11) + qJ(5);
	t11 = sin(t13);
	t12 = cos(t13);
	t21 = t11 * r_i_i_C(1) + t12 * r_i_i_C(2) + t24;
	t17 = sin(qJ(2));
	t7 = t12 * t27;
	t6 = -t17 * t23 + t20 * t19;
	t5 = t20 * t17 + t19 * t23;
	t4 = t17 * t22 + t18 * t19;
	t3 = t18 * t17 - t19 * t22;
	t2 = t5 * t11 + t12 * t29;
	t1 = -t11 * t29 + t5 * t12;
	t8 = [-t18 * pkin(1) + t7 * r_i_i_C(1) - t21 * t3 + (-t11 * r_i_i_C(2) + t30) * t27 + t25 * t4, t21 * t6 + t25 * t5, t5, t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t20 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t24 * t5 - t25 * t6 + t30 * t29, t21 * t4 + t25 * t3, t3, t4, (t11 * t27 + t3 * t12) * r_i_i_C(1) + (-t3 * t11 + t7) * r_i_i_C(2), 0; 0, (t21 * t17 - t25 * t19) * t15, -t28, t15 * t17, (-t26 * t11 - t12 * t28) * r_i_i_C(1) + (t11 * t28 - t26 * t12) * r_i_i_C(2), 0;];
	Ja_transl = t8;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (259->50), mult. (442->80), div. (0->0), fcn. (561->12), ass. (0->36)
	t46 = r_i_i_C(3) + pkin(10);
	t45 = pkin(2) + pkin(9) + qJ(4);
	t24 = sin(pkin(6));
	t27 = sin(qJ(2));
	t44 = t24 * t27;
	t28 = sin(qJ(1));
	t43 = t24 * t28;
	t30 = cos(qJ(2));
	t42 = t24 * t30;
	t31 = cos(qJ(1));
	t41 = t24 * t31;
	t40 = cos(pkin(6));
	t39 = t24 * (pkin(8) + cos(pkin(11)) * pkin(4) + pkin(3));
	t38 = sin(pkin(11)) * pkin(4) + qJ(3);
	t37 = t28 * t40;
	t36 = t31 * t40;
	t26 = sin(qJ(6));
	t29 = cos(qJ(6));
	t35 = t29 * r_i_i_C(1) - t26 * r_i_i_C(2) + pkin(5);
	t13 = t28 * t27 - t30 * t36;
	t22 = pkin(11) + qJ(5);
	t20 = sin(t22);
	t21 = cos(t22);
	t7 = -t13 * t20 + t21 * t41;
	t34 = t13 * t21 + t20 * t41;
	t33 = t26 * r_i_i_C(1) + t29 * r_i_i_C(2) + t45;
	t32 = t35 * t20 - t46 * t21 + t38;
	t16 = -t27 * t37 + t31 * t30;
	t15 = t31 * t27 + t30 * t37;
	t14 = t27 * t36 + t28 * t30;
	t12 = -t20 * t42 + t21 * t40;
	t4 = t15 * t20 + t21 * t43;
	t3 = -t15 * t21 + t20 * t43;
	t2 = t16 * t26 + t4 * t29;
	t1 = t16 * t29 - t4 * t26;
	t5 = [-t28 * pkin(1) - t38 * t13 - t33 * t14 + t31 * t39 + t46 * t34 + t35 * t7, -t15 * t33 + t16 * t32, t15, t16, -t35 * t3 + t46 * t4, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t31 * pkin(1) + t4 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t38 * t15 + t45 * t16 + t28 * t39 + t46 * t3, -t13 * t33 + t14 * t32, t13, t14, t35 * t34 - t46 * t7, (t14 * t29 + t7 * t26) * r_i_i_C(1) + (-t14 * t26 + t7 * t29) * r_i_i_C(2); 0, (t27 * t32 + t30 * t33) * t24, -t42, t44, t46 * t12 + t35 * (-t20 * t40 - t21 * t42), (-t12 * t26 + t29 * t44) * r_i_i_C(1) + (-t12 * t29 - t26 * t44) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end