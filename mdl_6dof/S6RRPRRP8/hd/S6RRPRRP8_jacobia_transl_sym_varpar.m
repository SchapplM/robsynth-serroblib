% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:41
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
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
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
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
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (29->11), mult. (65->16), div. (0->0), fcn. (72->6), ass. (0->13)
	t1 = sin(pkin(10));
	t2 = cos(pkin(10));
	t10 = r_i_i_C(1) * t2 - r_i_i_C(2) * t1 + pkin(2);
	t11 = r_i_i_C(3) + qJ(3);
	t3 = sin(qJ(2));
	t5 = cos(qJ(2));
	t8 = t10 * t5 + t11 * t3;
	t12 = pkin(1) + t8;
	t9 = t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + pkin(7);
	t7 = -t10 * t3 + t11 * t5;
	t6 = cos(qJ(1));
	t4 = sin(qJ(1));
	t13 = [-t12 * t4 + t9 * t6, t7 * t6, t6 * t3, 0, 0, 0; t12 * t6 + t9 * t4, t7 * t4, t4 * t3, 0, 0, 0; 0, t8, -t5, 0, 0, 0;];
	Ja_transl = t13;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (72->24), mult. (94->38), div. (0->0), fcn. (105->8), ass. (0->21)
	t13 = cos(qJ(2));
	t11 = sin(qJ(2));
	t22 = r_i_i_C(3) + pkin(8) + qJ(3);
	t17 = t22 * t11;
	t5 = cos(pkin(10)) * pkin(3) + pkin(2);
	t23 = t13 * t5 + pkin(1) + t17;
	t12 = sin(qJ(1));
	t21 = t12 * t13;
	t14 = cos(qJ(1));
	t20 = t13 * t14;
	t19 = sin(pkin(10)) * pkin(3) + pkin(7);
	t8 = pkin(10) + qJ(4);
	t6 = sin(t8);
	t7 = cos(t8);
	t16 = r_i_i_C(1) * t7 - r_i_i_C(2) * t6 + t5;
	t15 = -t16 * t11 + t22 * t13;
	t4 = t12 * t6 + t7 * t20;
	t3 = t12 * t7 - t6 * t20;
	t2 = t14 * t6 - t7 * t21;
	t1 = t14 * t7 + t6 * t21;
	t9 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t23 * t12 + t19 * t14, t15 * t14, t14 * t11, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t19 * t12 + t23 * t14, t15 * t12, t12 * t11, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t16 * t13 + t17, -t13, (-r_i_i_C(1) * t6 - r_i_i_C(2) * t7) * t11, 0, 0;];
	Ja_transl = t9;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (152->32), mult. (133->47), div. (0->0), fcn. (148->10), ass. (0->29)
	t20 = cos(qJ(2));
	t18 = sin(qJ(2));
	t30 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3);
	t25 = t30 * t18;
	t17 = pkin(10) + qJ(4);
	t14 = cos(t17);
	t9 = pkin(4) * t14 + cos(pkin(10)) * pkin(3) + pkin(2);
	t35 = t20 * t9 + pkin(1) + t25;
	t15 = qJ(5) + t17;
	t11 = sin(t15);
	t12 = cos(t15);
	t21 = cos(qJ(1));
	t27 = t21 * t12;
	t19 = sin(qJ(1));
	t29 = t19 * t20;
	t5 = t11 * t29 + t27;
	t28 = t21 * t11;
	t6 = -t12 * t29 + t28;
	t34 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t7 = t19 * t12 - t20 * t28;
	t8 = t19 * t11 + t20 * t27;
	t33 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t13 = sin(t17);
	t32 = pkin(4) * t13;
	t31 = pkin(7) + t32 + sin(pkin(10)) * pkin(3);
	t24 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
	t23 = r_i_i_C(1) * t12 - r_i_i_C(2) * t11 + t9;
	t22 = -t23 * t18 + t30 * t20;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t35 * t19 + t31 * t21, t22 * t21, t21 * t18, (-t13 * t20 * t21 + t14 * t19) * pkin(4) + t33, t33, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t31 * t19 + t35 * t21, t22 * t19, t19 * t18, (-t13 * t29 - t14 * t21) * pkin(4) + t34, t34, 0; 0, t23 * t20 + t25, -t20, (t24 - t32) * t18, t24 * t18, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:41:03
	% EndTime: 2019-10-10 10:41:03
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (266->35), mult. (212->48), div. (0->0), fcn. (244->10), ass. (0->32)
	t43 = pkin(5) + r_i_i_C(1);
	t35 = r_i_i_C(3) + qJ(6);
	t23 = pkin(10) + qJ(4);
	t20 = cos(t23);
	t11 = pkin(4) * t20 + cos(pkin(10)) * pkin(3) + pkin(2);
	t26 = cos(qJ(2));
	t24 = sin(qJ(2));
	t40 = r_i_i_C(2) + pkin(9) + pkin(8) + qJ(3);
	t32 = t40 * t24;
	t45 = t26 * t11 + pkin(1) + t32;
	t21 = qJ(5) + t23;
	t17 = sin(t21);
	t18 = cos(t21);
	t44 = t35 * t17 + t43 * t18 + t11;
	t19 = sin(t23);
	t42 = pkin(4) * t19;
	t41 = pkin(7) + t42 + sin(pkin(10)) * pkin(3);
	t25 = sin(qJ(1));
	t38 = t25 * t26;
	t27 = cos(qJ(1));
	t37 = t27 * t17;
	t36 = t27 * t18;
	t34 = t35 * t18 * t24;
	t33 = t43 * t17;
	t7 = t17 * t38 + t36;
	t8 = t18 * t38 - t37;
	t30 = t35 * t8 - t43 * t7;
	t10 = t25 * t17 + t26 * t36;
	t9 = -t25 * t18 + t26 * t37;
	t29 = t35 * t10 - t43 * t9;
	t28 = -t44 * t24 + t40 * t26;
	t1 = [-t45 * t25 + t41 * t27 - t35 * t7 - t43 * t8, t28 * t27, t27 * t24, (-t19 * t26 * t27 + t20 * t25) * pkin(4) + t29, t29, t9; t43 * t10 + t41 * t25 + t45 * t27 + t35 * t9, t28 * t25, t25 * t24, (-t19 * t38 - t20 * t27) * pkin(4) + t30, t30, t7; 0, t44 * t26 + t32, -t26, (-t33 - t42) * t24 + t34, -t24 * t33 + t34, t24 * t17;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end