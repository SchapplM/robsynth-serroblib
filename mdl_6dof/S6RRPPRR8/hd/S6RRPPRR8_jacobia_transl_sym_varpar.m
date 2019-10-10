% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:48
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
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
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
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
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (45->20), mult. (104->30), div. (0->0), fcn. (120->6), ass. (0->18)
	t15 = r_i_i_C(2) + qJ(3);
	t7 = sin(qJ(2));
	t12 = t15 * t7;
	t9 = cos(qJ(2));
	t20 = t9 * pkin(2) + pkin(1) + t12;
	t14 = r_i_i_C(3) + qJ(4);
	t17 = pkin(3) + r_i_i_C(1);
	t5 = sin(pkin(10));
	t6 = cos(pkin(10));
	t19 = t14 * t5 + t17 * t6 + pkin(2);
	t8 = sin(qJ(1));
	t18 = t8 * t9;
	t10 = cos(qJ(1));
	t16 = t10 * t9;
	t11 = t15 * t9 - t19 * t7;
	t3 = t5 * t16 - t8 * t6;
	t1 = t10 * t6 + t5 * t18;
	t2 = [pkin(7) * t10 + t17 * (t10 * t5 - t6 * t18) - t14 * t1 - t20 * t8, t11 * t10, t10 * t7, t3, 0, 0; t8 * pkin(7) + t17 * (t6 * t16 + t8 * t5) + t14 * t3 + t20 * t10, t11 * t8, t8 * t7, t1, 0, 0; 0, t19 * t9 + t12, -t9, t7 * t5, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (88->34), mult. (217->55), div. (0->0), fcn. (263->8), ass. (0->25)
	t14 = cos(qJ(2));
	t11 = sin(qJ(2));
	t21 = -r_i_i_C(3) - pkin(8) + qJ(3);
	t19 = t21 * t11;
	t26 = t14 * pkin(2) + pkin(1) + t19;
	t10 = sin(qJ(5));
	t13 = cos(qJ(5));
	t24 = pkin(3) + pkin(4);
	t17 = t13 * r_i_i_C(1) - t10 * r_i_i_C(2) + t24;
	t18 = t10 * r_i_i_C(1) + t13 * r_i_i_C(2) + qJ(4);
	t8 = sin(pkin(10));
	t9 = cos(pkin(10));
	t25 = t17 * t9 + t18 * t8 + pkin(2);
	t12 = sin(qJ(1));
	t23 = t12 * t14;
	t15 = cos(qJ(1));
	t22 = t14 * t15;
	t16 = -t25 * t11 + t21 * t14;
	t6 = t12 * t8 + t9 * t22;
	t5 = -t12 * t9 + t8 * t22;
	t4 = -t15 * t8 + t9 * t23;
	t3 = t15 * t9 + t8 * t23;
	t2 = t5 * t10 + t13 * t6;
	t1 = -t6 * t10 + t13 * t5;
	t7 = [t15 * pkin(7) - t26 * t12 - t17 * t4 - t18 * t3, t16 * t15, t15 * t11, t5, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t12 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * qJ(4) + t26 * t15 + t24 * t6, t16 * t12, t12 * t11, t3, (-t4 * t10 + t3 * t13) * r_i_i_C(1) + (-t3 * t10 - t4 * t13) * r_i_i_C(2), 0; 0, t25 * t14 + t19, -t14, t11 * t8, ((-t10 * t9 + t13 * t8) * r_i_i_C(1) + (-t10 * t8 - t13 * t9) * r_i_i_C(2)) * t11, 0;];
	Ja_transl = t7;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (176->43), mult. (315->66), div. (0->0), fcn. (381->10), ass. (0->33)
	t24 = cos(qJ(2));
	t21 = sin(qJ(2));
	t33 = r_i_i_C(3) - qJ(3) + pkin(9) + pkin(8);
	t30 = t33 * t21;
	t42 = -t24 * pkin(2) - pkin(1) + t30;
	t18 = sin(pkin(10));
	t19 = cos(pkin(10));
	t17 = qJ(5) + qJ(6);
	t15 = sin(t17);
	t16 = cos(t17);
	t20 = sin(qJ(5));
	t31 = t20 * pkin(5) + qJ(4);
	t28 = t15 * r_i_i_C(1) + t16 * r_i_i_C(2) + t31;
	t23 = cos(qJ(5));
	t37 = t23 * pkin(5) + pkin(3) + pkin(4);
	t29 = t16 * r_i_i_C(1) - t15 * r_i_i_C(2) + t37;
	t41 = t28 * t18 + t29 * t19 + pkin(2);
	t25 = cos(qJ(1));
	t35 = t25 * t18;
	t22 = sin(qJ(1));
	t36 = t22 * t24;
	t10 = t19 * t36 - t35;
	t34 = t25 * t19;
	t9 = t18 * t36 + t34;
	t40 = (-t10 * t15 + t9 * t16) * r_i_i_C(1) + (-t10 * t16 - t9 * t15) * r_i_i_C(2);
	t11 = -t22 * t19 + t24 * t35;
	t12 = t22 * t18 + t24 * t34;
	t5 = t11 * t16 - t12 * t15;
	t6 = t11 * t15 + t12 * t16;
	t39 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t38 = ((-t15 * t19 + t16 * t18) * r_i_i_C(1) + (-t15 * t18 - t16 * t19) * r_i_i_C(2)) * t21;
	t27 = -t41 * t21 - t33 * t24;
	t1 = [t25 * pkin(7) - t29 * t10 + t42 * t22 - t28 * t9, t27 * t25, t25 * t21, t11, (t11 * t23 - t12 * t20) * pkin(5) + t39, t39; t22 * pkin(7) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t31 * t11 + t37 * t12 - t42 * t25, t27 * t22, t22 * t21, t9, (-t10 * t20 + t23 * t9) * pkin(5) + t40, t40; 0, t41 * t24 - t30, -t24, t21 * t18, (t18 * t23 - t19 * t20) * t21 * pkin(5) + t38, t38;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end