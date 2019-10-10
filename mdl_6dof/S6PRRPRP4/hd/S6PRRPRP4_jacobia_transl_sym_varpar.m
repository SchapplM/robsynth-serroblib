% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:21
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(10));
	t1 = sin(pkin(10));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (35->19), mult. (92->37), div. (0->0), fcn. (112->8), ass. (0->19)
	t6 = sin(pkin(6));
	t9 = sin(qJ(3));
	t20 = t6 * t9;
	t19 = pkin(8) + r_i_i_C(3);
	t11 = cos(qJ(3));
	t18 = t11 * t6;
	t10 = sin(qJ(2));
	t5 = sin(pkin(10));
	t17 = t5 * t10;
	t12 = cos(qJ(2));
	t16 = t5 * t12;
	t7 = cos(pkin(10));
	t15 = t7 * t10;
	t14 = t7 * t12;
	t13 = t11 * r_i_i_C(1) - t9 * r_i_i_C(2) + pkin(2);
	t8 = cos(pkin(6));
	t4 = -t8 * t17 + t14;
	t2 = t8 * t15 + t16;
	t1 = [0, t19 * t4 + t13 * (-t8 * t16 - t15), (t5 * t18 - t4 * t9) * r_i_i_C(1) + (-t11 * t4 - t5 * t20) * r_i_i_C(2), 0, 0, 0; 0, t19 * t2 + t13 * (t8 * t14 - t17), (-t7 * t18 - t2 * t9) * r_i_i_C(1) + (-t11 * t2 + t7 * t20) * r_i_i_C(2), 0, 0, 0; 1, (t19 * t10 + t13 * t12) * t6, (-t10 * t20 + t11 * t8) * r_i_i_C(1) + (-t10 * t18 - t8 * t9) * r_i_i_C(2), 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (66->21), mult. (169->39), div. (0->0), fcn. (211->8), ass. (0->22)
	t26 = pkin(3) - r_i_i_C(2);
	t25 = pkin(8) + r_i_i_C(1);
	t12 = sin(pkin(6));
	t15 = sin(qJ(3));
	t24 = t12 * t15;
	t17 = cos(qJ(3));
	t23 = t12 * t17;
	t14 = cos(pkin(6));
	t16 = sin(qJ(2));
	t22 = t14 * t16;
	t18 = cos(qJ(2));
	t21 = t14 * t18;
	t20 = r_i_i_C(3) + qJ(4);
	t19 = t20 * t15 + t17 * t26 + pkin(2);
	t13 = cos(pkin(10));
	t11 = sin(pkin(10));
	t9 = -t14 * t17 + t16 * t24;
	t8 = -t11 * t22 + t13 * t18;
	t6 = t11 * t18 + t13 * t22;
	t3 = -t11 * t23 + t15 * t8;
	t1 = t13 * t23 + t15 * t6;
	t2 = [0, t25 * t8 + t19 * (-t11 * t21 - t13 * t16), t20 * (t11 * t24 + t17 * t8) - t26 * t3, t3, 0, 0; 0, t25 * t6 + t19 * (-t11 * t16 + t13 * t21), t20 * (-t13 * t24 + t17 * t6) - t26 * t1, t1, 0, 0; 1, (t16 * t25 + t18 * t19) * t12, -t26 * t9 + t20 * (t14 * t15 + t16 * t23), t9, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (125->33), mult. (325->63), div. (0->0), fcn. (411->10), ass. (0->28)
	t17 = sin(qJ(3));
	t20 = cos(qJ(3));
	t16 = sin(qJ(5));
	t19 = cos(qJ(5));
	t24 = t16 * r_i_i_C(1) + t19 * r_i_i_C(2) + qJ(4);
	t28 = pkin(3) + pkin(9) + r_i_i_C(3);
	t34 = t24 * t17 + t28 * t20 + pkin(2);
	t15 = sin(pkin(6));
	t33 = t15 * t17;
	t32 = t15 * t20;
	t21 = cos(qJ(2));
	t31 = t15 * t21;
	t30 = cos(pkin(6));
	t29 = cos(pkin(10));
	t14 = sin(pkin(10));
	t27 = t14 * t30;
	t26 = t15 * t29;
	t25 = t30 * t29;
	t23 = t19 * r_i_i_C(1) - t16 * r_i_i_C(2) + pkin(4) + pkin(8);
	t18 = sin(qJ(2));
	t9 = t18 * t33 - t30 * t20;
	t8 = -t18 * t27 + t29 * t21;
	t7 = t29 * t18 + t21 * t27;
	t6 = t14 * t21 + t18 * t25;
	t5 = t14 * t18 - t21 * t25;
	t3 = -t14 * t32 + t8 * t17;
	t1 = t6 * t17 + t20 * t26;
	t2 = [0, t23 * t8 - t34 * t7, t24 * (t14 * t33 + t8 * t20) - t28 * t3, t3, (-t7 * t16 + t3 * t19) * r_i_i_C(1) + (-t3 * t16 - t7 * t19) * r_i_i_C(2), 0; 0, t23 * t6 - t34 * t5, t24 * (-t17 * t26 + t6 * t20) - t28 * t1, t1, (t1 * t19 - t5 * t16) * r_i_i_C(1) + (-t1 * t16 - t5 * t19) * r_i_i_C(2), 0; 1, (t23 * t18 + t34 * t21) * t15, -t28 * t9 + t24 * (t30 * t17 + t18 * t32), t9, (t16 * t31 + t9 * t19) * r_i_i_C(1) + (-t9 * t16 + t19 * t31) * r_i_i_C(2), 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (165->34), mult. (414->63), div. (0->0), fcn. (522->10), ass. (0->32)
	t40 = pkin(5) + r_i_i_C(1);
	t19 = sin(qJ(3));
	t22 = cos(qJ(3));
	t18 = sin(qJ(5));
	t21 = cos(qJ(5));
	t25 = t21 * r_i_i_C(2) + t40 * t18 + qJ(4);
	t34 = pkin(3) + r_i_i_C(3) + qJ(6) + pkin(9);
	t24 = -t25 * t19 - t34 * t22 - pkin(2);
	t16 = sin(pkin(6));
	t39 = t16 * t19;
	t38 = t16 * t22;
	t23 = cos(qJ(2));
	t37 = t16 * t23;
	t36 = cos(pkin(6));
	t35 = cos(pkin(10));
	t15 = sin(pkin(10));
	t32 = t15 * t36;
	t31 = t16 * t35;
	t28 = t36 * t35;
	t26 = -t18 * r_i_i_C(2) + t40 * t21 + pkin(4) + pkin(8);
	t20 = sin(qJ(2));
	t10 = t36 * t19 + t20 * t38;
	t9 = t20 * t39 - t36 * t22;
	t8 = -t20 * t32 + t35 * t23;
	t7 = t35 * t20 + t23 * t32;
	t6 = t15 * t23 + t20 * t28;
	t5 = t15 * t20 - t23 * t28;
	t4 = t15 * t39 + t8 * t22;
	t3 = -t15 * t38 + t8 * t19;
	t2 = -t19 * t31 + t6 * t22;
	t1 = t6 * t19 + t22 * t31;
	t11 = [0, t24 * t7 + t26 * t8, t25 * t4 - t34 * t3, t3, (-t3 * t18 - t7 * t21) * r_i_i_C(2) + t40 * (-t7 * t18 + t3 * t21), t4; 0, t24 * t5 + t26 * t6, -t34 * t1 + t25 * t2, t1, (-t1 * t18 - t5 * t21) * r_i_i_C(2) + t40 * (t1 * t21 - t5 * t18), t2; 1, (t26 * t20 - t24 * t23) * t16, t25 * t10 - t34 * t9, t9, (-t9 * t18 + t21 * t37) * r_i_i_C(2) + t40 * (t18 * t37 + t9 * t21), t10;];
	Ja_transl = t11;
else
	Ja_transl=NaN(3,6);
end