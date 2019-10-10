% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:03
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:06
	% EndTime: 2019-10-09 22:03:06
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:06
	% EndTime: 2019-10-09 22:03:06
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:06
	% EndTime: 2019-10-09 22:03:06
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(11));
	t1 = sin(pkin(11));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:06
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (20->10), mult. (47->17), div. (0->0), fcn. (60->6), ass. (0->14)
	t8 = cos(pkin(6));
	t9 = sin(qJ(2));
	t15 = t8 * t9;
	t14 = pkin(2) - r_i_i_C(2);
	t10 = cos(qJ(2));
	t5 = sin(pkin(11));
	t13 = t5 * t10;
	t7 = cos(pkin(11));
	t12 = t7 * t10;
	t11 = r_i_i_C(3) + qJ(3);
	t6 = sin(pkin(6));
	t3 = t13 * t8 + t7 * t9;
	t1 = -t12 * t8 + t5 * t9;
	t2 = [0, t11 * (-t15 * t5 + t12) - t14 * t3, t3, 0, 0, 0; 0, t11 * (t15 * t7 + t13) - t14 * t1, t1, 0, 0, 0; 1, (t10 * t14 + t11 * t9) * t6, -t6 * t10, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:06
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (43->20), mult. (109->39), div. (0->0), fcn. (136->8), ass. (0->19)
	t10 = sin(qJ(4));
	t7 = sin(pkin(6));
	t21 = t10 * t7;
	t11 = sin(qJ(2));
	t6 = sin(pkin(11));
	t20 = t11 * t6;
	t12 = cos(qJ(4));
	t19 = t12 * t7;
	t13 = cos(qJ(2));
	t9 = cos(pkin(6));
	t18 = t13 * t9;
	t17 = t7 * t13;
	t8 = cos(pkin(11));
	t16 = t8 * t11;
	t15 = pkin(2) + pkin(8) + r_i_i_C(3);
	t14 = r_i_i_C(1) * t10 + r_i_i_C(2) * t12 + qJ(3);
	t3 = t6 * t18 + t16;
	t1 = -t8 * t18 + t20;
	t2 = [0, t14 * (t13 * t8 - t9 * t20) - t15 * t3, t3, (t12 * t3 - t6 * t21) * r_i_i_C(1) + (-t10 * t3 - t6 * t19) * r_i_i_C(2), 0, 0; 0, t14 * (t13 * t6 + t9 * t16) - t15 * t1, t1, (t1 * t12 + t8 * t21) * r_i_i_C(1) + (-t1 * t10 + t8 * t19) * r_i_i_C(2), 0, 0; 1, (t14 * t11 + t15 * t13) * t7, -t17, (-t9 * t10 - t12 * t17) * r_i_i_C(1) + (t10 * t17 - t9 * t12) * r_i_i_C(2), 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:06
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (110->34), mult. (286->63), div. (0->0), fcn. (362->10), ass. (0->28)
	t32 = pkin(9) + r_i_i_C(3);
	t14 = sin(pkin(6));
	t18 = sin(qJ(4));
	t31 = t14 * t18;
	t19 = sin(qJ(2));
	t30 = t14 * t19;
	t21 = cos(qJ(4));
	t29 = t14 * t21;
	t22 = cos(qJ(2));
	t28 = t14 * t22;
	t16 = cos(pkin(6));
	t27 = t16 * t19;
	t26 = t16 * t22;
	t17 = sin(qJ(5));
	t20 = cos(qJ(5));
	t25 = t20 * r_i_i_C(1) - t17 * r_i_i_C(2) + pkin(4);
	t24 = t17 * r_i_i_C(1) + t20 * r_i_i_C(2) + pkin(2) + pkin(8);
	t23 = t25 * t18 - t32 * t21 + qJ(3);
	t15 = cos(pkin(11));
	t13 = sin(pkin(11));
	t11 = t16 * t21 - t18 * t28;
	t9 = -t13 * t27 + t15 * t22;
	t8 = t13 * t26 + t15 * t19;
	t7 = t13 * t22 + t15 * t27;
	t6 = t13 * t19 - t15 * t26;
	t4 = t15 * t29 - t6 * t18;
	t2 = t13 * t29 + t8 * t18;
	t1 = [0, t23 * t9 - t24 * t8, t8, t32 * t2 + t25 * (-t13 * t31 + t8 * t21), (-t2 * t17 + t9 * t20) * r_i_i_C(1) + (-t9 * t17 - t2 * t20) * r_i_i_C(2), 0; 0, t23 * t7 - t24 * t6, t6, -t32 * t4 + t25 * (t15 * t31 + t6 * t21), (t4 * t17 + t7 * t20) * r_i_i_C(1) + (-t7 * t17 + t4 * t20) * r_i_i_C(2), 0; 1, (t23 * t19 + t24 * t22) * t14, -t28, t32 * t11 + t25 * (-t16 * t18 - t21 * t28), (-t11 * t17 + t20 * t30) * r_i_i_C(1) + (-t11 * t20 - t17 * t30) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:03:06
	% EndTime: 2019-10-09 22:03:07
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (205->43), mult. (401->74), div. (0->0), fcn. (506->12), ass. (0->34)
	t23 = sin(pkin(11));
	t25 = cos(pkin(11));
	t32 = cos(qJ(2));
	t26 = cos(pkin(6));
	t29 = sin(qJ(2));
	t38 = t26 * t29;
	t15 = -t23 * t38 + t25 * t32;
	t22 = qJ(5) + qJ(6);
	t20 = sin(t22);
	t21 = cos(t22);
	t37 = t26 * t32;
	t14 = t23 * t37 + t25 * t29;
	t28 = sin(qJ(4));
	t24 = sin(pkin(6));
	t31 = cos(qJ(4));
	t40 = t24 * t31;
	t8 = t14 * t28 + t23 * t40;
	t46 = (t15 * t21 - t8 * t20) * r_i_i_C(1) + (-t15 * t20 - t8 * t21) * r_i_i_C(2);
	t12 = t23 * t29 - t25 * t37;
	t10 = -t12 * t28 + t25 * t40;
	t13 = t23 * t32 + t25 * t38;
	t45 = (t10 * t20 + t13 * t21) * r_i_i_C(1) + (t10 * t21 - t13 * t20) * r_i_i_C(2);
	t39 = t24 * t32;
	t17 = t26 * t31 - t28 * t39;
	t41 = t24 * t29;
	t44 = (-t17 * t20 + t21 * t41) * r_i_i_C(1) + (-t17 * t21 - t20 * t41) * r_i_i_C(2);
	t43 = r_i_i_C(3) + pkin(10) + pkin(9);
	t42 = t24 * t28;
	t30 = cos(qJ(5));
	t36 = t30 * pkin(5) + r_i_i_C(1) * t21 - r_i_i_C(2) * t20 + pkin(4);
	t27 = sin(qJ(5));
	t35 = t27 * pkin(5) + t20 * r_i_i_C(1) + t21 * r_i_i_C(2) + pkin(2) + pkin(8);
	t34 = t36 * t28 - t43 * t31 + qJ(3);
	t1 = [0, -t35 * t14 + t34 * t15, t14, t43 * t8 + t36 * (t14 * t31 - t23 * t42), (t15 * t30 - t27 * t8) * pkin(5) + t46, t46; 0, -t35 * t12 + t34 * t13, t12, -t43 * t10 + t36 * (t12 * t31 + t25 * t42), (t10 * t27 + t13 * t30) * pkin(5) + t45, t45; 1, (t34 * t29 + t35 * t32) * t24, -t39, t43 * t17 + t36 * (-t26 * t28 - t31 * t39), (-t17 * t27 + t30 * t41) * pkin(5) + t44, t44;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end