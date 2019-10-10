% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRR5
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
% Datum: 2019-10-09 22:01
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
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
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
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
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.16s
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
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (102->28), mult. (179->51), div. (0->0), fcn. (220->10), ass. (0->25)
	t14 = qJ(4) + qJ(5);
	t12 = sin(t14);
	t13 = cos(t14);
	t15 = sin(pkin(11));
	t16 = sin(pkin(6));
	t31 = t15 * t16;
	t17 = cos(pkin(11));
	t20 = sin(qJ(2));
	t18 = cos(pkin(6));
	t22 = cos(qJ(2));
	t26 = t18 * t22;
	t9 = t15 * t26 + t17 * t20;
	t34 = (-t12 * t31 + t13 * t9) * r_i_i_C(1) + (-t12 * t9 - t13 * t31) * r_i_i_C(2);
	t30 = t16 * t17;
	t7 = t15 * t20 - t17 * t26;
	t33 = (t12 * t30 + t13 * t7) * r_i_i_C(1) + (-t12 * t7 + t13 * t30) * r_i_i_C(2);
	t28 = t16 * t22;
	t32 = (-t12 * t18 - t13 * t28) * r_i_i_C(1) + (t12 * t28 - t13 * t18) * r_i_i_C(2);
	t19 = sin(qJ(4));
	t29 = t16 * t19;
	t27 = t18 * t20;
	t25 = pkin(2) + r_i_i_C(3) + pkin(9) + pkin(8);
	t24 = pkin(4) * t19 + r_i_i_C(1) * t12 + r_i_i_C(2) * t13 + qJ(3);
	t21 = cos(qJ(4));
	t1 = [0, -t25 * t9 + t24 * (-t15 * t27 + t17 * t22), t9, (-t15 * t29 + t21 * t9) * pkin(4) + t34, t34, 0; 0, -t25 * t7 + t24 * (t15 * t22 + t17 * t27), t7, (t17 * t29 + t21 * t7) * pkin(4) + t33, t33, 0; 1, (t24 * t20 + t25 * t22) * t16, -t28, (-t18 * t19 - t21 * t28) * pkin(4) + t32, t32, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (250->43), mult. (416->74), div. (0->0), fcn. (521->12), ass. (0->35)
	t52 = pkin(10) + r_i_i_C(3);
	t30 = sin(qJ(6));
	t33 = cos(qJ(6));
	t51 = t33 * r_i_i_C(1) - t30 * r_i_i_C(2) + pkin(5);
	t26 = sin(pkin(11));
	t27 = sin(pkin(6));
	t48 = t26 * t27;
	t28 = cos(pkin(11));
	t47 = t27 * t28;
	t31 = sin(qJ(4));
	t46 = t27 * t31;
	t32 = sin(qJ(2));
	t45 = t27 * t32;
	t35 = cos(qJ(2));
	t44 = t27 * t35;
	t29 = cos(pkin(6));
	t43 = t29 * t32;
	t42 = t29 * t35;
	t41 = t30 * r_i_i_C(1) + t33 * r_i_i_C(2) + pkin(2) + pkin(8) + pkin(9);
	t19 = t26 * t42 + t28 * t32;
	t25 = qJ(4) + qJ(5);
	t23 = sin(t25);
	t24 = cos(t25);
	t8 = t19 * t23 + t24 * t48;
	t40 = t52 * t8 + t51 * (t19 * t24 - t23 * t48);
	t17 = t26 * t32 - t28 * t42;
	t10 = -t17 * t23 + t24 * t47;
	t39 = -t52 * t10 + t51 * (t17 * t24 + t23 * t47);
	t16 = -t23 * t44 + t29 * t24;
	t38 = t52 * t16 + t51 * (-t29 * t23 - t24 * t44);
	t37 = t31 * pkin(4) + t51 * t23 - t52 * t24 + qJ(3);
	t34 = cos(qJ(4));
	t20 = -t26 * t43 + t28 * t35;
	t18 = t26 * t35 + t28 * t43;
	t1 = [0, -t41 * t19 + t37 * t20, t19, (t19 * t34 - t26 * t46) * pkin(4) + t40, t40, (t20 * t33 - t8 * t30) * r_i_i_C(1) + (-t20 * t30 - t8 * t33) * r_i_i_C(2); 0, -t41 * t17 + t37 * t18, t17, (t17 * t34 + t28 * t46) * pkin(4) + t39, t39, (t10 * t30 + t18 * t33) * r_i_i_C(1) + (t10 * t33 - t18 * t30) * r_i_i_C(2); 1, (t37 * t32 + t41 * t35) * t27, -t44, (-t29 * t31 - t34 * t44) * pkin(4) + t38, t38, (-t16 * t30 + t33 * t45) * r_i_i_C(1) + (-t16 * t33 - t30 * t45) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end