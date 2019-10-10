% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:59
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:26
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:26
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:26
	% EndTime: 2019-10-09 21:59:27
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
	% StartTime: 2019-10-09 21:59:26
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (25->11), mult. (63->20), div. (0->0), fcn. (78->8), ass. (0->13)
	t11 = cos(pkin(6));
	t12 = sin(qJ(2));
	t17 = t11 * t12;
	t13 = cos(qJ(2));
	t16 = t11 * t13;
	t15 = r_i_i_C(3) + qJ(3);
	t14 = r_i_i_C(1) * cos(pkin(12)) - r_i_i_C(2) * sin(pkin(12)) + pkin(2);
	t10 = cos(pkin(11));
	t8 = sin(pkin(6));
	t7 = sin(pkin(11));
	t3 = t10 * t12 + t7 * t16;
	t1 = -t10 * t16 + t7 * t12;
	t2 = [0, t15 * (t10 * t13 - t7 * t17) - t14 * t3, t3, 0, 0, 0; 0, t15 * (t10 * t17 + t7 * t13) - t14 * t1, t1, 0, 0, 0; 1, (t15 * t12 + t14 * t13) * t8, -t8 * t13, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:26
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (62->22), mult. (102->42), div. (0->0), fcn. (127->9), ass. (0->21)
	t23 = r_i_i_C(3) + pkin(8) + qJ(3);
	t10 = sin(pkin(11));
	t11 = sin(pkin(6));
	t22 = t10 * t11;
	t12 = cos(pkin(11));
	t21 = t11 * t12;
	t15 = sin(qJ(2));
	t20 = t11 * t15;
	t13 = cos(pkin(6));
	t19 = t13 * t15;
	t16 = cos(qJ(2));
	t18 = t13 * t16;
	t9 = pkin(12) + qJ(4);
	t7 = sin(t9);
	t8 = cos(t9);
	t17 = r_i_i_C(1) * t8 - r_i_i_C(2) * t7 + cos(pkin(12)) * pkin(3) + pkin(2);
	t4 = -t10 * t19 + t12 * t16;
	t3 = t10 * t18 + t12 * t15;
	t2 = t10 * t16 + t12 * t19;
	t1 = t10 * t15 - t12 * t18;
	t5 = [0, -t17 * t3 + t23 * t4, t3, (t8 * t22 - t4 * t7) * r_i_i_C(1) + (-t7 * t22 - t4 * t8) * r_i_i_C(2), 0, 0; 0, -t17 * t1 + t23 * t2, t1, (-t2 * t7 - t8 * t21) * r_i_i_C(1) + (-t2 * t8 + t7 * t21) * r_i_i_C(2), 0, 0; 1, (t23 * t15 + t17 * t16) * t11, -t11 * t16, (t13 * t8 - t7 * t20) * r_i_i_C(1) + (-t13 * t7 - t8 * t20) * r_i_i_C(2), 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:26
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (168->36), mult. (279->64), div. (0->0), fcn. (353->11), ass. (0->29)
	t15 = pkin(12) + qJ(4);
	t13 = sin(t15);
	t14 = cos(t15);
	t19 = sin(qJ(5));
	t21 = cos(qJ(5));
	t25 = t21 * r_i_i_C(1) - t19 * r_i_i_C(2) + pkin(4);
	t34 = pkin(9) + r_i_i_C(3);
	t35 = t34 * t13 + t25 * t14 + cos(pkin(12)) * pkin(3) + pkin(2);
	t16 = sin(pkin(11));
	t17 = sin(pkin(6));
	t33 = t16 * t17;
	t20 = sin(qJ(2));
	t32 = t17 * t20;
	t22 = cos(qJ(2));
	t31 = t17 * t22;
	t30 = cos(pkin(6));
	t29 = cos(pkin(11));
	t28 = t16 * t30;
	t27 = t17 * t29;
	t26 = t30 * t29;
	t24 = t19 * r_i_i_C(1) + t21 * r_i_i_C(2) + pkin(8) + qJ(3);
	t10 = -t20 * t28 + t29 * t22;
	t9 = t29 * t20 + t22 * t28;
	t8 = t16 * t22 + t20 * t26;
	t7 = t16 * t20 - t22 * t26;
	t6 = t30 * t13 + t14 * t32;
	t4 = t10 * t14 + t13 * t33;
	t2 = -t13 * t27 + t8 * t14;
	t1 = [0, t24 * t10 - t35 * t9, t9, t34 * t4 + t25 * (-t10 * t13 + t14 * t33), (-t4 * t19 + t9 * t21) * r_i_i_C(1) + (-t9 * t19 - t4 * t21) * r_i_i_C(2), 0; 0, t24 * t8 - t35 * t7, t7, t34 * t2 + t25 * (-t8 * t13 - t14 * t27), (-t2 * t19 + t7 * t21) * r_i_i_C(1) + (-t7 * t19 - t2 * t21) * r_i_i_C(2), 0; 1, (t24 * t20 + t35 * t22) * t17, -t31, t34 * t6 + t25 * (-t13 * t32 + t30 * t14), (-t6 * t19 - t21 * t31) * r_i_i_C(1) + (t19 * t31 - t6 * t21) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:59:26
	% EndTime: 2019-10-09 21:59:27
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (281->45), mult. (394->75), div. (0->0), fcn. (497->13), ass. (0->35)
	t24 = pkin(12) + qJ(4);
	t20 = sin(t24);
	t21 = cos(t24);
	t25 = qJ(5) + qJ(6);
	t22 = sin(t25);
	t23 = cos(t25);
	t31 = cos(qJ(5));
	t36 = t31 * pkin(5) + r_i_i_C(1) * t23 - r_i_i_C(2) * t22 + pkin(4);
	t45 = r_i_i_C(3) + pkin(10) + pkin(9);
	t49 = t45 * t20 + t36 * t21 + cos(pkin(12)) * pkin(3) + pkin(2);
	t26 = sin(pkin(11));
	t30 = sin(qJ(2));
	t32 = cos(qJ(2));
	t40 = cos(pkin(11));
	t41 = cos(pkin(6));
	t37 = t41 * t40;
	t13 = t26 * t30 - t32 * t37;
	t14 = t26 * t32 + t30 * t37;
	t27 = sin(pkin(6));
	t38 = t27 * t40;
	t8 = t14 * t21 - t20 * t38;
	t48 = (t13 * t23 - t8 * t22) * r_i_i_C(1) + (-t13 * t22 - t8 * t23) * r_i_i_C(2);
	t39 = t26 * t41;
	t16 = -t30 * t39 + t40 * t32;
	t44 = t26 * t27;
	t10 = t16 * t21 + t20 * t44;
	t15 = t40 * t30 + t32 * t39;
	t47 = (-t10 * t22 + t15 * t23) * r_i_i_C(1) + (-t10 * t23 - t15 * t22) * r_i_i_C(2);
	t43 = t27 * t30;
	t12 = t41 * t20 + t21 * t43;
	t42 = t27 * t32;
	t46 = (-t12 * t22 - t23 * t42) * r_i_i_C(1) + (-t12 * t23 + t22 * t42) * r_i_i_C(2);
	t29 = sin(qJ(5));
	t35 = t29 * pkin(5) + t22 * r_i_i_C(1) + t23 * r_i_i_C(2) + pkin(8) + qJ(3);
	t1 = [0, -t15 * t49 + t35 * t16, t15, t45 * t10 + t36 * (-t16 * t20 + t21 * t44), (-t10 * t29 + t15 * t31) * pkin(5) + t47, t47; 0, -t13 * t49 + t35 * t14, t13, t45 * t8 + t36 * (-t14 * t20 - t21 * t38), (t13 * t31 - t29 * t8) * pkin(5) + t48, t48; 1, (t35 * t30 + t32 * t49) * t27, -t42, t45 * t12 + t36 * (-t20 * t43 + t41 * t21), (-t12 * t29 - t31 * t42) * pkin(5) + t46, t46;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end