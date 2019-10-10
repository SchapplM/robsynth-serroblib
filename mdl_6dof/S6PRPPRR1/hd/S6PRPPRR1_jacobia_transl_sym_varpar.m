% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:24
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPPRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (23->15), mult. (60->31), div. (0->0), fcn. (77->8), ass. (0->14)
	t10 = cos(pkin(6));
	t12 = cos(qJ(2));
	t13 = t10 * t12;
	t11 = sin(qJ(2));
	t5 = sin(pkin(11));
	t8 = cos(pkin(11));
	t4 = -t11 * t8 - t12 * t5;
	t3 = t11 * t5 - t12 * t8;
	t9 = cos(pkin(10));
	t7 = sin(pkin(6));
	t6 = sin(pkin(10));
	t2 = t4 * t10;
	t1 = t3 * t10;
	t14 = [0, (t6 * t1 + t9 * t4) * r_i_i_C(1) + (-t6 * t2 + t9 * t3) * r_i_i_C(2) + (-t9 * t11 - t6 * t13) * pkin(2), t6 * t7, 0, 0, 0; 0, (-t9 * t1 + t6 * t4) * r_i_i_C(1) + (t9 * t2 + t6 * t3) * r_i_i_C(2) + (-t6 * t11 + t9 * t13) * pkin(2), -t9 * t7, 0, 0, 0; 1, (t12 * pkin(2) - t3 * r_i_i_C(1) + t4 * r_i_i_C(2)) * t7, t10, 0, 0, 0;];
	Ja_transl = t14;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (60->20), mult. (155->34), div. (0->0), fcn. (203->10), ass. (0->19)
	t13 = sin(pkin(11));
	t17 = cos(pkin(11));
	t20 = sin(qJ(2));
	t21 = cos(qJ(2));
	t9 = t20 * t13 - t17 * t21;
	t19 = cos(pkin(6));
	t26 = t19 * t21;
	t24 = r_i_i_C(3) + qJ(4);
	t23 = t13 * t21 + t20 * t17;
	t22 = r_i_i_C(1) * cos(pkin(12)) - r_i_i_C(2) * sin(pkin(12)) + pkin(3);
	t18 = cos(pkin(10));
	t15 = sin(pkin(6));
	t14 = sin(pkin(10));
	t8 = t23 * t19;
	t7 = t9 * t19;
	t5 = t9 * t15;
	t4 = t14 * t7 - t18 * t23;
	t2 = -t14 * t23 - t18 * t7;
	t1 = [0, -t24 * (t14 * t8 + t18 * t9) + (-t14 * t26 - t18 * t20) * pkin(2) + t22 * t4, t14 * t15, -t4, 0, 0; 0, -t24 * (t14 * t9 - t18 * t8) + (-t14 * t20 + t18 * t26) * pkin(2) + t22 * t2, -t18 * t15, -t2, 0, 0; 1, -t22 * t5 + (pkin(2) * t21 + t23 * t24) * t15, t19, t5, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (107->31), mult. (220->54), div. (0->0), fcn. (288->11), ass. (0->27)
	t19 = sin(pkin(11));
	t22 = cos(pkin(11));
	t26 = sin(qJ(2));
	t27 = cos(qJ(2));
	t11 = t26 * t19 - t27 * t22;
	t36 = r_i_i_C(3) + pkin(8) + qJ(4);
	t20 = sin(pkin(10));
	t21 = sin(pkin(6));
	t35 = t20 * t21;
	t23 = cos(pkin(10));
	t34 = t23 * t21;
	t24 = cos(pkin(6));
	t33 = t24 * t27;
	t29 = t27 * t19 + t26 * t22;
	t10 = t29 * t24;
	t3 = t23 * t10 - t20 * t11;
	t30 = t20 * t10 + t23 * t11;
	t18 = pkin(12) + qJ(5);
	t16 = sin(t18);
	t17 = cos(t18);
	t28 = t17 * r_i_i_C(1) - t16 * r_i_i_C(2) + cos(pkin(12)) * pkin(4) + pkin(3);
	t9 = t11 * t24;
	t8 = t29 * t21;
	t7 = t11 * t21;
	t5 = t20 * t9 - t23 * t29;
	t2 = -t20 * t29 - t23 * t9;
	t1 = [0, -t36 * t30 + (-t20 * t33 - t23 * t26) * pkin(2) + t28 * t5, t35, -t5, (t16 * t30 + t17 * t35) * r_i_i_C(1) + (-t16 * t35 + t17 * t30) * r_i_i_C(2), 0; 0, t36 * t3 + (-t20 * t26 + t23 * t33) * pkin(2) + t28 * t2, -t34, -t2, (-t3 * t16 - t17 * t34) * r_i_i_C(1) + (t16 * t34 - t3 * t17) * r_i_i_C(2), 0; 1, t21 * t27 * pkin(2) - t28 * t7 + t36 * t8, t24, t7, (-t8 * t16 + t24 * t17) * r_i_i_C(1) + (-t24 * t16 - t8 * t17) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (268->46), mult. (540->76), div. (0->0), fcn. (712->13), ass. (0->34)
	t26 = sin(pkin(11));
	t33 = sin(qJ(2));
	t35 = cos(qJ(2));
	t43 = cos(pkin(11));
	t38 = -t33 * t26 + t35 * t43;
	t25 = pkin(12) + qJ(5);
	t23 = sin(t25);
	t24 = cos(t25);
	t32 = sin(qJ(6));
	t34 = cos(qJ(6));
	t40 = t34 * r_i_i_C(1) - t32 * r_i_i_C(2) + pkin(5);
	t48 = pkin(9) + r_i_i_C(3);
	t36 = t48 * t23 + t40 * t24 + cos(pkin(12)) * pkin(4) + pkin(3);
	t27 = sin(pkin(10));
	t28 = sin(pkin(6));
	t47 = t27 * t28;
	t29 = cos(pkin(10));
	t46 = t29 * t28;
	t30 = cos(pkin(6));
	t45 = t30 * t35;
	t19 = -t35 * t26 - t33 * t43;
	t17 = t19 * t30;
	t7 = -t29 * t17 + t27 * t38;
	t41 = -t27 * t17 - t29 * t38;
	t39 = t32 * r_i_i_C(1) + t34 * r_i_i_C(2) + pkin(8) + qJ(4);
	t37 = t38 * t30;
	t16 = t19 * t28;
	t15 = t38 * t28;
	t12 = -t16 * t24 + t30 * t23;
	t9 = t29 * t19 - t27 * t37;
	t6 = t27 * t19 + t29 * t37;
	t4 = t23 * t47 - t24 * t41;
	t2 = -t23 * t46 + t7 * t24;
	t1 = [0, (-t27 * t45 - t29 * t33) * pkin(2) - t39 * t41 + t36 * t9, t47, -t9, t48 * t4 + t40 * (t23 * t41 + t24 * t47), (-t4 * t32 - t9 * t34) * r_i_i_C(1) + (t9 * t32 - t4 * t34) * r_i_i_C(2); 0, (-t27 * t33 + t29 * t45) * pkin(2) + t39 * t7 + t36 * t6, -t46, -t6, t48 * t2 + t40 * (-t7 * t23 - t24 * t46), (-t2 * t32 - t6 * t34) * r_i_i_C(1) + (-t2 * t34 + t6 * t32) * r_i_i_C(2); 1, t28 * t35 * pkin(2) + t36 * t15 - t39 * t16, t30, -t15, t48 * t12 + t40 * (t16 * t23 + t30 * t24), (-t12 * t32 - t15 * t34) * r_i_i_C(1) + (-t12 * t34 + t15 * t32) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end