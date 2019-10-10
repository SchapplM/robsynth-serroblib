% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (35->19), mult. (92->37), div. (0->0), fcn. (112->8), ass. (0->19)
	t6 = sin(pkin(6));
	t9 = sin(qJ(3));
	t20 = t6 * t9;
	t19 = pkin(8) + r_i_i_C(3);
	t11 = cos(qJ(3));
	t18 = t11 * t6;
	t10 = sin(qJ(2));
	t5 = sin(pkin(11));
	t17 = t5 * t10;
	t12 = cos(qJ(2));
	t16 = t5 * t12;
	t7 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (66->21), mult. (169->39), div. (0->0), fcn. (211->8), ass. (0->22)
	t28 = pkin(3) + r_i_i_C(1);
	t27 = pkin(8) + r_i_i_C(2);
	t14 = sin(pkin(6));
	t17 = sin(qJ(3));
	t26 = t14 * t17;
	t19 = cos(qJ(3));
	t25 = t14 * t19;
	t16 = cos(pkin(6));
	t18 = sin(qJ(2));
	t24 = t16 * t18;
	t20 = cos(qJ(2));
	t23 = t16 * t20;
	t22 = r_i_i_C(3) + qJ(4);
	t21 = t22 * t17 + t19 * t28 + pkin(2);
	t15 = cos(pkin(11));
	t13 = sin(pkin(11));
	t9 = -t16 * t19 + t18 * t26;
	t8 = -t13 * t24 + t15 * t20;
	t6 = t13 * t20 + t15 * t24;
	t3 = -t13 * t25 + t17 * t8;
	t1 = t15 * t25 + t17 * t6;
	t2 = [0, t27 * t8 + t21 * (-t13 * t23 - t15 * t18), t22 * (t13 * t26 + t19 * t8) - t28 * t3, t3, 0, 0; 0, t27 * t6 + t21 * (-t13 * t18 + t15 * t23), t22 * (-t15 * t26 + t19 * t6) - t28 * t1, t1, 0, 0; 1, (t18 * t27 + t20 * t21) * t14, -t28 * t9 + t22 * (t16 * t17 + t18 * t25), t9, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (139->33), mult. (362->61), div. (0->0), fcn. (465->10), ass. (0->27)
	t15 = sin(pkin(6));
	t19 = sin(qJ(3));
	t31 = t15 * t19;
	t22 = cos(qJ(3));
	t30 = t15 * t22;
	t17 = cos(pkin(6));
	t20 = sin(qJ(2));
	t29 = t17 * t20;
	t23 = cos(qJ(2));
	t28 = t17 * t23;
	t27 = -r_i_i_C(3) - pkin(9) + pkin(8);
	t18 = sin(qJ(5));
	t21 = cos(qJ(5));
	t26 = t18 * r_i_i_C(1) + t21 * r_i_i_C(2) + qJ(4);
	t25 = t21 * r_i_i_C(1) - t18 * r_i_i_C(2) + pkin(3) + pkin(4);
	t24 = t26 * t19 + t25 * t22 + pkin(2);
	t16 = cos(pkin(11));
	t14 = sin(pkin(11));
	t10 = t17 * t19 + t20 * t30;
	t9 = -t17 * t22 + t20 * t31;
	t8 = -t14 * t29 + t16 * t23;
	t6 = t14 * t23 + t16 * t29;
	t4 = t14 * t31 + t8 * t22;
	t3 = -t14 * t30 + t8 * t19;
	t2 = -t16 * t31 + t6 * t22;
	t1 = t16 * t30 + t6 * t19;
	t5 = [0, t27 * t8 + t24 * (-t14 * t28 - t16 * t20), -t25 * t3 + t26 * t4, t3, (-t4 * t18 + t3 * t21) * r_i_i_C(1) + (-t3 * t18 - t4 * t21) * r_i_i_C(2), 0; 0, t27 * t6 + t24 * (-t14 * t20 + t16 * t28), -t25 * t1 + t26 * t2, t1, (t1 * t21 - t2 * t18) * r_i_i_C(1) + (-t1 * t18 - t2 * t21) * r_i_i_C(2), 0; 1, (t27 * t20 + t24 * t23) * t15, t26 * t10 - t25 * t9, t9, (-t10 * t18 + t9 * t21) * r_i_i_C(1) + (-t10 * t21 - t9 * t18) * r_i_i_C(2), 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (309->51), mult. (815->88), div. (0->0), fcn. (1061->12), ass. (0->42)
	t35 = sin(qJ(5));
	t36 = sin(qJ(3));
	t39 = cos(qJ(5));
	t40 = cos(qJ(3));
	t34 = sin(qJ(6));
	t38 = cos(qJ(6));
	t46 = t38 * r_i_i_C(1) - t34 * r_i_i_C(2) + pkin(5);
	t61 = r_i_i_C(3) + pkin(10);
	t67 = t61 * (t35 * t40 - t36 * t39) + (t35 * t36 + t39 * t40) * t46;
	t62 = pkin(3) + pkin(4);
	t44 = t36 * qJ(4) + t62 * t40 + pkin(2);
	t66 = t44 + t67;
	t37 = sin(qJ(2));
	t57 = cos(pkin(6));
	t32 = sin(pkin(6));
	t60 = t32 * t36;
	t27 = t37 * t60 - t57 * t40;
	t59 = t32 * t40;
	t28 = t57 * t36 + t37 * t59;
	t16 = t27 * t35 + t28 * t39;
	t65 = t46 * (t27 * t39 - t28 * t35) + t61 * t16;
	t33 = cos(pkin(11));
	t41 = cos(qJ(2));
	t56 = sin(pkin(11));
	t47 = t57 * t56;
	t26 = t33 * t41 - t37 * t47;
	t52 = t32 * t56;
	t19 = t26 * t36 - t40 * t52;
	t20 = t26 * t40 + t36 * t52;
	t8 = t19 * t35 + t20 * t39;
	t64 = t46 * (t19 * t39 - t20 * t35) + t61 * t8;
	t51 = t33 * t57;
	t24 = t37 * t51 + t56 * t41;
	t17 = t24 * t36 + t33 * t59;
	t18 = t24 * t40 - t33 * t60;
	t4 = t17 * t35 + t18 * t39;
	t63 = t46 * (t17 * t39 - t18 * t35) + t61 * t4;
	t58 = t32 * t41;
	t45 = -t34 * r_i_i_C(1) - t38 * r_i_i_C(2) + pkin(8) - pkin(9);
	t25 = -t33 * t37 - t41 * t47;
	t23 = -t56 * t37 + t41 * t51;
	t1 = [0, t66 * t25 + t45 * t26, t20 * qJ(4) - t62 * t19 - t64, t19, t64, (t25 * t38 - t8 * t34) * r_i_i_C(1) + (-t25 * t34 - t8 * t38) * r_i_i_C(2); 0, t66 * t23 + t45 * t24, t18 * qJ(4) - t62 * t17 - t63, t17, t63, (t23 * t38 - t4 * t34) * r_i_i_C(1) + (-t23 * t34 - t4 * t38) * r_i_i_C(2); 1, (t45 * t37 + t44 * t41) * t32 + t67 * t58, t28 * qJ(4) - t62 * t27 - t65, t27, t65, (-t16 * t34 + t38 * t58) * r_i_i_C(1) + (-t16 * t38 - t34 * t58) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end