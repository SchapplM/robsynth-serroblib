% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
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
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
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
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (70->28), mult. (123->52), div. (0->0), fcn. (150->10), ass. (0->24)
	t26 = r_i_i_C(3) + qJ(4) + pkin(8);
	t10 = sin(pkin(11));
	t11 = sin(pkin(6));
	t25 = t10 * t11;
	t12 = cos(pkin(11));
	t24 = t11 * t12;
	t16 = sin(qJ(2));
	t23 = t11 * t16;
	t17 = cos(qJ(3));
	t22 = t11 * t17;
	t13 = cos(pkin(6));
	t21 = t13 * t16;
	t18 = cos(qJ(2));
	t20 = t13 * t18;
	t9 = qJ(3) + pkin(12);
	t7 = sin(t9);
	t8 = cos(t9);
	t19 = t17 * pkin(3) + t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + pkin(2);
	t15 = sin(qJ(3));
	t4 = -t10 * t21 + t12 * t18;
	t3 = t10 * t20 + t12 * t16;
	t2 = t10 * t18 + t12 * t21;
	t1 = t10 * t16 - t12 * t20;
	t5 = [0, -t19 * t3 + t26 * t4, (t8 * t25 - t4 * t7) * r_i_i_C(1) + (-t7 * t25 - t4 * t8) * r_i_i_C(2) + (t10 * t22 - t15 * t4) * pkin(3), t3, 0, 0; 0, -t19 * t1 + t26 * t2, (-t2 * t7 - t8 * t24) * r_i_i_C(1) + (-t2 * t8 + t7 * t24) * r_i_i_C(2) + (-t12 * t22 - t15 * t2) * pkin(3), t1, 0, 0; 1, (t26 * t16 + t19 * t18) * t11, (t13 * t8 - t7 * t23) * r_i_i_C(1) + (-t13 * t7 - t8 * t23) * r_i_i_C(2) + (t13 * t17 - t15 * t23) * pkin(3), -t11 * t18, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (146->31), mult. (168->51), div. (0->0), fcn. (205->12), ass. (0->27)
	t21 = qJ(3) + pkin(12);
	t18 = qJ(5) + t21;
	t16 = sin(t18);
	t17 = cos(t18);
	t23 = sin(pkin(6));
	t24 = cos(pkin(11));
	t32 = t23 * t24;
	t22 = sin(pkin(11));
	t27 = cos(qJ(2));
	t25 = cos(pkin(6));
	t26 = sin(qJ(2));
	t30 = t25 * t26;
	t8 = t22 * t27 + t24 * t30;
	t37 = (-t16 * t8 - t17 * t32) * r_i_i_C(1) + (t16 * t32 - t17 * t8) * r_i_i_C(2);
	t10 = -t22 * t30 + t24 * t27;
	t33 = t22 * t23;
	t36 = (-t10 * t16 + t17 * t33) * r_i_i_C(1) + (-t10 * t17 - t16 * t33) * r_i_i_C(2);
	t31 = t23 * t26;
	t35 = (-t16 * t31 + t17 * t25) * r_i_i_C(1) + (-t16 * t25 - t17 * t31) * r_i_i_C(2);
	t34 = r_i_i_C(3) + pkin(9) + qJ(4) + pkin(8);
	t29 = t25 * t27;
	t13 = pkin(4) * cos(t21) + cos(qJ(3)) * pkin(3);
	t28 = r_i_i_C(1) * t17 - r_i_i_C(2) * t16 + pkin(2) + t13;
	t12 = -pkin(4) * sin(t21) - sin(qJ(3)) * pkin(3);
	t9 = t22 * t29 + t24 * t26;
	t7 = t22 * t26 - t24 * t29;
	t1 = [0, t34 * t10 - t28 * t9, t10 * t12 + t13 * t33 + t36, t9, t36, 0; 0, -t28 * t7 + t34 * t8, t12 * t8 - t13 * t32 + t37, t7, t37, 0; 1, (t34 * t26 + t28 * t27) * t23, t12 * t31 + t13 * t25 + t35, -t23 * t27, t35, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (351->46), mult. (405->71), div. (0->0), fcn. (506->14), ass. (0->37)
	t59 = pkin(10) + r_i_i_C(3);
	t35 = sin(qJ(6));
	t37 = cos(qJ(6));
	t58 = t37 * r_i_i_C(1) - t35 * r_i_i_C(2) + pkin(5);
	t31 = qJ(3) + pkin(12);
	t23 = pkin(4) * cos(t31) + cos(qJ(3)) * pkin(3);
	t28 = qJ(5) + t31;
	t26 = sin(t28);
	t27 = cos(t28);
	t57 = t59 * t26 + t58 * t27 + pkin(2) + t23;
	t32 = sin(pkin(11));
	t33 = sin(pkin(6));
	t53 = t32 * t33;
	t36 = sin(qJ(2));
	t52 = t32 * t36;
	t38 = cos(qJ(2));
	t51 = t32 * t38;
	t50 = t33 * t36;
	t49 = t33 * t38;
	t48 = cos(pkin(11));
	t47 = t33 * t48;
	t46 = t48 * t36;
	t45 = t48 * t38;
	t43 = t35 * r_i_i_C(1) + t37 * r_i_i_C(2) + pkin(8) + pkin(9) + qJ(4);
	t34 = cos(pkin(6));
	t17 = t34 * t46 + t51;
	t8 = t17 * t27 - t26 * t47;
	t42 = t59 * t8 + t58 * (-t17 * t26 - t27 * t47);
	t19 = -t34 * t52 + t45;
	t10 = t19 * t27 + t26 * t53;
	t41 = t59 * t10 + t58 * (-t19 * t26 + t27 * t53);
	t15 = t34 * t26 + t27 * t50;
	t40 = t59 * t15 + t58 * (-t26 * t50 + t34 * t27);
	t22 = -pkin(4) * sin(t31) - sin(qJ(3)) * pkin(3);
	t18 = t34 * t51 + t46;
	t16 = -t34 * t45 + t52;
	t1 = [0, -t18 * t57 + t43 * t19, t19 * t22 + t23 * t53 + t41, t18, t41, (-t10 * t35 + t18 * t37) * r_i_i_C(1) + (-t10 * t37 - t18 * t35) * r_i_i_C(2); 0, -t16 * t57 + t43 * t17, t17 * t22 - t23 * t47 + t42, t16, t42, (t16 * t37 - t8 * t35) * r_i_i_C(1) + (-t16 * t35 - t8 * t37) * r_i_i_C(2); 1, (t43 * t36 + t57 * t38) * t33, t22 * t50 + t34 * t23 + t40, -t49, t40, (-t15 * t35 - t37 * t49) * r_i_i_C(1) + (-t15 * t37 + t35 * t49) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end