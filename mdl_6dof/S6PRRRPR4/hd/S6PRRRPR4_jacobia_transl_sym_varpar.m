% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.15s
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
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (102->32), mult. (269->63), div. (0->0), fcn. (338->10), ass. (0->28)
	t15 = sin(qJ(3));
	t18 = cos(qJ(3));
	t14 = sin(qJ(4));
	t17 = cos(qJ(4));
	t22 = t17 * r_i_i_C(1) - t14 * r_i_i_C(2) + pkin(3);
	t31 = pkin(9) + r_i_i_C(3);
	t32 = t31 * t15 + t22 * t18 + pkin(2);
	t13 = sin(pkin(6));
	t30 = t13 * t15;
	t29 = t13 * t18;
	t19 = cos(qJ(2));
	t28 = t13 * t19;
	t27 = cos(pkin(6));
	t26 = cos(pkin(11));
	t12 = sin(pkin(11));
	t25 = t12 * t27;
	t24 = t13 * t26;
	t23 = t27 * t26;
	t21 = t14 * r_i_i_C(1) + t17 * r_i_i_C(2) + pkin(8);
	t16 = sin(qJ(2));
	t10 = t27 * t15 + t16 * t29;
	t8 = -t16 * t25 + t26 * t19;
	t7 = t26 * t16 + t19 * t25;
	t6 = t12 * t19 + t16 * t23;
	t5 = t12 * t16 - t19 * t23;
	t4 = t12 * t30 + t8 * t18;
	t2 = -t15 * t24 + t6 * t18;
	t1 = [0, t21 * t8 - t32 * t7, t31 * t4 + t22 * (t12 * t29 - t8 * t15), (-t4 * t14 + t7 * t17) * r_i_i_C(1) + (-t7 * t14 - t4 * t17) * r_i_i_C(2), 0, 0; 0, t21 * t6 - t32 * t5, t31 * t2 + t22 * (-t6 * t15 - t18 * t24), (-t2 * t14 + t5 * t17) * r_i_i_C(1) + (-t5 * t14 - t2 * t17) * r_i_i_C(2), 0, 0; 1, (t21 * t16 + t32 * t19) * t13, t31 * t10 + t22 * (-t16 * t30 + t27 * t18), (-t10 * t14 - t17 * t28) * r_i_i_C(1) + (-t10 * t17 + t14 * t28) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (167->41), mult. (337->74), div. (0->0), fcn. (423->12), ass. (0->34)
	t22 = sin(qJ(3));
	t25 = cos(qJ(3));
	t17 = qJ(4) + pkin(12);
	t15 = sin(t17);
	t16 = cos(t17);
	t24 = cos(qJ(4));
	t29 = t24 * pkin(4) + t16 * r_i_i_C(1) - t15 * r_i_i_C(2) + pkin(3);
	t38 = r_i_i_C(3) + qJ(5) + pkin(9);
	t39 = t38 * t22 + t29 * t25 + pkin(2);
	t19 = sin(pkin(6));
	t37 = t19 * t22;
	t36 = t19 * t25;
	t26 = cos(qJ(2));
	t35 = t19 * t26;
	t34 = cos(pkin(6));
	t33 = cos(pkin(11));
	t18 = sin(pkin(11));
	t32 = t18 * t34;
	t31 = t19 * t33;
	t30 = t34 * t33;
	t21 = sin(qJ(4));
	t28 = t21 * pkin(4) + t15 * r_i_i_C(1) + t16 * r_i_i_C(2) + pkin(8);
	t23 = sin(qJ(2));
	t10 = t34 * t22 + t23 * t36;
	t9 = t23 * t37 - t34 * t25;
	t8 = -t23 * t32 + t33 * t26;
	t7 = t33 * t23 + t26 * t32;
	t6 = t18 * t26 + t23 * t30;
	t5 = t18 * t23 - t26 * t30;
	t4 = t18 * t37 + t8 * t25;
	t3 = -t18 * t36 + t8 * t22;
	t2 = -t22 * t31 + t6 * t25;
	t1 = t6 * t22 + t25 * t31;
	t11 = [0, t28 * t8 - t39 * t7, -t29 * t3 + t38 * t4, (-t4 * t15 + t7 * t16) * r_i_i_C(1) + (-t7 * t15 - t4 * t16) * r_i_i_C(2) + (-t4 * t21 + t7 * t24) * pkin(4), t3, 0; 0, t28 * t6 - t39 * t5, -t29 * t1 + t38 * t2, (-t2 * t15 + t5 * t16) * r_i_i_C(1) + (-t5 * t15 - t2 * t16) * r_i_i_C(2) + (-t2 * t21 + t5 * t24) * pkin(4), t1, 0; 1, (t28 * t23 + t39 * t26) * t19, t38 * t10 - t29 * t9, (-t10 * t15 - t16 * t35) * r_i_i_C(1) + (-t10 * t16 + t15 * t35) * r_i_i_C(2) + (-t10 * t21 - t24 * t35) * pkin(4), t9, 0;];
	Ja_transl = t11;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:52:15
	% EndTime: 2019-10-09 22:52:15
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (280->44), mult. (414->73), div. (0->0), fcn. (520->14), ass. (0->38)
	t32 = sin(qJ(3));
	t34 = cos(qJ(3));
	t29 = qJ(4) + pkin(12);
	t19 = pkin(5) * cos(t29) + cos(qJ(4)) * pkin(4);
	t26 = qJ(6) + t29;
	t24 = sin(t26);
	t25 = cos(t26);
	t38 = r_i_i_C(1) * t25 - r_i_i_C(2) * t24 + pkin(3) + t19;
	t47 = r_i_i_C(3) + pkin(10) + qJ(5) + pkin(9);
	t51 = t47 * t32 + t38 * t34 + pkin(2);
	t30 = sin(pkin(11));
	t33 = sin(qJ(2));
	t35 = cos(qJ(2));
	t42 = cos(pkin(11));
	t43 = cos(pkin(6));
	t39 = t43 * t42;
	t11 = t30 * t33 - t35 * t39;
	t12 = t30 * t35 + t33 * t39;
	t31 = sin(pkin(6));
	t40 = t31 * t42;
	t8 = t12 * t34 - t32 * t40;
	t50 = (t11 * t25 - t8 * t24) * r_i_i_C(1) + (-t11 * t24 - t8 * t25) * r_i_i_C(2);
	t41 = t30 * t43;
	t14 = -t33 * t41 + t42 * t35;
	t46 = t31 * t32;
	t10 = t14 * t34 + t30 * t46;
	t13 = t42 * t33 + t35 * t41;
	t49 = (-t10 * t24 + t13 * t25) * r_i_i_C(1) + (-t10 * t25 - t13 * t24) * r_i_i_C(2);
	t45 = t31 * t34;
	t16 = t43 * t32 + t33 * t45;
	t44 = t31 * t35;
	t48 = (-t16 * t24 - t25 * t44) * r_i_i_C(1) + (-t16 * t25 + t24 * t44) * r_i_i_C(2);
	t18 = pkin(5) * sin(t29) + sin(qJ(4)) * pkin(4);
	t37 = t24 * r_i_i_C(1) + t25 * r_i_i_C(2) + pkin(8) + t18;
	t15 = t33 * t46 - t43 * t34;
	t9 = t14 * t32 - t30 * t45;
	t7 = t12 * t32 + t34 * t40;
	t1 = [0, -t13 * t51 + t37 * t14, t47 * t10 - t38 * t9, -t10 * t18 + t13 * t19 + t49, t9, t49; 0, -t11 * t51 + t37 * t12, -t38 * t7 + t47 * t8, t11 * t19 - t8 * t18 + t50, t7, t50; 1, (t37 * t33 + t51 * t35) * t31, -t38 * t15 + t47 * t16, -t16 * t18 - t19 * t44 + t48, t15, t48;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end