% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
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
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.19s
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
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.26s
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
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (174->44), mult. (453->82), div. (0->0), fcn. (579->10), ass. (0->35)
	t29 = sin(qJ(3));
	t32 = cos(qJ(3));
	t49 = pkin(9) + r_i_i_C(2);
	t51 = pkin(3) * t32 + t49 * t29 + pkin(2);
	t50 = pkin(4) + r_i_i_C(1);
	t27 = sin(pkin(6));
	t48 = t27 * t29;
	t47 = t27 * t32;
	t33 = cos(qJ(2));
	t46 = t27 * t33;
	t28 = sin(qJ(4));
	t45 = t28 * t32;
	t31 = cos(qJ(4));
	t44 = t31 * t32;
	t43 = t32 * t33;
	t42 = r_i_i_C(3) + qJ(5);
	t41 = cos(pkin(6));
	t40 = cos(pkin(11));
	t26 = sin(pkin(11));
	t38 = t26 * t41;
	t37 = t27 * t40;
	t36 = t41 * t40;
	t34 = t42 * t28 + t50 * t31 + pkin(3);
	t30 = sin(qJ(2));
	t24 = t41 * t29 + t30 * t47;
	t22 = -t30 * t38 + t40 * t33;
	t21 = t40 * t30 + t33 * t38;
	t20 = t26 * t33 + t30 * t36;
	t19 = t26 * t30 - t33 * t36;
	t13 = t24 * t28 + t31 * t46;
	t12 = t22 * t32 + t26 * t48;
	t10 = t20 * t32 - t29 * t37;
	t3 = t12 * t28 - t21 * t31;
	t1 = t10 * t28 - t19 * t31;
	t2 = [0, t22 * pkin(8) + t50 * (-t21 * t44 + t22 * t28) + t42 * (-t21 * t45 - t22 * t31) - t51 * t21, t49 * t12 + t34 * (-t22 * t29 + t26 * t47), t42 * (t12 * t31 + t21 * t28) - t50 * t3, t3, 0; 0, t20 * pkin(8) + t50 * (-t19 * t44 + t20 * t28) + t42 * (-t19 * t45 - t20 * t31) - t51 * t19, t49 * t10 + t34 * (-t20 * t29 - t32 * t37), t42 * (t10 * t31 + t19 * t28) - t50 * t1, t1, 0; 1, (t50 * (t28 * t30 + t31 * t43) + t42 * (t28 * t43 - t30 * t31) + pkin(8) * t30 + t51 * t33) * t27, t49 * t24 + t34 * (-t30 * t48 + t41 * t32), t42 * (t24 * t31 - t28 * t46) - t50 * t13, t13, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:56:07
	% EndTime: 2019-10-09 22:56:07
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (329->56), mult. (866->104), div. (0->0), fcn. (1129->12), ass. (0->40)
	t31 = sin(qJ(3));
	t35 = cos(qJ(3));
	t45 = -r_i_i_C(3) - pkin(10) + pkin(9);
	t54 = pkin(3) * t35 + t45 * t31 + pkin(2);
	t28 = sin(pkin(6));
	t53 = t28 * t31;
	t52 = t28 * t35;
	t36 = cos(qJ(2));
	t51 = t28 * t36;
	t30 = sin(qJ(4));
	t50 = t30 * t35;
	t34 = cos(qJ(4));
	t49 = t34 * t35;
	t48 = t35 * t36;
	t47 = cos(pkin(6));
	t46 = cos(pkin(11));
	t27 = sin(pkin(11));
	t43 = t27 * t47;
	t42 = t28 * t46;
	t41 = t47 * t46;
	t29 = sin(qJ(6));
	t33 = cos(qJ(6));
	t40 = r_i_i_C(1) * t29 + r_i_i_C(2) * t33 + qJ(5);
	t39 = r_i_i_C(1) * t33 - r_i_i_C(2) * t29 + pkin(4) + pkin(5);
	t37 = t40 * t30 + t39 * t34 + pkin(3);
	t32 = sin(qJ(2));
	t24 = t47 * t31 + t32 * t52;
	t22 = -t32 * t43 + t46 * t36;
	t21 = t46 * t32 + t36 * t43;
	t20 = t27 * t36 + t32 * t41;
	t19 = t27 * t32 - t36 * t41;
	t14 = t24 * t34 - t30 * t51;
	t13 = t24 * t30 + t34 * t51;
	t12 = t22 * t35 + t27 * t53;
	t10 = t20 * t35 - t31 * t42;
	t4 = t12 * t34 + t21 * t30;
	t3 = t12 * t30 - t21 * t34;
	t2 = t10 * t34 + t19 * t30;
	t1 = t10 * t30 - t19 * t34;
	t5 = [0, t22 * pkin(8) + t40 * (-t21 * t50 - t22 * t34) + t39 * (-t21 * t49 + t22 * t30) - t54 * t21, t45 * t12 + t37 * (-t22 * t31 + t27 * t52), -t39 * t3 + t40 * t4, t3, (-t29 * t4 + t3 * t33) * r_i_i_C(1) + (-t29 * t3 - t33 * t4) * r_i_i_C(2); 0, t20 * pkin(8) + t40 * (-t19 * t50 - t20 * t34) + t39 * (-t19 * t49 + t20 * t30) - t54 * t19, t45 * t10 + t37 * (-t20 * t31 - t35 * t42), -t39 * t1 + t40 * t2, t1, (t1 * t33 - t2 * t29) * r_i_i_C(1) + (-t1 * t29 - t2 * t33) * r_i_i_C(2); 1, (t40 * (t30 * t48 - t32 * t34) + t39 * (t30 * t32 + t34 * t48) + t32 * pkin(8) + t54 * t36) * t28, t45 * t24 + t37 * (-t32 * t53 + t47 * t35), -t39 * t13 + t40 * t14, t13, (t13 * t33 - t14 * t29) * r_i_i_C(1) + (-t13 * t29 - t14 * t33) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end