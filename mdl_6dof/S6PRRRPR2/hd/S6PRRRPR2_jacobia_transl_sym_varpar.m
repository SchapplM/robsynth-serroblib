% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR2
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
% Datum: 2019-10-09 22:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (92->27), mult. (152->51), div. (0->0), fcn. (184->10), ass. (0->25)
	t14 = qJ(3) + qJ(4);
	t12 = sin(t14);
	t13 = cos(t14);
	t16 = sin(pkin(6));
	t17 = cos(pkin(11));
	t29 = t16 * t17;
	t15 = sin(pkin(11));
	t22 = cos(qJ(2));
	t18 = cos(pkin(6));
	t20 = sin(qJ(2));
	t26 = t18 * t20;
	t8 = t15 * t22 + t17 * t26;
	t34 = (-t12 * t8 - t13 * t29) * r_i_i_C(1) + (t12 * t29 - t13 * t8) * r_i_i_C(2);
	t10 = -t15 * t26 + t17 * t22;
	t30 = t15 * t16;
	t33 = (-t10 * t12 + t13 * t30) * r_i_i_C(1) + (-t10 * t13 - t12 * t30) * r_i_i_C(2);
	t28 = t16 * t20;
	t32 = (-t12 * t28 + t13 * t18) * r_i_i_C(1) + (-t12 * t18 - t13 * t28) * r_i_i_C(2);
	t31 = r_i_i_C(3) + pkin(9) + pkin(8);
	t21 = cos(qJ(3));
	t27 = t16 * t21;
	t25 = t18 * t22;
	t24 = pkin(3) * t21 + r_i_i_C(1) * t13 - r_i_i_C(2) * t12 + pkin(2);
	t19 = sin(qJ(3));
	t1 = [0, t31 * t10 + t24 * (-t15 * t25 - t17 * t20), (-t10 * t19 + t15 * t27) * pkin(3) + t33, t33, 0, 0; 0, t31 * t8 + t24 * (-t15 * t20 + t17 * t25), (-t17 * t27 - t19 * t8) * pkin(3) + t34, t34, 0, 0; 1, (t31 * t20 + t24 * t22) * t16, (t18 * t21 - t19 * t28) * pkin(3) + t32, t32, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (216->32), mult. (342->55), div. (0->0), fcn. (426->12), ass. (0->32)
	t52 = r_i_i_C(3) + qJ(5);
	t27 = sin(pkin(12));
	t30 = cos(pkin(12));
	t51 = r_i_i_C(1) * t30 - r_i_i_C(2) * t27 + pkin(4);
	t28 = sin(pkin(11));
	t29 = sin(pkin(6));
	t48 = t28 * t29;
	t31 = cos(pkin(11));
	t47 = t29 * t31;
	t34 = sin(qJ(2));
	t46 = t29 * t34;
	t35 = cos(qJ(3));
	t45 = t29 * t35;
	t32 = cos(pkin(6));
	t44 = t32 * t34;
	t36 = cos(qJ(2));
	t43 = t32 * t36;
	t42 = t27 * r_i_i_C(1) + t30 * r_i_i_C(2) + pkin(8) + pkin(9);
	t17 = t28 * t36 + t31 * t44;
	t26 = qJ(3) + qJ(4);
	t24 = sin(t26);
	t25 = cos(t26);
	t7 = t17 * t24 + t25 * t47;
	t41 = t52 * (t17 * t25 - t24 * t47) - t51 * t7;
	t19 = -t28 * t44 + t31 * t36;
	t9 = t19 * t24 - t25 * t48;
	t40 = t52 * (t19 * t25 + t24 * t48) - t51 * t9;
	t14 = t24 * t46 - t32 * t25;
	t39 = t52 * (t32 * t24 + t25 * t46) - t51 * t14;
	t38 = t35 * pkin(3) + t52 * t24 + t51 * t25 + pkin(2);
	t33 = sin(qJ(3));
	t1 = [0, t42 * t19 + t38 * (-t28 * t43 - t31 * t34), (-t19 * t33 + t28 * t45) * pkin(3) + t40, t40, t9, 0; 0, t42 * t17 + t38 * (-t28 * t34 + t31 * t43), (-t17 * t33 - t31 * t45) * pkin(3) + t41, t41, t7, 0; 1, (t42 * t34 + t38 * t36) * t29, (t32 * t35 - t33 * t46) * pkin(3) + t39, t39, t14, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (310->44), mult. (426->76), div. (0->0), fcn. (532->14), ass. (0->39)
	t63 = r_i_i_C(3) + pkin(10) + qJ(5);
	t32 = pkin(12) + qJ(6);
	t28 = sin(t32);
	t29 = cos(t32);
	t62 = cos(pkin(12)) * pkin(5) + pkin(4) + t29 * r_i_i_C(1) - t28 * r_i_i_C(2);
	t33 = qJ(3) + qJ(4);
	t30 = sin(t33);
	t31 = cos(t33);
	t41 = cos(qJ(3));
	t61 = t41 * pkin(3) + t63 * t30 + t62 * t31 + pkin(2);
	t35 = sin(pkin(11));
	t36 = sin(pkin(6));
	t57 = t35 * t36;
	t37 = cos(pkin(11));
	t56 = t36 * t37;
	t40 = sin(qJ(2));
	t55 = t36 * t40;
	t54 = t36 * t41;
	t42 = cos(qJ(2));
	t53 = t36 * t42;
	t52 = cos(pkin(6));
	t51 = t40 * t52;
	t50 = t42 * t52;
	t20 = t35 * t42 + t37 * t51;
	t10 = t20 * t31 - t30 * t56;
	t9 = t20 * t30 + t31 * t56;
	t48 = t63 * t10 - t62 * t9;
	t22 = -t35 * t51 + t37 * t42;
	t11 = t22 * t30 - t31 * t57;
	t12 = t22 * t31 + t30 * t57;
	t47 = -t62 * t11 + t63 * t12;
	t17 = t30 * t55 - t52 * t31;
	t18 = t52 * t30 + t31 * t55;
	t46 = -t62 * t17 + t63 * t18;
	t45 = sin(pkin(12)) * pkin(5) + t28 * r_i_i_C(1) + t29 * r_i_i_C(2) + pkin(9) + pkin(8);
	t39 = sin(qJ(3));
	t21 = t35 * t50 + t37 * t40;
	t19 = t35 * t40 - t37 * t50;
	t1 = [0, -t21 * t61 + t45 * t22, (-t22 * t39 + t35 * t54) * pkin(3) + t47, t47, t11, (-t12 * t28 + t21 * t29) * r_i_i_C(1) + (-t12 * t29 - t21 * t28) * r_i_i_C(2); 0, -t19 * t61 + t45 * t20, (-t20 * t39 - t37 * t54) * pkin(3) + t48, t48, t9, (-t10 * t28 + t19 * t29) * r_i_i_C(1) + (-t10 * t29 - t19 * t28) * r_i_i_C(2); 1, (t45 * t40 + t61 * t42) * t36, (-t39 * t55 + t52 * t41) * pkin(3) + t46, t46, t17, (-t18 * t28 - t29 * t53) * r_i_i_C(1) + (-t18 * t29 + t28 * t53) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end