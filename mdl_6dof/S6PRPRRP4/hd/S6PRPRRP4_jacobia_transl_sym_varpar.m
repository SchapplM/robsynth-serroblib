% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:48:18
	% EndTime: 2019-10-09 21:48:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:48:18
	% EndTime: 2019-10-09 21:48:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:48:18
	% EndTime: 2019-10-09 21:48:18
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
	% StartTime: 2019-10-09 21:48:18
	% EndTime: 2019-10-09 21:48:18
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (25->11), mult. (63->20), div. (0->0), fcn. (78->8), ass. (0->13)
	t11 = cos(pkin(6));
	t12 = sin(qJ(2));
	t17 = t11 * t12;
	t13 = cos(qJ(2));
	t16 = t11 * t13;
	t15 = r_i_i_C(3) + qJ(3);
	t14 = r_i_i_C(1) * cos(pkin(11)) - r_i_i_C(2) * sin(pkin(11)) + pkin(2);
	t10 = cos(pkin(10));
	t8 = sin(pkin(6));
	t7 = sin(pkin(10));
	t3 = t10 * t12 + t7 * t16;
	t1 = -t10 * t16 + t7 * t12;
	t2 = [0, t15 * (t10 * t13 - t7 * t17) - t14 * t3, t3, 0, 0, 0; 0, t15 * (t10 * t17 + t7 * t13) - t14 * t1, t1, 0, 0, 0; 1, (t15 * t12 + t14 * t13) * t8, -t8 * t13, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:48:18
	% EndTime: 2019-10-09 21:48:18
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (62->22), mult. (102->42), div. (0->0), fcn. (127->9), ass. (0->21)
	t23 = r_i_i_C(3) + pkin(8) + qJ(3);
	t10 = sin(pkin(10));
	t11 = sin(pkin(6));
	t22 = t10 * t11;
	t12 = cos(pkin(10));
	t21 = t11 * t12;
	t15 = sin(qJ(2));
	t20 = t11 * t15;
	t13 = cos(pkin(6));
	t19 = t13 * t15;
	t16 = cos(qJ(2));
	t18 = t13 * t16;
	t9 = pkin(11) + qJ(4);
	t7 = sin(t9);
	t8 = cos(t9);
	t17 = r_i_i_C(1) * t8 - r_i_i_C(2) * t7 + cos(pkin(11)) * pkin(3) + pkin(2);
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
	% StartTime: 2019-10-09 21:48:18
	% EndTime: 2019-10-09 21:48:18
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (168->36), mult. (279->64), div. (0->0), fcn. (353->11), ass. (0->29)
	t15 = pkin(11) + qJ(4);
	t13 = sin(t15);
	t14 = cos(t15);
	t19 = sin(qJ(5));
	t21 = cos(qJ(5));
	t25 = t21 * r_i_i_C(1) - t19 * r_i_i_C(2) + pkin(4);
	t34 = pkin(9) + r_i_i_C(3);
	t35 = t34 * t13 + t25 * t14 + cos(pkin(11)) * pkin(3) + pkin(2);
	t16 = sin(pkin(10));
	t17 = sin(pkin(6));
	t33 = t16 * t17;
	t20 = sin(qJ(2));
	t32 = t17 * t20;
	t22 = cos(qJ(2));
	t31 = t17 * t22;
	t30 = cos(pkin(6));
	t29 = cos(pkin(10));
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
	% StartTime: 2019-10-09 21:48:18
	% EndTime: 2019-10-09 21:48:18
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (276->48), mult. (463->84), div. (0->0), fcn. (594->11), ass. (0->37)
	t29 = pkin(11) + qJ(4);
	t27 = sin(t29);
	t28 = cos(t29);
	t52 = pkin(9) + r_i_i_C(2);
	t54 = t52 * t27 + pkin(4) * t28 + cos(pkin(11)) * pkin(3) + pkin(2);
	t53 = pkin(5) + r_i_i_C(1);
	t33 = sin(qJ(5));
	t51 = t28 * t33;
	t35 = cos(qJ(5));
	t50 = t28 * t35;
	t30 = sin(pkin(10));
	t31 = sin(pkin(6));
	t49 = t30 * t31;
	t34 = sin(qJ(2));
	t48 = t31 * t34;
	t36 = cos(qJ(2));
	t47 = t33 * t36;
	t46 = t35 * t36;
	t45 = r_i_i_C(3) + qJ(6);
	t44 = cos(pkin(6));
	t43 = cos(pkin(10));
	t41 = t30 * t44;
	t40 = t31 * t43;
	t39 = t44 * t43;
	t37 = t45 * t33 + t53 * t35 + pkin(4);
	t32 = -pkin(8) - qJ(3);
	t24 = -t34 * t41 + t43 * t36;
	t23 = t43 * t34 + t36 * t41;
	t22 = t30 * t36 + t34 * t39;
	t21 = t30 * t34 - t36 * t39;
	t18 = t44 * t27 + t28 * t48;
	t13 = t18 * t33 + t31 * t46;
	t12 = t24 * t28 + t27 * t49;
	t10 = t22 * t28 - t27 * t40;
	t3 = t12 * t33 - t23 * t35;
	t1 = t10 * t33 - t21 * t35;
	t2 = [0, -t24 * t32 + t53 * (-t23 * t50 + t24 * t33) + t45 * (-t23 * t51 - t24 * t35) - t54 * t23, t23, t52 * t12 + t37 * (-t24 * t27 + t28 * t49), t45 * (t12 * t35 + t23 * t33) - t53 * t3, t3; 0, -t22 * t32 + t53 * (-t21 * t50 + t22 * t33) + t45 * (-t21 * t51 - t22 * t35) - t54 * t21, t21, t52 * t10 + t37 * (-t22 * t27 - t28 * t40), t45 * (t10 * t35 + t21 * t33) - t53 * t1, t1; 1, (t53 * (t28 * t46 + t33 * t34) + t45 * (t28 * t47 - t34 * t35) - t32 * t34 + t54 * t36) * t31, -t31 * t36, t52 * t18 + t37 * (-t27 * t48 + t44 * t28), t45 * (t18 * t35 - t31 * t47) - t53 * t13, t13;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end