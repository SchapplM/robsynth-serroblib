% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPPRR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (20->10), mult. (47->16), div. (0->0), fcn. (60->6), ass. (0->15)
	t17 = pkin(2) + r_i_i_C(1);
	t10 = sin(qJ(2));
	t6 = sin(pkin(10));
	t16 = t6 * t10;
	t11 = cos(qJ(2));
	t15 = t6 * t11;
	t8 = cos(pkin(10));
	t14 = t8 * t10;
	t13 = t8 * t11;
	t12 = r_i_i_C(3) + qJ(3);
	t9 = cos(pkin(6));
	t7 = sin(pkin(6));
	t3 = t15 * t9 + t14;
	t1 = -t13 * t9 + t16;
	t2 = [0, t12 * (-t16 * t9 + t13) - t17 * t3, t3, 0, 0, 0; 0, t12 * (t14 * t9 + t15) - t17 * t1, t1, 0, 0, 0; 1, (t10 * t12 + t11 * t17) * t7, -t7 * t11, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (37->15), mult. (85->24), div. (0->0), fcn. (113->8), ass. (0->15)
	t11 = cos(pkin(6));
	t12 = sin(qJ(2));
	t17 = t11 * t12;
	t13 = cos(qJ(2));
	t16 = t11 * t13;
	t6 = sin(pkin(11));
	t9 = cos(pkin(11));
	t15 = r_i_i_C(1) * t6 + r_i_i_C(2) * t9 + qJ(3);
	t14 = r_i_i_C(1) * t9 - r_i_i_C(2) * t6 + pkin(2) + pkin(3);
	t10 = cos(pkin(10));
	t8 = sin(pkin(6));
	t7 = sin(pkin(10));
	t3 = t10 * t12 + t7 * t16;
	t1 = -t10 * t16 + t12 * t7;
	t2 = [0, t15 * (t10 * t13 - t7 * t17) - t14 * t3, t3, -t7 * t8, 0, 0; 0, t15 * (t10 * t17 + t13 * t7) - t14 * t1, t1, t10 * t8, 0, 0; 1, (t15 * t12 + t14 * t13) * t8, -t8 * t13, -t11, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (90->35), mult. (222->59), div. (0->0), fcn. (291->10), ass. (0->25)
	t32 = pkin(2) + pkin(3);
	t31 = pkin(8) + r_i_i_C(3);
	t17 = sin(pkin(10));
	t18 = sin(pkin(6));
	t30 = t17 * t18;
	t20 = cos(pkin(10));
	t29 = t20 * t18;
	t21 = cos(pkin(6));
	t23 = sin(qJ(2));
	t28 = t21 * t23;
	t25 = cos(qJ(2));
	t27 = t21 * t25;
	t11 = t17 * t23 - t20 * t27;
	t12 = t17 * t25 + t20 * t28;
	t16 = sin(pkin(11));
	t19 = cos(pkin(11));
	t3 = t11 * t16 + t12 * t19;
	t13 = t17 * t27 + t20 * t23;
	t14 = -t17 * t28 + t20 * t25;
	t6 = t13 * t16 + t14 * t19;
	t22 = sin(qJ(5));
	t24 = cos(qJ(5));
	t26 = t24 * r_i_i_C(1) - t22 * r_i_i_C(2) + pkin(4);
	t10 = (-t16 * t25 + t19 * t23) * t18;
	t1 = [0, t14 * qJ(3) - t31 * t6 - t32 * t13 + t26 * (-t13 * t19 + t14 * t16), t13, -t30, (-t6 * t22 - t24 * t30) * r_i_i_C(1) + (t22 * t30 - t6 * t24) * r_i_i_C(2), 0; 0, t12 * qJ(3) - t32 * t11 - t31 * t3 + t26 * (-t11 * t19 + t12 * t16), t11, t29, (-t3 * t22 + t24 * t29) * r_i_i_C(1) + (-t22 * t29 - t3 * t24) * r_i_i_C(2), 0; 1, -t31 * t10 + (t26 * (t16 * t23 + t19 * t25) + qJ(3) * t23 + t32 * t25) * t18, -t18 * t25, -t21, (-t10 * t22 - t21 * t24) * r_i_i_C(1) + (-t10 * t24 + t21 * t22) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (212->48), mult. (542->82), div. (0->0), fcn. (715->12), ass. (0->36)
	t48 = -pkin(3) - pkin(2);
	t47 = pkin(9) + r_i_i_C(3);
	t28 = sin(pkin(10));
	t29 = sin(pkin(6));
	t46 = t28 * t29;
	t38 = cos(qJ(2));
	t45 = t29 * t38;
	t31 = cos(pkin(10));
	t44 = t31 * t29;
	t32 = cos(pkin(6));
	t35 = sin(qJ(2));
	t43 = t32 * t35;
	t42 = t32 * t38;
	t21 = t28 * t35 - t31 * t42;
	t22 = t28 * t38 + t31 * t43;
	t27 = sin(pkin(11));
	t30 = cos(pkin(11));
	t7 = -t21 * t30 + t22 * t27;
	t8 = t21 * t27 + t22 * t30;
	t23 = t28 * t42 + t31 * t35;
	t24 = -t28 * t43 + t31 * t38;
	t11 = -t23 * t30 + t24 * t27;
	t12 = t23 * t27 + t24 * t30;
	t33 = sin(qJ(6));
	t36 = cos(qJ(6));
	t41 = t36 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(5);
	t40 = t33 * r_i_i_C(1) + t36 * r_i_i_C(2) + pkin(8);
	t34 = sin(qJ(5));
	t37 = cos(qJ(5));
	t39 = t47 * t34 + t41 * t37 + pkin(4);
	t20 = t29 * t35 * t30 - t27 * t45;
	t19 = (t27 * t35 + t30 * t38) * t29;
	t14 = t20 * t37 - t32 * t34;
	t4 = t12 * t37 - t34 * t46;
	t2 = t34 * t44 + t8 * t37;
	t1 = [0, t24 * qJ(3) + t39 * t11 - t40 * t12 + t48 * t23, t23, -t46, t47 * t4 + t41 * (-t12 * t34 - t37 * t46), (t11 * t36 - t4 * t33) * r_i_i_C(1) + (-t11 * t33 - t4 * t36) * r_i_i_C(2); 0, t22 * qJ(3) + t48 * t21 + t39 * t7 - t40 * t8, t21, t44, t47 * t2 + t41 * (-t8 * t34 + t37 * t44), (-t2 * t33 + t7 * t36) * r_i_i_C(1) + (-t2 * t36 - t7 * t33) * r_i_i_C(2); 1, (t35 * qJ(3) - t48 * t38) * t29 - t40 * t20 + t39 * t19, -t45, -t32, t47 * t14 + t41 * (-t20 * t34 - t32 * t37), (-t14 * t33 + t19 * t36) * r_i_i_C(1) + (-t14 * t36 - t19 * t33) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end