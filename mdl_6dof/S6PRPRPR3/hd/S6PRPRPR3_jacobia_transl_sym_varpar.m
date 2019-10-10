% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
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
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
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
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (76->27), mult. (197->52), div. (0->0), fcn. (255->10), ass. (0->23)
	t30 = pkin(8) + r_i_i_C(3);
	t15 = sin(pkin(10));
	t16 = sin(pkin(6));
	t29 = t15 * t16;
	t18 = cos(pkin(10));
	t28 = t18 * t16;
	t19 = cos(pkin(6));
	t23 = cos(qJ(2));
	t27 = t19 * t23;
	t14 = sin(pkin(11));
	t17 = cos(pkin(11));
	t21 = sin(qJ(2));
	t25 = t14 * t23 + t21 * t17;
	t10 = t25 * t19;
	t11 = t14 * t21 - t23 * t17;
	t3 = t10 * t18 - t11 * t15;
	t26 = t10 * t15 + t11 * t18;
	t20 = sin(qJ(4));
	t22 = cos(qJ(4));
	t24 = r_i_i_C(1) * t22 - r_i_i_C(2) * t20 + pkin(3);
	t9 = t11 * t19;
	t8 = t25 * t16;
	t1 = [0, -t30 * t26 + (-t15 * t27 - t18 * t21) * pkin(2) + t24 * (t15 * t9 - t18 * t25), t29, (t20 * t26 + t22 * t29) * r_i_i_C(1) + (-t20 * t29 + t22 * t26) * r_i_i_C(2), 0, 0; 0, t30 * t3 + (-t15 * t21 + t18 * t27) * pkin(2) + t24 * (-t15 * t25 - t18 * t9), -t28, (-t20 * t3 - t22 * t28) * r_i_i_C(1) + (t20 * t28 - t22 * t3) * r_i_i_C(2), 0, 0; 1, t30 * t8 + (pkin(2) * t23 - t11 * t24) * t16, t19, (t19 * t22 - t20 * t8) * r_i_i_C(1) + (-t19 * t20 - t22 * t8) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (132->29), mult. (339->52), div. (0->0), fcn. (444->10), ass. (0->28)
	t38 = pkin(4) - r_i_i_C(2);
	t37 = pkin(8) + r_i_i_C(1);
	t21 = sin(pkin(10));
	t22 = sin(pkin(6));
	t36 = t21 * t22;
	t24 = cos(pkin(10));
	t35 = t24 * t22;
	t25 = cos(pkin(6));
	t29 = cos(qJ(2));
	t34 = t25 * t29;
	t33 = r_i_i_C(3) + qJ(5);
	t20 = sin(pkin(11));
	t23 = cos(pkin(11));
	t27 = sin(qJ(2));
	t31 = t29 * t20 + t27 * t23;
	t16 = t31 * t25;
	t17 = t27 * t20 - t29 * t23;
	t7 = t24 * t16 - t21 * t17;
	t32 = t21 * t16 + t24 * t17;
	t26 = sin(qJ(4));
	t28 = cos(qJ(4));
	t30 = t33 * t26 + t38 * t28 + pkin(3);
	t15 = t17 * t25;
	t14 = t31 * t22;
	t11 = t14 * t26 - t25 * t28;
	t3 = -t26 * t32 - t28 * t36;
	t1 = t7 * t26 + t28 * t35;
	t2 = [0, -t37 * t32 + (-t21 * t34 - t24 * t27) * pkin(2) + t30 * (t21 * t15 - t24 * t31), t36, t33 * (t26 * t36 - t28 * t32) - t38 * t3, t3, 0; 0, t37 * t7 + (-t21 * t27 + t24 * t34) * pkin(2) + t30 * (-t24 * t15 - t21 * t31), -t35, t33 * (-t26 * t35 + t7 * t28) - t38 * t1, t1, 0; 1, t37 * t14 + (pkin(2) * t29 - t17 * t30) * t22, t25, t33 * (t14 * t28 + t25 * t26) - t38 * t11, t11, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (241->41), mult. (625->75), div. (0->0), fcn. (824->12), ass. (0->33)
	t24 = sin(pkin(11));
	t31 = sin(qJ(2));
	t34 = cos(qJ(2));
	t43 = cos(pkin(11));
	t37 = -t31 * t24 + t34 * t43;
	t30 = sin(qJ(4));
	t33 = cos(qJ(4));
	t29 = sin(qJ(6));
	t32 = cos(qJ(6));
	t39 = t29 * r_i_i_C(1) + t32 * r_i_i_C(2) + qJ(5);
	t42 = pkin(4) + pkin(9) + r_i_i_C(3);
	t35 = t30 * t39 + t33 * t42 + pkin(3);
	t25 = sin(pkin(10));
	t26 = sin(pkin(6));
	t47 = t25 * t26;
	t27 = cos(pkin(10));
	t46 = t27 * t26;
	t28 = cos(pkin(6));
	t45 = t28 * t34;
	t19 = -t34 * t24 - t31 * t43;
	t17 = t19 * t28;
	t7 = -t27 * t17 + t25 * t37;
	t40 = -t25 * t17 - t27 * t37;
	t38 = t32 * r_i_i_C(1) - t29 * r_i_i_C(2) + pkin(5) + pkin(8);
	t36 = t37 * t28;
	t16 = t19 * t26;
	t15 = t37 * t26;
	t11 = -t16 * t30 - t28 * t33;
	t9 = t27 * t19 - t25 * t36;
	t6 = t25 * t19 + t27 * t36;
	t3 = -t30 * t40 - t33 * t47;
	t1 = t7 * t30 + t33 * t46;
	t2 = [0, (-t25 * t45 - t27 * t31) * pkin(2) - t38 * t40 + t35 * t9, t47, t39 * (t30 * t47 - t33 * t40) - t42 * t3, t3, (t9 * t29 + t3 * t32) * r_i_i_C(1) + (-t3 * t29 + t9 * t32) * r_i_i_C(2); 0, (-t25 * t31 + t27 * t45) * pkin(2) + t38 * t7 + t35 * t6, -t46, t39 * (-t30 * t46 + t7 * t33) - t42 * t1, t1, (t1 * t32 + t6 * t29) * r_i_i_C(1) + (-t1 * t29 + t6 * t32) * r_i_i_C(2); 1, t26 * t34 * pkin(2) + t35 * t15 - t38 * t16, t28, t39 * (-t16 * t33 + t28 * t30) - t42 * t11, t11, (t11 * t32 + t15 * t29) * r_i_i_C(1) + (-t11 * t29 + t15 * t32) * r_i_i_C(2);];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end