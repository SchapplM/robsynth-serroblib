% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:58
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:58
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:58
	% EndTime: 2019-10-09 21:29:59
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
	% StartTime: 2019-10-09 21:29:58
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-09 21:29:58
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.22s
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
	% StartTime: 2019-10-09 21:29:58
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (120->37), mult. (254->63), div. (0->0), fcn. (329->12), ass. (0->29)
	t19 = sin(pkin(11));
	t22 = cos(pkin(11));
	t27 = sin(qJ(2));
	t29 = cos(qJ(2));
	t11 = t27 * t19 - t29 * t22;
	t38 = r_i_i_C(3) + qJ(5) + pkin(8);
	t20 = sin(pkin(10));
	t21 = sin(pkin(6));
	t37 = t20 * t21;
	t23 = cos(pkin(10));
	t36 = t23 * t21;
	t24 = cos(pkin(6));
	t35 = t24 * t29;
	t31 = t29 * t19 + t27 * t22;
	t10 = t31 * t24;
	t3 = t23 * t10 - t20 * t11;
	t32 = t20 * t10 + t23 * t11;
	t18 = qJ(4) + pkin(12);
	t16 = sin(t18);
	t17 = cos(t18);
	t28 = cos(qJ(4));
	t30 = t28 * pkin(4) + t17 * r_i_i_C(1) - t16 * r_i_i_C(2) + pkin(3);
	t26 = sin(qJ(4));
	t9 = t11 * t24;
	t8 = t31 * t21;
	t7 = t11 * t21;
	t5 = t20 * t9 - t23 * t31;
	t2 = -t20 * t31 - t23 * t9;
	t1 = [0, -t38 * t32 + (-t20 * t35 - t23 * t27) * pkin(2) + t30 * t5, t37, (t16 * t32 + t17 * t37) * r_i_i_C(1) + (-t16 * t37 + t17 * t32) * r_i_i_C(2) + (t26 * t32 + t28 * t37) * pkin(4), -t5, 0; 0, t38 * t3 + (-t20 * t27 + t23 * t35) * pkin(2) + t30 * t2, -t36, (-t3 * t16 - t17 * t36) * r_i_i_C(1) + (t16 * t36 - t3 * t17) * r_i_i_C(2) + (-t3 * t26 - t28 * t36) * pkin(4), -t2, 0; 1, t21 * t29 * pkin(2) - t30 * t7 + t38 * t8, t24, (-t8 * t16 + t24 * t17) * r_i_i_C(1) + (-t24 * t16 - t8 * t17) * r_i_i_C(2) + (t24 * t28 - t8 * t26) * pkin(4), t7, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:58
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (281->52), mult. (574->85), div. (0->0), fcn. (753->14), ass. (0->36)
	t26 = sin(pkin(11));
	t34 = sin(qJ(2));
	t37 = cos(qJ(2));
	t45 = cos(pkin(11));
	t40 = -t34 * t26 + t37 * t45;
	t25 = qJ(4) + pkin(12);
	t23 = sin(t25);
	t24 = cos(t25);
	t36 = cos(qJ(4));
	t32 = sin(qJ(6));
	t35 = cos(qJ(6));
	t42 = t35 * r_i_i_C(1) - t32 * r_i_i_C(2) + pkin(5);
	t50 = pkin(9) + r_i_i_C(3);
	t38 = t36 * pkin(4) + t50 * t23 + t42 * t24 + pkin(3);
	t27 = sin(pkin(10));
	t28 = sin(pkin(6));
	t49 = t27 * t28;
	t29 = cos(pkin(10));
	t48 = t29 * t28;
	t30 = cos(pkin(6));
	t47 = t30 * t37;
	t19 = -t37 * t26 - t34 * t45;
	t17 = t19 * t30;
	t7 = -t29 * t17 + t27 * t40;
	t43 = -t27 * t17 - t29 * t40;
	t41 = t32 * r_i_i_C(1) + t35 * r_i_i_C(2) + pkin(8) + qJ(5);
	t39 = t40 * t30;
	t33 = sin(qJ(4));
	t16 = t19 * t28;
	t15 = t40 * t28;
	t12 = -t16 * t24 + t30 * t23;
	t9 = t29 * t19 - t27 * t39;
	t6 = t27 * t19 + t29 * t39;
	t4 = t23 * t49 - t24 * t43;
	t2 = -t23 * t48 + t7 * t24;
	t1 = [0, (-t27 * t47 - t29 * t34) * pkin(2) - t41 * t43 + t38 * t9, t49, t50 * t4 + (t33 * t43 + t36 * t49) * pkin(4) + t42 * (t23 * t43 + t24 * t49), -t9, (-t4 * t32 - t9 * t35) * r_i_i_C(1) + (t9 * t32 - t4 * t35) * r_i_i_C(2); 0, (-t27 * t34 + t29 * t47) * pkin(2) + t41 * t7 + t38 * t6, -t48, t50 * t2 + (-t33 * t7 - t36 * t48) * pkin(4) + t42 * (-t7 * t23 - t24 * t48), -t6, (-t2 * t32 - t6 * t35) * r_i_i_C(1) + (-t2 * t35 + t6 * t32) * r_i_i_C(2); 1, t28 * t37 * pkin(2) + t15 * t38 - t41 * t16, t30, t50 * t12 + (t16 * t33 + t30 * t36) * pkin(4) + t42 * (t16 * t23 + t30 * t24), -t15, (-t12 * t32 - t15 * t35) * r_i_i_C(1) + (-t12 * t35 + t15 * t32) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end