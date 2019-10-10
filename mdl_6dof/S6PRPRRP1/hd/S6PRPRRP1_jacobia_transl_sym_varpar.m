% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRP1
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
% Datum: 2019-10-09 21:42
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRP1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:42:44
	% EndTime: 2019-10-09 21:42:44
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:42:44
	% EndTime: 2019-10-09 21:42:44
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:42:44
	% EndTime: 2019-10-09 21:42:44
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
	% StartTime: 2019-10-09 21:42:44
	% EndTime: 2019-10-09 21:42:44
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-09 21:42:44
	% EndTime: 2019-10-09 21:42:45
	% DurationCPUTime: 0.17s
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
	% StartTime: 2019-10-09 21:42:44
	% EndTime: 2019-10-09 21:42:45
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (198->40), mult. (517->75), div. (0->0), fcn. (679->12), ass. (0->33)
	t22 = sin(pkin(11));
	t29 = sin(qJ(2));
	t32 = cos(qJ(2));
	t40 = cos(pkin(11));
	t35 = -t29 * t22 + t32 * t40;
	t28 = sin(qJ(4));
	t31 = cos(qJ(4));
	t27 = sin(qJ(5));
	t30 = cos(qJ(5));
	t37 = t30 * r_i_i_C(1) - t27 * r_i_i_C(2) + pkin(4);
	t45 = pkin(9) + r_i_i_C(3);
	t33 = t28 * t45 + t37 * t31 + pkin(3);
	t23 = sin(pkin(10));
	t24 = sin(pkin(6));
	t44 = t23 * t24;
	t25 = cos(pkin(10));
	t43 = t25 * t24;
	t26 = cos(pkin(6));
	t42 = t26 * t32;
	t19 = -t32 * t22 - t29 * t40;
	t17 = t19 * t26;
	t7 = -t25 * t17 + t23 * t35;
	t38 = -t23 * t17 - t25 * t35;
	t36 = t27 * r_i_i_C(1) + t30 * r_i_i_C(2) + pkin(8);
	t34 = t35 * t26;
	t16 = t19 * t24;
	t15 = t35 * t24;
	t12 = -t16 * t31 + t26 * t28;
	t9 = t25 * t19 - t23 * t34;
	t6 = t23 * t19 + t25 * t34;
	t4 = t28 * t44 - t31 * t38;
	t2 = -t28 * t43 + t7 * t31;
	t1 = [0, (-t23 * t42 - t25 * t29) * pkin(2) - t36 * t38 + t33 * t9, t44, t45 * t4 + t37 * (t28 * t38 + t31 * t44), (-t4 * t27 - t9 * t30) * r_i_i_C(1) + (t9 * t27 - t4 * t30) * r_i_i_C(2), 0; 0, (-t23 * t29 + t25 * t42) * pkin(2) + t36 * t7 + t33 * t6, -t43, t45 * t2 + t37 * (-t7 * t28 - t31 * t43), (-t2 * t27 - t6 * t30) * r_i_i_C(1) + (-t2 * t30 + t6 * t27) * r_i_i_C(2), 0; 1, t24 * t32 * pkin(2) + t33 * t15 - t36 * t16, t26, t45 * t12 + t37 * (t16 * t28 + t26 * t31), (-t12 * t27 - t15 * t30) * r_i_i_C(1) + (-t12 * t30 + t15 * t27) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:42:44
	% EndTime: 2019-10-09 21:42:45
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (253->41), mult. (637->75), div. (0->0), fcn. (836->12), ass. (0->37)
	t54 = pkin(5) + r_i_i_C(1);
	t25 = sin(pkin(11));
	t33 = sin(qJ(2));
	t36 = cos(qJ(2));
	t48 = cos(pkin(11));
	t40 = -t33 * t25 + t36 * t48;
	t31 = sin(qJ(5));
	t34 = cos(qJ(5));
	t55 = t34 * r_i_i_C(2) + t54 * t31 + pkin(8);
	t32 = sin(qJ(4));
	t35 = cos(qJ(4));
	t41 = -r_i_i_C(2) * t31 + t54 * t34 + pkin(4);
	t53 = r_i_i_C(3) + qJ(6) + pkin(9);
	t37 = t53 * t32 + t41 * t35 + pkin(3);
	t26 = sin(pkin(10));
	t27 = sin(pkin(6));
	t52 = t26 * t27;
	t28 = cos(pkin(10));
	t51 = t28 * t27;
	t29 = cos(pkin(6));
	t50 = t29 * t36;
	t19 = -t36 * t25 - t33 * t48;
	t17 = t19 * t29;
	t7 = -t17 * t28 + t26 * t40;
	t42 = -t26 * t17 - t28 * t40;
	t38 = t29 * t40;
	t16 = t19 * t27;
	t15 = t40 * t27;
	t12 = -t16 * t35 + t29 * t32;
	t11 = -t16 * t32 - t29 * t35;
	t9 = t19 * t28 - t26 * t38;
	t6 = t26 * t19 + t28 * t38;
	t4 = t32 * t52 - t35 * t42;
	t3 = -t32 * t42 - t35 * t52;
	t2 = -t32 * t51 + t35 * t7;
	t1 = t32 * t7 + t35 * t51;
	t5 = [0, (-t26 * t50 - t28 * t33) * pkin(2) - t55 * t42 + t37 * t9, t52, -t41 * t3 + t53 * t4, (t31 * t9 - t34 * t4) * r_i_i_C(2) + t54 * (-t31 * t4 - t34 * t9), t3; 0, (-t26 * t33 + t28 * t50) * pkin(2) + t55 * t7 + t37 * t6, -t51, -t41 * t1 + t53 * t2, (-t2 * t34 + t31 * t6) * r_i_i_C(2) + t54 * (-t2 * t31 - t34 * t6), t1; 1, t27 * t36 * pkin(2) + t37 * t15 - t55 * t16, t29, -t41 * t11 + t53 * t12, (-t12 * t34 + t15 * t31) * r_i_i_C(2) + t54 * (-t12 * t31 - t15 * t34), t11;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end