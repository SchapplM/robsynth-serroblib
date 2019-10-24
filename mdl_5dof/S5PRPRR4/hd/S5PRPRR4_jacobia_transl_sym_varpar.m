% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRR4
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:26
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRPRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:44
	% EndTime: 2019-10-24 10:26:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:44
	% EndTime: 2019-10-24 10:26:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:44
	% EndTime: 2019-10-24 10:26:44
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(5));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(9));
	t1 = sin(pkin(9));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(5)), 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:44
	% EndTime: 2019-10-24 10:26:44
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (23->15), mult. (60->31), div. (0->0), fcn. (77->8), ass. (0->14)
	t10 = cos(pkin(5));
	t12 = cos(qJ(2));
	t13 = t10 * t12;
	t11 = sin(qJ(2));
	t5 = sin(pkin(10));
	t8 = cos(pkin(10));
	t4 = -t11 * t8 - t12 * t5;
	t3 = t11 * t5 - t12 * t8;
	t9 = cos(pkin(9));
	t7 = sin(pkin(5));
	t6 = sin(pkin(9));
	t2 = t4 * t10;
	t1 = t3 * t10;
	t14 = [0, (t6 * t1 + t9 * t4) * r_i_i_C(1) + (-t6 * t2 + t9 * t3) * r_i_i_C(2) + (-t9 * t11 - t6 * t13) * pkin(2), t6 * t7, 0, 0; 0, (-t9 * t1 + t6 * t4) * r_i_i_C(1) + (t9 * t2 + t6 * t3) * r_i_i_C(2) + (-t6 * t11 + t9 * t13) * pkin(2), -t9 * t7, 0, 0; 1, (t12 * pkin(2) - t3 * r_i_i_C(1) + t4 * r_i_i_C(2)) * t7, t10, 0, 0;];
	Ja_transl = t14;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:44
	% EndTime: 2019-10-24 10:26:45
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (76->27), mult. (197->52), div. (0->0), fcn. (255->10), ass. (0->23)
	t30 = pkin(7) + r_i_i_C(3);
	t15 = sin(pkin(9));
	t16 = sin(pkin(5));
	t29 = t15 * t16;
	t18 = cos(pkin(9));
	t28 = t18 * t16;
	t19 = cos(pkin(5));
	t23 = cos(qJ(2));
	t27 = t19 * t23;
	t14 = sin(pkin(10));
	t17 = cos(pkin(10));
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
	t1 = [0, -t30 * t26 + (-t15 * t27 - t18 * t21) * pkin(2) + t24 * (t15 * t9 - t18 * t25), t29, (t20 * t26 + t22 * t29) * r_i_i_C(1) + (-t20 * t29 + t22 * t26) * r_i_i_C(2), 0; 0, t30 * t3 + (-t15 * t21 + t18 * t27) * pkin(2) + t24 * (-t15 * t25 - t18 * t9), -t28, (-t20 * t3 - t22 * t28) * r_i_i_C(1) + (t20 * t28 - t22 * t3) * r_i_i_C(2), 0; 1, t30 * t8 + (pkin(2) * t23 - t11 * t24) * t16, t19, (t19 * t22 - t20 * t8) * r_i_i_C(1) + (-t19 * t20 - t22 * t8) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:26:44
	% EndTime: 2019-10-24 10:26:45
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (198->40), mult. (517->75), div. (0->0), fcn. (679->12), ass. (0->33)
	t22 = sin(pkin(10));
	t29 = sin(qJ(2));
	t32 = cos(qJ(2));
	t40 = cos(pkin(10));
	t35 = -t29 * t22 + t32 * t40;
	t28 = sin(qJ(4));
	t31 = cos(qJ(4));
	t27 = sin(qJ(5));
	t30 = cos(qJ(5));
	t37 = t30 * r_i_i_C(1) - t27 * r_i_i_C(2) + pkin(4);
	t45 = pkin(8) + r_i_i_C(3);
	t33 = t28 * t45 + t37 * t31 + pkin(3);
	t23 = sin(pkin(9));
	t24 = sin(pkin(5));
	t44 = t23 * t24;
	t25 = cos(pkin(9));
	t43 = t25 * t24;
	t26 = cos(pkin(5));
	t42 = t26 * t32;
	t19 = -t32 * t22 - t29 * t40;
	t17 = t19 * t26;
	t7 = -t25 * t17 + t23 * t35;
	t38 = -t23 * t17 - t25 * t35;
	t36 = t27 * r_i_i_C(1) + t30 * r_i_i_C(2) + pkin(7);
	t34 = t35 * t26;
	t16 = t19 * t24;
	t15 = t35 * t24;
	t12 = -t16 * t31 + t26 * t28;
	t9 = t25 * t19 - t23 * t34;
	t6 = t23 * t19 + t25 * t34;
	t4 = t28 * t44 - t31 * t38;
	t2 = -t28 * t43 + t7 * t31;
	t1 = [0, (-t23 * t42 - t25 * t29) * pkin(2) - t36 * t38 + t33 * t9, t44, t45 * t4 + t37 * (t28 * t38 + t31 * t44), (-t4 * t27 - t9 * t30) * r_i_i_C(1) + (t9 * t27 - t4 * t30) * r_i_i_C(2); 0, (-t23 * t29 + t25 * t42) * pkin(2) + t36 * t7 + t33 * t6, -t43, t45 * t2 + t37 * (-t7 * t28 - t31 * t43), (-t2 * t27 - t6 * t30) * r_i_i_C(1) + (-t2 * t30 + t6 * t27) * r_i_i_C(2); 1, t24 * t32 * pkin(2) + t33 * t15 - t36 * t16, t26, t45 * t12 + t37 * (t16 * t28 + t26 * t31), (-t12 * t27 - t15 * t30) * r_i_i_C(1) + (-t12 * t30 + t15 * t27) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end