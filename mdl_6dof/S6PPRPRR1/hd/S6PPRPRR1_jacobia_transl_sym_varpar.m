% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:08
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRPRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRPRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRPRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t1 = sin(pkin(6));
	t2 = [0, sin(pkin(11)) * t1, 0, 0, 0, 0; 0, -cos(pkin(11)) * t1, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (24->18), mult. (70->40), div. (0->0), fcn. (93->10), ass. (0->20)
	t6 = sin(pkin(11));
	t8 = sin(pkin(6));
	t21 = t6 * t8;
	t10 = cos(pkin(11));
	t20 = t10 * t8;
	t11 = cos(pkin(7));
	t9 = cos(pkin(12));
	t19 = t11 * t9;
	t12 = cos(pkin(6));
	t18 = t12 * t6;
	t17 = t10 * t12;
	t5 = sin(pkin(12));
	t7 = sin(pkin(7));
	t16 = t11 * (-t10 * t5 - t9 * t18) + t7 * t21;
	t15 = -(t9 * t17 - t5 * t6) * t11 + t7 * t20;
	t14 = cos(qJ(3));
	t13 = sin(qJ(3));
	t4 = t10 * t9 - t5 * t18;
	t2 = t5 * t17 + t6 * t9;
	t1 = [0, t21, (-t4 * t13 + t16 * t14) * r_i_i_C(1) + (-t16 * t13 - t14 * t4) * r_i_i_C(2), 0, 0, 0; 0, -t20, (-t2 * t13 - t15 * t14) * r_i_i_C(1) + (t15 * t13 - t14 * t2) * r_i_i_C(2), 0, 0, 0; 1, t12, (r_i_i_C(1) * t14 - r_i_i_C(2) * t13) * t7 * t12 + ((-t13 * t5 + t14 * t19) * r_i_i_C(1) + (-t13 * t19 - t14 * t5) * r_i_i_C(2)) * t8, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (60->34), mult. (173->69), div. (0->0), fcn. (226->12), ass. (0->28)
	t13 = sin(pkin(11));
	t15 = sin(pkin(6));
	t27 = t13 * t15;
	t20 = cos(pkin(6));
	t26 = t13 * t20;
	t14 = sin(pkin(7));
	t25 = t14 * t15;
	t18 = cos(pkin(11));
	t24 = t18 * t15;
	t23 = t18 * t20;
	t11 = sin(pkin(13));
	t16 = cos(pkin(13));
	t21 = sin(qJ(3));
	t22 = cos(qJ(3));
	t10 = -t11 * t22 - t21 * t16;
	t9 = t21 * t11 - t16 * t22;
	t19 = cos(pkin(7));
	t17 = cos(pkin(12));
	t12 = sin(pkin(12));
	t8 = -t12 * t26 + t17 * t18;
	t7 = -t12 * t18 - t17 * t26;
	t6 = t12 * t23 + t13 * t17;
	t5 = -t12 * t13 + t17 * t23;
	t4 = t10 * t19;
	t3 = t9 * t19;
	t2 = t10 * t14;
	t1 = t9 * t14;
	t28 = [0, t27, (-t1 * t27 + t10 * t8 - t3 * t7) * r_i_i_C(1) + (t2 * t27 + t4 * t7 + t8 * t9) * r_i_i_C(2) + (-t8 * t21 + (t13 * t25 + t19 * t7) * t22) * pkin(3), -t14 * t7 + t19 * t27, 0, 0; 0, -t24, (t1 * t24 + t10 * t6 - t3 * t5) * r_i_i_C(1) + (-t2 * t24 + t4 * t5 + t6 * t9) * r_i_i_C(2) + (-t6 * t21 + (-t14 * t24 + t19 * t5) * t22) * pkin(3), -t14 * t5 - t19 * t24, 0, 0; 1, t20, (pkin(3) * t14 * t22 - r_i_i_C(1) * t1 + r_i_i_C(2) * t2) * t20 + ((t10 * t12 - t17 * t3) * r_i_i_C(1) + (t12 * t9 + t17 * t4) * r_i_i_C(2) + (t17 * t19 * t22 - t12 * t21) * pkin(3)) * t15, -t17 * t25 + t19 * t20, 0, 0;];
	Ja_transl = t28;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (183->46), mult. (511->90), div. (0->0), fcn. (675->14), ass. (0->38)
	t27 = sin(pkin(7));
	t24 = sin(pkin(13));
	t29 = cos(pkin(13));
	t35 = sin(qJ(3));
	t37 = cos(qJ(3));
	t40 = t37 * t24 + t35 * t29;
	t13 = t40 * t27;
	t32 = cos(pkin(7));
	t15 = t40 * t32;
	t21 = t35 * t24 - t37 * t29;
	t25 = sin(pkin(12));
	t28 = sin(pkin(6));
	t30 = cos(pkin(12));
	t33 = cos(pkin(6));
	t9 = t33 * t13 + (t15 * t30 - t21 * t25) * t28;
	t48 = -pkin(9) - r_i_i_C(3);
	t26 = sin(pkin(11));
	t47 = t26 * t28;
	t46 = t26 * t33;
	t45 = t27 * t28;
	t31 = cos(pkin(11));
	t44 = t31 * t28;
	t43 = t31 * t33;
	t34 = sin(qJ(5));
	t36 = cos(qJ(5));
	t39 = t36 * r_i_i_C(1) - t34 * r_i_i_C(2) + pkin(4);
	t17 = -t26 * t25 + t30 * t43;
	t18 = t25 * t43 + t26 * t30;
	t38 = t13 * t44 - t17 * t15 + t18 * t21;
	t19 = -t31 * t25 - t30 * t46;
	t20 = -t25 * t46 + t31 * t30;
	t6 = t13 * t47 + t19 * t15 - t20 * t21;
	t16 = -t30 * t45 + t33 * t32;
	t14 = t21 * t32;
	t12 = t21 * t27;
	t11 = -t19 * t27 + t32 * t47;
	t10 = -t17 * t27 - t32 * t44;
	t1 = [0, t47, -t48 * t6 + t39 * (-t12 * t47 - t19 * t14 - t20 * t40) + (-t20 * t35 + (t19 * t32 + t26 * t45) * t37) * pkin(3), t11, (t11 * t36 - t6 * t34) * r_i_i_C(1) + (-t11 * t34 - t6 * t36) * r_i_i_C(2), 0; 0, -t44, t48 * t38 + t39 * (t12 * t44 - t17 * t14 - t18 * t40) + (-t18 * t35 + (t17 * t32 - t27 * t44) * t37) * pkin(3), t10, (t10 * t36 + t34 * t38) * r_i_i_C(1) + (-t10 * t34 + t36 * t38) * r_i_i_C(2), 0; 1, t33, -t48 * t9 + t39 * (-t33 * t12 + (-t14 * t30 - t25 * t40) * t28) + (t27 * t33 * t37 + (t30 * t32 * t37 - t25 * t35) * t28) * pkin(3), t16, (t16 * t36 - t9 * t34) * r_i_i_C(1) + (-t16 * t34 - t9 * t36) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (462->59), mult. (1281->112), div. (0->0), fcn. (1706->16), ass. (0->48)
	t35 = sin(pkin(7));
	t32 = sin(pkin(13));
	t37 = cos(pkin(13));
	t44 = sin(qJ(3));
	t47 = cos(qJ(3));
	t52 = t47 * t32 + t44 * t37;
	t21 = t52 * t35;
	t40 = cos(pkin(7));
	t23 = t52 * t40;
	t29 = t44 * t32 - t47 * t37;
	t33 = sin(pkin(12));
	t36 = sin(pkin(6));
	t38 = cos(pkin(12));
	t41 = cos(pkin(6));
	t15 = t41 * t21 + (t23 * t38 - t29 * t33) * t36;
	t60 = pkin(10) + r_i_i_C(3);
	t34 = sin(pkin(11));
	t59 = t34 * t36;
	t58 = t34 * t41;
	t57 = t35 * t36;
	t39 = cos(pkin(11));
	t56 = t39 * t36;
	t55 = t39 * t41;
	t42 = sin(qJ(6));
	t45 = cos(qJ(6));
	t51 = t45 * r_i_i_C(1) - t42 * r_i_i_C(2) + pkin(5);
	t50 = -t42 * r_i_i_C(1) - t45 * r_i_i_C(2) - pkin(9);
	t25 = -t34 * t33 + t38 * t55;
	t26 = t33 * t55 + t34 * t38;
	t49 = t21 * t56 - t25 * t23 + t26 * t29;
	t27 = -t39 * t33 - t38 * t58;
	t28 = -t33 * t58 + t39 * t38;
	t10 = t21 * t59 + t27 * t23 - t28 * t29;
	t43 = sin(qJ(5));
	t46 = cos(qJ(5));
	t48 = t60 * t43 + t51 * t46 + pkin(4);
	t24 = -t38 * t57 + t41 * t40;
	t22 = t29 * t40;
	t20 = t29 * t35;
	t17 = -t27 * t35 + t40 * t59;
	t16 = -t25 * t35 - t40 * t56;
	t14 = -t41 * t20 + (-t22 * t38 - t33 * t52) * t36;
	t12 = t15 * t46 + t24 * t43;
	t9 = -t20 * t59 - t27 * t22 - t28 * t52;
	t6 = t20 * t56 - t25 * t22 - t26 * t52;
	t4 = t10 * t46 + t17 * t43;
	t2 = t16 * t43 - t46 * t49;
	t1 = [0, t59, -t50 * t10 + (-t28 * t44 + (t27 * t40 + t34 * t57) * t47) * pkin(3) + t48 * t9, t17, t60 * t4 + t51 * (-t10 * t43 + t17 * t46), (-t4 * t42 - t9 * t45) * r_i_i_C(1) + (-t4 * t45 + t9 * t42) * r_i_i_C(2); 0, -t56, t50 * t49 + (-t26 * t44 + (t25 * t40 - t35 * t56) * t47) * pkin(3) + t48 * t6, t16, t60 * t2 + t51 * (t16 * t46 + t43 * t49), (-t2 * t42 - t6 * t45) * r_i_i_C(1) + (-t2 * t45 + t6 * t42) * r_i_i_C(2); 1, t41, -t50 * t15 + (t41 * t35 * t47 + (t38 * t40 * t47 - t33 * t44) * t36) * pkin(3) + t48 * t14, t24, t60 * t12 + t51 * (-t15 * t43 + t24 * t46), (-t12 * t42 - t14 * t45) * r_i_i_C(1) + (-t12 * t45 + t14 * t42) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end