% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t1 = sin(pkin(6));
	t2 = [0, sin(pkin(12)) * t1, 0, 0, 0, 0; 0, -cos(pkin(12)) * t1, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (24->18), mult. (70->40), div. (0->0), fcn. (93->10), ass. (0->20)
	t6 = sin(pkin(12));
	t8 = sin(pkin(6));
	t21 = t6 * t8;
	t10 = cos(pkin(12));
	t20 = t10 * t8;
	t11 = cos(pkin(7));
	t9 = cos(pkin(13));
	t19 = t11 * t9;
	t12 = cos(pkin(6));
	t18 = t12 * t6;
	t17 = t10 * t12;
	t5 = sin(pkin(13));
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
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (102->33), mult. (288->67), div. (0->0), fcn. (377->12), ass. (0->34)
	t36 = pkin(9) + r_i_i_C(3);
	t15 = sin(pkin(12));
	t17 = sin(pkin(6));
	t35 = t15 * t17;
	t21 = cos(pkin(6));
	t34 = t15 * t21;
	t16 = sin(pkin(7));
	t33 = t16 * t17;
	t32 = t16 * t21;
	t18 = cos(pkin(13));
	t20 = cos(pkin(7));
	t31 = t18 * t20;
	t19 = cos(pkin(12));
	t30 = t19 * t17;
	t29 = t19 * t21;
	t22 = sin(qJ(4));
	t24 = cos(qJ(4));
	t28 = t24 * r_i_i_C(1) - t22 * r_i_i_C(2) + pkin(3);
	t14 = sin(pkin(13));
	t10 = -t15 * t14 + t18 * t29;
	t27 = t10 * t20 - t16 * t30;
	t12 = -t19 * t14 - t18 * t34;
	t26 = t12 * t20 + t15 * t33;
	t25 = cos(qJ(3));
	t23 = sin(qJ(3));
	t13 = -t14 * t34 + t19 * t18;
	t11 = t14 * t29 + t15 * t18;
	t9 = -t18 * t33 + t21 * t20;
	t8 = -t12 * t16 + t20 * t35;
	t7 = -t10 * t16 - t20 * t30;
	t6 = t23 * t32 + (t14 * t25 + t23 * t31) * t17;
	t4 = t13 * t25 + t26 * t23;
	t2 = t11 * t25 + t27 * t23;
	t1 = [0, t35, t36 * t4 + t28 * (-t13 * t23 + t26 * t25), (-t4 * t22 + t8 * t24) * r_i_i_C(1) + (-t8 * t22 - t4 * t24) * r_i_i_C(2), 0, 0; 0, -t30, t36 * t2 + t28 * (-t11 * t23 + t27 * t25), (-t2 * t22 + t7 * t24) * r_i_i_C(1) + (-t2 * t24 - t7 * t22) * r_i_i_C(2), 0, 0; 1, t21, t36 * t6 + t28 * (t25 * t32 + (-t14 * t23 + t25 * t31) * t17), (-t6 * t22 + t9 * t24) * r_i_i_C(1) + (-t9 * t22 - t6 * t24) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (198->41), mult. (459->77), div. (0->0), fcn. (599->14), ass. (0->40)
	t24 = sin(pkin(13));
	t25 = sin(pkin(12));
	t28 = cos(pkin(13));
	t29 = cos(pkin(12));
	t31 = cos(pkin(6));
	t40 = t29 * t31;
	t16 = -t25 * t24 + t28 * t40;
	t26 = sin(pkin(7));
	t30 = cos(pkin(7));
	t27 = sin(pkin(6));
	t41 = t29 * t27;
	t13 = -t16 * t26 - t30 * t41;
	t23 = qJ(4) + qJ(5);
	t21 = sin(t23);
	t22 = cos(t23);
	t17 = t24 * t40 + t25 * t28;
	t33 = sin(qJ(3));
	t35 = cos(qJ(3));
	t38 = t16 * t30 - t26 * t41;
	t8 = t17 * t35 + t38 * t33;
	t50 = (t13 * t22 - t8 * t21) * r_i_i_C(1) + (-t13 * t21 - t8 * t22) * r_i_i_C(2);
	t45 = t25 * t31;
	t19 = -t24 * t45 + t29 * t28;
	t18 = -t29 * t24 - t28 * t45;
	t44 = t26 * t27;
	t37 = t18 * t30 + t25 * t44;
	t10 = t19 * t35 + t37 * t33;
	t46 = t25 * t27;
	t14 = -t18 * t26 + t30 * t46;
	t49 = (-t10 * t21 + t14 * t22) * r_i_i_C(1) + (-t10 * t22 - t14 * t21) * r_i_i_C(2);
	t42 = t28 * t30;
	t43 = t26 * t31;
	t12 = t33 * t43 + (t24 * t35 + t33 * t42) * t27;
	t15 = -t28 * t44 + t31 * t30;
	t48 = (-t12 * t21 + t15 * t22) * r_i_i_C(1) + (-t12 * t22 - t15 * t21) * r_i_i_C(2);
	t47 = r_i_i_C(3) + pkin(10) + pkin(9);
	t34 = cos(qJ(4));
	t39 = t34 * pkin(4) + r_i_i_C(1) * t22 - r_i_i_C(2) * t21 + pkin(3);
	t32 = sin(qJ(4));
	t1 = [0, t46, t47 * t10 + t39 * (-t19 * t33 + t37 * t35), (-t10 * t32 + t14 * t34) * pkin(4) + t49, t49, 0; 0, -t41, t47 * t8 + t39 * (-t17 * t33 + t38 * t35), (t13 * t34 - t32 * t8) * pkin(4) + t50, t50, 0; 1, t31, t47 * t12 + t39 * (t35 * t43 + (-t24 * t33 + t35 * t42) * t27), (-t12 * t32 + t15 * t34) * pkin(4) + t48, t48, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:18:20
	% EndTime: 2019-10-09 21:18:20
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (498->55), mult. (1136->97), div. (0->0), fcn. (1492->16), ass. (0->52)
	t76 = pkin(11) + r_i_i_C(3);
	t38 = sin(qJ(6));
	t41 = cos(qJ(6));
	t75 = t41 * r_i_i_C(1) - t38 * r_i_i_C(2) + pkin(5);
	t37 = cos(pkin(12));
	t63 = sin(pkin(13));
	t64 = sin(pkin(12));
	t54 = t64 * t63;
	t66 = cos(pkin(13));
	t61 = t37 * t66;
	t68 = cos(pkin(6));
	t46 = -t68 * t61 + t54;
	t36 = sin(pkin(6));
	t65 = sin(pkin(7));
	t62 = t36 * t65;
	t67 = cos(pkin(7));
	t74 = t37 * t62 + t46 * t67;
	t55 = t64 * t66;
	t60 = t37 * t63;
	t47 = t68 * t55 + t60;
	t59 = t64 * t36;
	t73 = t47 * t67 - t65 * t59;
	t70 = cos(qJ(3));
	t69 = t37 * t36;
	t57 = t68 * t65;
	t56 = t67 * t66;
	t53 = t38 * r_i_i_C(1) + t41 * r_i_i_C(2) + pkin(9) + pkin(10);
	t27 = t68 * t60 + t55;
	t40 = sin(qJ(3));
	t17 = t27 * t70 - t74 * t40;
	t22 = t46 * t65 - t67 * t69;
	t35 = qJ(4) + qJ(5);
	t33 = sin(t35);
	t34 = cos(t35);
	t8 = t17 * t34 + t22 * t33;
	t51 = t76 * t8 + t75 * (-t17 * t33 + t22 * t34);
	t28 = -t68 * t54 + t61;
	t19 = t28 * t70 - t73 * t40;
	t23 = t47 * t65 + t67 * t59;
	t10 = t19 * t34 + t23 * t33;
	t50 = t76 * t10 + t75 * (-t19 * t33 + t23 * t34);
	t21 = t40 * t57 + (t40 * t56 + t70 * t63) * t36;
	t26 = -t66 * t62 + t68 * t67;
	t15 = t21 * t34 + t26 * t33;
	t49 = t76 * t15 + t75 * (-t21 * t33 + t26 * t34);
	t42 = cos(qJ(4));
	t48 = -t42 * pkin(4) - t76 * t33 - t75 * t34 - pkin(3);
	t39 = sin(qJ(4));
	t20 = -t70 * t57 + (t40 * t63 - t56 * t70) * t36;
	t18 = t28 * t40 + t73 * t70;
	t16 = t27 * t40 + t74 * t70;
	t1 = [0, t59, t48 * t18 + t53 * t19, (-t19 * t39 + t23 * t42) * pkin(4) + t50, t50, (-t10 * t38 + t18 * t41) * r_i_i_C(1) + (-t10 * t41 - t18 * t38) * r_i_i_C(2); 0, -t69, t48 * t16 + t53 * t17, (-t17 * t39 + t22 * t42) * pkin(4) + t51, t51, (t16 * t41 - t8 * t38) * r_i_i_C(1) + (-t16 * t38 - t8 * t41) * r_i_i_C(2); 1, t68, t48 * t20 + t53 * t21, (-t21 * t39 + t26 * t42) * pkin(4) + t49, t49, (-t15 * t38 + t20 * t41) * r_i_i_C(1) + (-t15 * t41 - t20 * t38) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end