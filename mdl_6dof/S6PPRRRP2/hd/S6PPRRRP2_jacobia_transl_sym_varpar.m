% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:16
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRRP2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRP2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRP2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t1 = sin(pkin(6));
	t2 = [0, sin(pkin(11)) * t1, 0, 0, 0, 0; 0, -cos(pkin(11)) * t1, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
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
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (102->33), mult. (288->67), div. (0->0), fcn. (377->12), ass. (0->34)
	t36 = pkin(9) + r_i_i_C(3);
	t15 = sin(pkin(11));
	t17 = sin(pkin(6));
	t35 = t15 * t17;
	t21 = cos(pkin(6));
	t34 = t15 * t21;
	t16 = sin(pkin(7));
	t33 = t16 * t17;
	t32 = t16 * t21;
	t18 = cos(pkin(12));
	t20 = cos(pkin(7));
	t31 = t18 * t20;
	t19 = cos(pkin(11));
	t30 = t19 * t17;
	t29 = t19 * t21;
	t22 = sin(qJ(4));
	t24 = cos(qJ(4));
	t28 = t24 * r_i_i_C(1) - t22 * r_i_i_C(2) + pkin(3);
	t14 = sin(pkin(12));
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
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (282->46), mult. (794->87), div. (0->0), fcn. (1045->14), ass. (0->46)
	t24 = cos(pkin(11));
	t47 = sin(pkin(12));
	t48 = sin(pkin(11));
	t38 = t48 * t47;
	t50 = cos(pkin(12));
	t45 = t24 * t50;
	t52 = cos(pkin(6));
	t32 = -t52 * t45 + t38;
	t23 = sin(pkin(6));
	t49 = sin(pkin(7));
	t46 = t23 * t49;
	t51 = cos(pkin(7));
	t57 = t24 * t46 + t32 * t51;
	t39 = t48 * t50;
	t44 = t24 * t47;
	t33 = t52 * t39 + t44;
	t43 = t48 * t23;
	t56 = t33 * t51 - t49 * t43;
	t55 = pkin(10) + r_i_i_C(3);
	t54 = cos(qJ(3));
	t53 = t24 * t23;
	t41 = t52 * t49;
	t40 = t51 * t50;
	t25 = sin(qJ(5));
	t28 = cos(qJ(5));
	t37 = t28 * r_i_i_C(1) - t25 * r_i_i_C(2) + pkin(4);
	t36 = t25 * r_i_i_C(1) + t28 * r_i_i_C(2) + pkin(9);
	t26 = sin(qJ(4));
	t29 = cos(qJ(4));
	t34 = -t55 * t26 - t37 * t29 - pkin(3);
	t27 = sin(qJ(3));
	t19 = -t52 * t38 + t45;
	t18 = t52 * t44 + t39;
	t17 = -t50 * t46 + t52 * t51;
	t14 = t33 * t49 + t51 * t43;
	t13 = t32 * t49 - t51 * t53;
	t12 = t27 * t41 + (t27 * t40 + t54 * t47) * t23;
	t11 = -t54 * t41 + (t27 * t47 - t40 * t54) * t23;
	t10 = t12 * t29 + t17 * t26;
	t8 = t19 * t54 - t56 * t27;
	t7 = t19 * t27 + t56 * t54;
	t6 = t18 * t54 - t57 * t27;
	t5 = t18 * t27 + t57 * t54;
	t4 = t14 * t26 + t8 * t29;
	t2 = t13 * t26 + t6 * t29;
	t1 = [0, t43, t34 * t7 + t36 * t8, t55 * t4 + t37 * (t14 * t29 - t8 * t26), (-t4 * t25 + t7 * t28) * r_i_i_C(1) + (-t7 * t25 - t4 * t28) * r_i_i_C(2), 0; 0, -t53, t34 * t5 + t36 * t6, t55 * t2 + t37 * (t13 * t29 - t6 * t26), (-t2 * t25 + t5 * t28) * r_i_i_C(1) + (-t2 * t28 - t5 * t25) * r_i_i_C(2), 0; 1, t52, t34 * t11 + t36 * t12, t55 * t10 + t37 * (-t12 * t26 + t17 * t29), (-t10 * t25 + t11 * t28) * r_i_i_C(1) + (-t10 * t28 - t11 * t25) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:16:21
	% EndTime: 2019-10-09 21:16:21
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (475->58), mult. (1331->105), div. (0->0), fcn. (1760->14), ass. (0->52)
	t39 = cos(pkin(11));
	t61 = sin(pkin(12));
	t62 = sin(pkin(11));
	t52 = t62 * t61;
	t64 = cos(pkin(12));
	t59 = t39 * t64;
	t66 = cos(pkin(6));
	t47 = -t66 * t59 + t52;
	t38 = sin(pkin(6));
	t63 = sin(pkin(7));
	t60 = t38 * t63;
	t65 = cos(pkin(7));
	t75 = t39 * t60 + t47 * t65;
	t53 = t62 * t64;
	t58 = t39 * t61;
	t48 = t66 * t53 + t58;
	t57 = t62 * t38;
	t74 = t48 * t65 - t63 * t57;
	t73 = pkin(5) + r_i_i_C(1);
	t72 = pkin(10) + r_i_i_C(2);
	t71 = cos(qJ(3));
	t70 = t39 * t38;
	t40 = sin(qJ(5));
	t44 = cos(qJ(4));
	t69 = t40 * t44;
	t43 = cos(qJ(5));
	t68 = t43 * t44;
	t67 = r_i_i_C(3) + qJ(6);
	t55 = t66 * t63;
	t54 = t65 * t64;
	t41 = sin(qJ(4));
	t50 = -pkin(4) * t44 - t72 * t41 - pkin(3);
	t49 = t67 * t40 + t73 * t43 + pkin(4);
	t42 = sin(qJ(3));
	t34 = -t66 * t52 + t59;
	t33 = t66 * t58 + t53;
	t32 = -t64 * t60 + t66 * t65;
	t29 = t48 * t63 + t65 * t57;
	t28 = t47 * t63 - t65 * t70;
	t27 = t42 * t55 + (t42 * t54 + t71 * t61) * t38;
	t26 = -t71 * t55 + (t42 * t61 - t54 * t71) * t38;
	t24 = t27 * t44 + t32 * t41;
	t22 = t34 * t71 - t74 * t42;
	t21 = t34 * t42 + t74 * t71;
	t20 = t33 * t71 - t75 * t42;
	t19 = t33 * t42 + t75 * t71;
	t13 = t24 * t40 - t26 * t43;
	t12 = t22 * t44 + t29 * t41;
	t10 = t20 * t44 + t28 * t41;
	t3 = t12 * t40 - t21 * t43;
	t1 = t10 * t40 - t19 * t43;
	t2 = [0, t57, t22 * pkin(9) + t73 * (-t21 * t68 + t22 * t40) + t67 * (-t21 * t69 - t22 * t43) + t50 * t21, t72 * t12 + t49 * (-t22 * t41 + t29 * t44), t67 * (t12 * t43 + t21 * t40) - t73 * t3, t3; 0, -t70, t20 * pkin(9) + t73 * (-t19 * t68 + t20 * t40) + t67 * (-t19 * t69 - t20 * t43) + t50 * t19, t72 * t10 + t49 * (-t20 * t41 + t28 * t44), t67 * (t10 * t43 + t19 * t40) - t73 * t1, t1; 1, t66, t27 * pkin(9) + t73 * (-t26 * t68 + t27 * t40) + t67 * (-t26 * t69 - t27 * t43) + t50 * t26, t72 * t24 + t49 * (-t27 * t41 + t32 * t44), t67 * (t24 * t43 + t26 * t40) - t73 * t13, t13;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end