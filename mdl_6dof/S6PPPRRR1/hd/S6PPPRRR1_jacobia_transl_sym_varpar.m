% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPPRRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPPRRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPPRRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobia_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t1 = sin(pkin(6));
	t2 = [0, sin(pkin(12)) * t1, 0, 0, 0, 0; 0, -cos(pkin(12)) * t1, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t2 = sin(pkin(12));
	t4 = sin(pkin(6));
	t11 = t2 * t4;
	t5 = cos(pkin(13));
	t8 = cos(pkin(6));
	t10 = t5 * t8;
	t6 = cos(pkin(12));
	t9 = t6 * t4;
	t7 = cos(pkin(7));
	t3 = sin(pkin(7));
	t1 = sin(pkin(13));
	t12 = [0, t11, -(-t6 * t1 - t2 * t10) * t3 + t7 * t11, 0, 0, 0; 0, -t9, -(-t2 * t1 + t6 * t10) * t3 - t7 * t9, 0, 0, 0; 1, t8, -t4 * t5 * t3 + t8 * t7, 0, 0, 0;];
	Ja_transl = t12;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (71->30), mult. (207->65), div. (0->0), fcn. (280->14), ass. (0->37)
	t16 = sin(pkin(12));
	t19 = sin(pkin(6));
	t39 = t16 * t19;
	t25 = cos(pkin(6));
	t38 = t16 * t25;
	t18 = sin(pkin(7));
	t37 = t18 * t19;
	t36 = t18 * t25;
	t21 = cos(pkin(13));
	t24 = cos(pkin(7));
	t35 = t21 * t24;
	t22 = cos(pkin(12));
	t34 = t22 * t19;
	t33 = t22 * t25;
	t15 = sin(pkin(13));
	t11 = t15 * t33 + t16 * t21;
	t14 = sin(pkin(14));
	t17 = sin(pkin(8));
	t20 = cos(pkin(14));
	t23 = cos(pkin(8));
	t10 = -t15 * t16 + t21 * t33;
	t29 = t10 * t24 - t18 * t34;
	t7 = -t10 * t18 - t24 * t34;
	t32 = (-t11 * t14 + t29 * t20) * t23 + t17 * t7;
	t13 = -t15 * t38 + t21 * t22;
	t12 = -t15 * t22 - t21 * t38;
	t28 = t12 * t24 + t16 * t37;
	t8 = -t12 * t18 + t24 * t39;
	t31 = t17 * t8 + t23 * (-t13 * t14 + t28 * t20);
	t9 = -t21 * t37 + t24 * t25;
	t30 = t17 * t9 + t23 * (t20 * t36 + (-t14 * t15 + t20 * t35) * t19);
	t27 = cos(qJ(4));
	t26 = sin(qJ(4));
	t6 = t15 * t19 * t20 + (t19 * t35 + t36) * t14;
	t4 = t13 * t20 + t28 * t14;
	t2 = t11 * t20 + t29 * t14;
	t1 = [0, t39, t8, (-t4 * t26 + t31 * t27) * r_i_i_C(1) + (-t31 * t26 - t27 * t4) * r_i_i_C(2), 0, 0; 0, -t34, t7, (-t2 * t26 + t32 * t27) * r_i_i_C(1) + (-t2 * t27 - t32 * t26) * r_i_i_C(2), 0, 0; 1, t25, t9, (-t6 * t26 + t30 * t27) * r_i_i_C(1) + (-t30 * t26 - t27 * t6) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (280->45), mult. (804->91), div. (0->0), fcn. (1074->16), ass. (0->50)
	t52 = pkin(10) + r_i_i_C(3);
	t25 = sin(pkin(12));
	t28 = sin(pkin(6));
	t51 = t25 * t28;
	t34 = cos(pkin(6));
	t50 = t25 * t34;
	t27 = sin(pkin(7));
	t49 = t27 * t28;
	t48 = t27 * t34;
	t30 = cos(pkin(13));
	t33 = cos(pkin(7));
	t47 = t30 * t33;
	t31 = cos(pkin(12));
	t46 = t31 * t28;
	t45 = t31 * t34;
	t24 = sin(pkin(13));
	t19 = -t25 * t24 + t30 * t45;
	t16 = -t19 * t27 - t33 * t46;
	t26 = sin(pkin(8));
	t32 = cos(pkin(8));
	t20 = t24 * t45 + t25 * t30;
	t23 = sin(pkin(14));
	t29 = cos(pkin(14));
	t40 = t19 * t33 - t27 * t46;
	t9 = -t20 * t23 + t40 * t29;
	t44 = t16 * t26 + t32 * t9;
	t22 = -t24 * t50 + t31 * t30;
	t21 = -t31 * t24 - t30 * t50;
	t39 = t21 * t33 + t25 * t49;
	t11 = -t22 * t23 + t39 * t29;
	t17 = -t21 * t27 + t33 * t51;
	t43 = t11 * t32 + t17 * t26;
	t14 = t29 * t48 + (-t23 * t24 + t29 * t47) * t28;
	t18 = -t30 * t49 + t34 * t33;
	t42 = t14 * t32 + t18 * t26;
	t35 = sin(qJ(5));
	t37 = cos(qJ(5));
	t41 = t37 * r_i_i_C(1) - t35 * r_i_i_C(2) + pkin(4);
	t38 = cos(qJ(4));
	t36 = sin(qJ(4));
	t15 = t28 * t24 * t29 + (t28 * t47 + t48) * t23;
	t13 = -t14 * t26 + t18 * t32;
	t12 = t22 * t29 + t39 * t23;
	t10 = t20 * t29 + t40 * t23;
	t8 = -t11 * t26 + t17 * t32;
	t7 = t16 * t32 - t9 * t26;
	t6 = t15 * t38 + t42 * t36;
	t4 = t12 * t38 + t43 * t36;
	t2 = t10 * t38 + t44 * t36;
	t1 = [0, t51, t17, t52 * t4 + t41 * (-t12 * t36 + t43 * t38), (-t4 * t35 + t8 * t37) * r_i_i_C(1) + (-t8 * t35 - t4 * t37) * r_i_i_C(2), 0; 0, -t46, t16, t52 * t2 + t41 * (-t10 * t36 + t44 * t38), (-t2 * t35 + t7 * t37) * r_i_i_C(1) + (-t2 * t37 - t7 * t35) * r_i_i_C(2), 0; 1, t34, t18, t52 * t6 + t41 * (-t15 * t36 + t42 * t38), (t13 * t37 - t6 * t35) * r_i_i_C(1) + (-t13 * t35 - t6 * t37) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (756->60), mult. (2166->115), div. (0->0), fcn. (2894->18), ass. (0->66)
	t65 = sin(pkin(13));
	t66 = sin(pkin(12));
	t54 = t66 * t65;
	t71 = cos(pkin(13));
	t72 = cos(pkin(12));
	t63 = t72 * t71;
	t75 = cos(pkin(6));
	t44 = t75 * t63 - t54;
	t74 = cos(pkin(7));
	t41 = t44 * t74;
	t56 = t66 * t71;
	t61 = t72 * t65;
	t45 = t75 * t61 + t56;
	t68 = sin(pkin(7));
	t70 = cos(pkin(14));
	t58 = t68 * t70;
	t69 = sin(pkin(6));
	t50 = t69 * t58;
	t64 = sin(pkin(14));
	t32 = -t70 * t41 + t45 * t64 + t72 * t50;
	t62 = t72 * t69;
	t38 = -t44 * t68 - t74 * t62;
	t67 = sin(pkin(8));
	t73 = cos(pkin(8));
	t80 = t32 * t73 - t38 * t67;
	t46 = -t75 * t56 - t61;
	t42 = t46 * t74;
	t47 = -t75 * t54 + t63;
	t33 = -t70 * t42 + t47 * t64 - t66 * t50;
	t55 = t66 * t69;
	t39 = -t46 * t68 + t74 * t55;
	t79 = t33 * t73 - t39 * t67;
	t60 = t69 * t71;
	t51 = t74 * t60;
	t59 = t69 * t65;
	t37 = -t70 * t51 - t75 * t58 + t64 * t59;
	t43 = -t68 * t60 + t75 * t74;
	t78 = t37 * t73 - t43 * t67;
	t77 = pkin(11) + r_i_i_C(3);
	t76 = cos(qJ(4));
	t57 = t68 * t64;
	t25 = sin(qJ(6));
	t28 = cos(qJ(6));
	t53 = t28 * r_i_i_C(1) - t25 * r_i_i_C(2) + pkin(5);
	t52 = t25 * r_i_i_C(1) + t28 * r_i_i_C(2) + pkin(10);
	t49 = t69 * t57;
	t26 = sin(qJ(5));
	t29 = cos(qJ(5));
	t48 = -t77 * t26 - t53 * t29 - pkin(4);
	t27 = sin(qJ(4));
	t23 = t64 * t51 + t75 * t57 + t70 * t59;
	t19 = t37 * t67 + t43 * t73;
	t18 = t64 * t42 + t47 * t70 + t66 * t49;
	t17 = t64 * t41 + t45 * t70 - t72 * t49;
	t14 = t33 * t67 + t39 * t73;
	t13 = t32 * t67 + t38 * t73;
	t12 = t23 * t76 - t78 * t27;
	t11 = t23 * t27 + t78 * t76;
	t10 = t18 * t76 - t79 * t27;
	t9 = t18 * t27 + t79 * t76;
	t8 = t17 * t76 - t80 * t27;
	t7 = t17 * t27 + t80 * t76;
	t6 = t12 * t29 + t19 * t26;
	t4 = t10 * t29 + t14 * t26;
	t2 = t13 * t26 + t8 * t29;
	t1 = [0, t55, t39, t52 * t10 + t48 * t9, t77 * t4 + t53 * (-t10 * t26 + t14 * t29), (-t4 * t25 + t9 * t28) * r_i_i_C(1) + (-t9 * t25 - t4 * t28) * r_i_i_C(2); 0, -t62, t38, t48 * t7 + t52 * t8, t77 * t2 + t53 * (t13 * t29 - t8 * t26), (-t2 * t25 + t7 * t28) * r_i_i_C(1) + (-t2 * t28 - t7 * t25) * r_i_i_C(2); 1, t75, t43, t48 * t11 + t52 * t12, t77 * t6 + t53 * (-t12 * t26 + t19 * t29), (t11 * t28 - t6 * t25) * r_i_i_C(1) + (-t11 * t25 - t6 * t28) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end