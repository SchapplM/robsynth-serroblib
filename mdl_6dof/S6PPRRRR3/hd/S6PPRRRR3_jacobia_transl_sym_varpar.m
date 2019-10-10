% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRRR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobia_transl_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t1 = sin(pkin(6));
	t2 = [0, sin(pkin(13)) * t1, 0, 0, 0, 0; 0, -cos(pkin(13)) * t1, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (24->18), mult. (70->40), div. (0->0), fcn. (93->10), ass. (0->20)
	t6 = sin(pkin(13));
	t8 = sin(pkin(6));
	t21 = t6 * t8;
	t10 = cos(pkin(13));
	t20 = t10 * t8;
	t11 = cos(pkin(7));
	t9 = cos(pkin(14));
	t19 = t11 * t9;
	t12 = cos(pkin(6));
	t18 = t12 * t6;
	t17 = t10 * t12;
	t5 = sin(pkin(14));
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
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (154->47), mult. (448->94), div. (0->0), fcn. (589->14), ass. (0->40)
	t16 = sin(pkin(13));
	t19 = sin(pkin(6));
	t43 = t16 * t19;
	t24 = cos(pkin(6));
	t42 = t16 * t24;
	t18 = sin(pkin(7));
	t41 = t18 * t19;
	t40 = t18 * t24;
	t21 = cos(pkin(13));
	t39 = t21 * t19;
	t38 = t21 * t24;
	t22 = cos(pkin(8));
	t25 = sin(qJ(4));
	t37 = t22 * t25;
	t27 = cos(qJ(4));
	t36 = t22 * t27;
	t23 = cos(pkin(7));
	t26 = sin(qJ(3));
	t35 = t23 * t26;
	t34 = t18 * t39;
	t17 = sin(pkin(8));
	t33 = t17 * (pkin(10) + r_i_i_C(3));
	t15 = sin(pkin(14));
	t20 = cos(pkin(14));
	t10 = -t16 * t15 + t20 * t38;
	t11 = t15 * t38 + t16 * t20;
	t28 = cos(qJ(3));
	t1 = -t11 * t26 + (t10 * t23 - t34) * t28;
	t32 = t1 * t22 + t17 * (-t10 * t18 - t23 * t39);
	t12 = -t21 * t15 - t20 * t42;
	t13 = -t15 * t42 + t21 * t20;
	t29 = t12 * t23 + t16 * t41;
	t3 = -t13 * t26 + t29 * t28;
	t31 = t17 * (-t12 * t18 + t23 * t43) + t22 * t3;
	t5 = t28 * t40 + (t20 * t23 * t28 - t15 * t26) * t19;
	t30 = t17 * (-t20 * t41 + t24 * t23) + t22 * t5;
	t6 = t26 * t40 + (t15 * t28 + t20 * t35) * t19;
	t4 = t13 * t28 + t29 * t26;
	t2 = t10 * t35 + t11 * t28 - t26 * t34;
	t7 = [0, t43, (t3 * t27 - t4 * t37) * r_i_i_C(1) + (-t3 * t25 - t4 * t36) * r_i_i_C(2) + t3 * pkin(3) + t4 * t33, (-t4 * t25 + t31 * t27) * r_i_i_C(1) + (-t31 * t25 - t4 * t27) * r_i_i_C(2), 0, 0; 0, -t39, (t1 * t27 - t2 * t37) * r_i_i_C(1) + (-t1 * t25 - t2 * t36) * r_i_i_C(2) + t1 * pkin(3) + t2 * t33, (-t2 * t25 + t32 * t27) * r_i_i_C(1) + (-t2 * t27 - t32 * t25) * r_i_i_C(2), 0, 0; 1, t24, (t5 * t27 - t6 * t37) * r_i_i_C(1) + (-t5 * t25 - t6 * t36) * r_i_i_C(2) + t5 * pkin(3) + t6 * t33, (-t6 * t25 + t30 * t27) * r_i_i_C(1) + (-t30 * t25 - t6 * t27) * r_i_i_C(2), 0, 0;];
	Ja_transl = t7;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (454->73), mult. (1310->140), div. (0->0), fcn. (1733->16), ass. (0->58)
	t64 = r_i_i_C(3) + pkin(11);
	t32 = sin(pkin(8));
	t63 = t32 * pkin(10);
	t31 = sin(pkin(13));
	t34 = sin(pkin(6));
	t62 = t31 * t34;
	t39 = cos(pkin(6));
	t61 = t31 * t39;
	t40 = sin(qJ(5));
	t60 = t32 * t40;
	t43 = cos(qJ(5));
	t59 = t32 * t43;
	t33 = sin(pkin(7));
	t58 = t33 * t34;
	t57 = t33 * t39;
	t36 = cos(pkin(13));
	t56 = t36 * t34;
	t55 = t36 * t39;
	t37 = cos(pkin(8));
	t41 = sin(qJ(4));
	t54 = t37 * t41;
	t44 = cos(qJ(4));
	t53 = t37 * t44;
	t38 = cos(pkin(7));
	t42 = sin(qJ(3));
	t52 = t38 * t42;
	t51 = t33 * t56;
	t30 = sin(pkin(14));
	t35 = cos(pkin(14));
	t25 = -t31 * t30 + t35 * t55;
	t26 = t30 * t55 + t31 * t35;
	t45 = cos(qJ(3));
	t15 = -t26 * t42 + (t25 * t38 - t51) * t45;
	t22 = -t25 * t33 - t38 * t56;
	t50 = t15 * t37 + t22 * t32;
	t28 = -t30 * t61 + t36 * t35;
	t27 = -t36 * t30 - t35 * t61;
	t46 = t27 * t38 + t31 * t58;
	t17 = -t28 * t42 + t46 * t45;
	t23 = -t27 * t33 + t38 * t62;
	t49 = t17 * t37 + t23 * t32;
	t20 = t45 * t57 + (t35 * t38 * t45 - t30 * t42) * t34;
	t24 = -t35 * t58 + t39 * t38;
	t48 = t20 * t37 + t24 * t32;
	t47 = t43 * r_i_i_C(1) - t40 * r_i_i_C(2) + pkin(4);
	t21 = t42 * t57 + (t30 * t45 + t35 * t52) * t34;
	t19 = -t20 * t32 + t24 * t37;
	t18 = t28 * t45 + t46 * t42;
	t16 = t25 * t52 + t26 * t45 - t42 * t51;
	t14 = t20 * t44 - t21 * t54;
	t12 = -t17 * t32 + t23 * t37;
	t11 = -t15 * t32 + t22 * t37;
	t10 = t21 * t44 + t48 * t41;
	t8 = t17 * t44 - t18 * t54;
	t6 = t15 * t44 - t16 * t54;
	t4 = t18 * t44 + t49 * t41;
	t2 = t16 * t44 + t50 * t41;
	t1 = [0, t62, (t18 * t60 + t8 * t43) * r_i_i_C(1) + (t18 * t59 - t8 * t40) * r_i_i_C(2) + t8 * pkin(4) + t17 * pkin(3) + t18 * t63 + t64 * (t17 * t41 + t18 * t53), t64 * t4 + t47 * (-t18 * t41 + t49 * t44), (t12 * t43 - t4 * t40) * r_i_i_C(1) + (-t12 * t40 - t4 * t43) * r_i_i_C(2), 0; 0, -t56, (t16 * t60 + t6 * t43) * r_i_i_C(1) + (t16 * t59 - t6 * t40) * r_i_i_C(2) + t6 * pkin(4) + t15 * pkin(3) + t16 * t63 + t64 * (t15 * t41 + t16 * t53), t64 * t2 + t47 * (-t16 * t41 + t50 * t44), (t11 * t43 - t2 * t40) * r_i_i_C(1) + (-t11 * t40 - t2 * t43) * r_i_i_C(2), 0; 1, t39, (t14 * t43 + t21 * t60) * r_i_i_C(1) + (-t14 * t40 + t21 * t59) * r_i_i_C(2) + t14 * pkin(4) + t20 * pkin(3) + t21 * t63 + t64 * (t20 * t41 + t21 * t53), t64 * t10 + t47 * (-t21 * t41 + t48 * t44), (-t10 * t40 + t19 * t43) * r_i_i_C(1) + (-t10 * t43 - t19 * t40) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.64s
	% Computational Cost: add. (1073->96), mult. (3094->175), div. (0->0), fcn. (4112->18), ass. (0->74)
	t90 = r_i_i_C(3) + pkin(12);
	t42 = sin(pkin(8));
	t41 = sin(pkin(13));
	t43 = cos(pkin(14));
	t44 = cos(pkin(6));
	t40 = sin(pkin(14));
	t82 = cos(pkin(13));
	t77 = t82 * t40;
	t37 = t41 * t43 + t44 * t77;
	t48 = sin(qJ(3));
	t51 = cos(qJ(3));
	t76 = t82 * t43;
	t70 = -t41 * t40 + t44 * t76;
	t80 = sin(pkin(7));
	t81 = sin(pkin(6));
	t73 = t80 * t81;
	t84 = cos(pkin(7));
	t91 = t70 * t84 - t82 * t73;
	t54 = t37 * t48 - t51 * t91;
	t74 = t81 * t84;
	t58 = -t70 * t80 - t82 * t74;
	t83 = cos(pkin(8));
	t94 = -t58 * t42 + t54 * t83;
	t87 = t41 * t44;
	t38 = -t40 * t87 + t76;
	t69 = -t43 * t87 - t77;
	t62 = t41 * t73 + t69 * t84;
	t55 = t38 * t48 - t62 * t51;
	t61 = t41 * t74 - t69 * t80;
	t93 = -t61 * t42 + t55 * t83;
	t65 = t43 * t74 + t80 * t44;
	t79 = t40 * t81;
	t60 = t48 * t79 - t65 * t51;
	t66 = -t43 * t73 + t44 * t84;
	t92 = -t66 * t42 + t60 * t83;
	t89 = pkin(10) * t42;
	t88 = cos(qJ(4));
	t46 = sin(qJ(5));
	t86 = t42 * t46;
	t50 = cos(qJ(5));
	t85 = t42 * t50;
	t47 = sin(qJ(4));
	t78 = t47 * t83;
	t75 = t83 * t88;
	t45 = sin(qJ(6));
	t49 = cos(qJ(6));
	t72 = t49 * r_i_i_C(1) - t45 * r_i_i_C(2) + pkin(5);
	t71 = t45 * r_i_i_C(1) + t49 * r_i_i_C(2) + pkin(11);
	t64 = -t90 * t46 - t72 * t50 - pkin(4);
	t35 = t65 * t48 + t51 * t79;
	t31 = t60 * t42 + t66 * t83;
	t30 = t38 * t51 + t62 * t48;
	t29 = t37 * t51 + t48 * t91;
	t26 = -t35 * t78 - t60 * t88;
	t25 = t35 * t75 - t60 * t47;
	t24 = t55 * t42 + t61 * t83;
	t23 = t54 * t42 + t58 * t83;
	t22 = t35 * t88 - t47 * t92;
	t21 = t35 * t47 + t88 * t92;
	t20 = t26 * t50 + t35 * t86;
	t18 = -t30 * t78 - t55 * t88;
	t17 = t30 * t75 - t55 * t47;
	t16 = -t29 * t78 - t54 * t88;
	t15 = t29 * t75 - t54 * t47;
	t14 = t22 * t50 + t31 * t46;
	t12 = t30 * t88 - t47 * t93;
	t11 = t30 * t47 + t88 * t93;
	t10 = t29 * t88 - t94 * t47;
	t9 = t29 * t47 + t94 * t88;
	t8 = t18 * t50 + t30 * t86;
	t6 = t16 * t50 + t29 * t86;
	t4 = t12 * t50 + t24 * t46;
	t2 = t10 * t50 + t23 * t46;
	t1 = [0, t41 * t81, (t17 * t45 + t8 * t49) * r_i_i_C(1) + (t17 * t49 - t8 * t45) * r_i_i_C(2) + t8 * pkin(5) + t18 * pkin(4) + t17 * pkin(11) - t55 * pkin(3) + t30 * t89 + t90 * (t18 * t46 - t30 * t85), t64 * t11 + t71 * t12, t90 * t4 + t72 * (-t12 * t46 + t24 * t50), (t11 * t49 - t4 * t45) * r_i_i_C(1) + (-t11 * t45 - t4 * t49) * r_i_i_C(2); 0, -t82 * t81, (t15 * t45 + t6 * t49) * r_i_i_C(1) + (t15 * t49 - t6 * t45) * r_i_i_C(2) + t6 * pkin(5) + t16 * pkin(4) + t15 * pkin(11) - t54 * pkin(3) + t29 * t89 + t90 * (t16 * t46 - t29 * t85), t71 * t10 + t64 * t9, t90 * t2 + t72 * (-t10 * t46 + t23 * t50), (-t2 * t45 + t9 * t49) * r_i_i_C(1) + (-t2 * t49 - t9 * t45) * r_i_i_C(2); 1, t44, (t20 * t49 + t25 * t45) * r_i_i_C(1) + (-t20 * t45 + t25 * t49) * r_i_i_C(2) + t20 * pkin(5) + t26 * pkin(4) + t25 * pkin(11) - t60 * pkin(3) + t35 * t89 + t90 * (t26 * t46 - t35 * t85), t64 * t21 + t71 * t22, t90 * t14 + t72 * (-t22 * t46 + t31 * t50), (-t14 * t45 + t21 * t49) * r_i_i_C(1) + (-t14 * t49 - t21 * t45) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end