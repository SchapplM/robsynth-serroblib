% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR12_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR12_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR12_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:58:31
	% EndTime: 2019-12-31 22:58:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:58:31
	% EndTime: 2019-12-31 22:58:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:58:31
	% EndTime: 2019-12-31 22:58:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (19->14), mult. (46->24), div. (0->0), fcn. (54->6), ass. (0->16)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t15 = t8 * t7;
	t9 = cos(qJ(2));
	t14 = t8 * t9;
	t10 = cos(qJ(1));
	t13 = t10 * t7;
	t12 = t10 * t9;
	t5 = sin(pkin(5));
	t11 = (pkin(8) + r_i_i_C(3)) * t5;
	t6 = cos(pkin(5));
	t4 = -t6 * t15 + t12;
	t3 = -t6 * t14 - t13;
	t2 = -t6 * t13 - t14;
	t1 = -t6 * t12 + t15;
	t16 = [-t8 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t10 * t11, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0; t10 * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t8 * t11, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0; 0, (r_i_i_C(1) * t9 - r_i_i_C(2) * t7) * t5, 0, 0, 0;];
	Ja_transl = t16;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:58:31
	% EndTime: 2019-12-31 22:58:32
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (93->46), mult. (252->84), div. (0->0), fcn. (316->10), ass. (0->33)
	t39 = pkin(9) + r_i_i_C(3);
	t12 = cos(pkin(6));
	t13 = sin(qJ(3));
	t16 = cos(qJ(3));
	t10 = sin(pkin(6));
	t11 = sin(pkin(5));
	t18 = cos(qJ(1));
	t35 = t11 * t18;
	t27 = t10 * t35;
	t14 = sin(qJ(2));
	t15 = sin(qJ(1));
	t17 = cos(qJ(2));
	t28 = cos(pkin(5));
	t24 = t18 * t28;
	t3 = t14 * t15 - t17 * t24;
	t4 = t14 * t24 + t15 * t17;
	t38 = (t12 * t3 + t27) * t16 + t4 * t13;
	t36 = t11 * t15;
	t34 = t12 * t13;
	t33 = t12 * t16;
	t32 = t13 * t14;
	t31 = t13 * t17;
	t30 = t14 * t16;
	t29 = t16 * t17;
	t26 = t10 * t39;
	t25 = t15 * t28;
	t5 = -t18 * t14 - t17 * t25;
	t22 = t10 * t36 + t12 * t5;
	t19 = t13 * t27 - t16 * t4 + t3 * t34;
	t6 = -t14 * t25 + t17 * t18;
	t2 = t22 * t13 + t16 * t6;
	t1 = -t13 * t6 + t22 * t16;
	t7 = [t19 * r_i_i_C(1) + t38 * r_i_i_C(2) - t4 * pkin(2) - t15 * pkin(1) + pkin(8) * t35 + t39 * (-t3 * t10 + t12 * t35), (t16 * t5 - t6 * t34) * r_i_i_C(1) + (-t13 * t5 - t6 * t33) * r_i_i_C(2) + t5 * pkin(2) + t6 * t26, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0; t18 * pkin(1) + t6 * pkin(2) + pkin(8) * t36 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t39 * (-t10 * t5 + t12 * t36), (-t16 * t3 - t4 * t34) * r_i_i_C(1) + (t13 * t3 - t4 * t33) * r_i_i_C(2) - t3 * pkin(2) + t4 * t26, -t38 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0; 0, ((-t12 * t32 + t29) * r_i_i_C(1) + (-t12 * t30 - t31) * r_i_i_C(2) + t17 * pkin(2) + t14 * t26) * t11, ((t12 * t29 - t32) * r_i_i_C(1) + (-t12 * t31 - t30) * r_i_i_C(2)) * t11 + (t16 * r_i_i_C(1) - t13 * r_i_i_C(2)) * t10 * t28, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:58:31
	% EndTime: 2019-12-31 22:58:32
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (242->73), mult. (663->129), div. (0->0), fcn. (854->12), ass. (0->49)
	t35 = sin(qJ(2));
	t36 = sin(qJ(1));
	t39 = cos(qJ(2));
	t40 = cos(qJ(1));
	t48 = cos(pkin(5));
	t44 = t40 * t48;
	t22 = t36 * t35 - t39 * t44;
	t30 = sin(pkin(6));
	t32 = cos(pkin(6));
	t31 = sin(pkin(5));
	t55 = t31 * t40;
	t15 = -t22 * t30 + t32 * t55;
	t33 = sin(qJ(4));
	t37 = cos(qJ(4));
	t23 = t35 * t44 + t36 * t39;
	t34 = sin(qJ(3));
	t38 = cos(qJ(3));
	t46 = t30 * t55;
	t54 = t32 * t34;
	t6 = t22 * t54 - t23 * t38 + t34 * t46;
	t62 = t15 * t37 - t6 * t33;
	t61 = t15 * t33 + t6 * t37;
	t60 = r_i_i_C(3) + pkin(10);
	t59 = t30 * pkin(9);
	t58 = t30 * t33;
	t57 = t30 * t37;
	t56 = t31 * t36;
	t53 = t32 * t38;
	t52 = t34 * t35;
	t51 = t34 * t39;
	t50 = t35 * t38;
	t49 = t38 * t39;
	t47 = t30 * t56;
	t45 = t36 * t48;
	t43 = t48 * t30;
	t42 = t37 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(3);
	t24 = -t40 * t35 - t39 * t45;
	t17 = -t24 * t30 + t32 * t56;
	t41 = -t23 * t34 + (-t22 * t32 - t46) * t38;
	t25 = -t35 * t45 + t40 * t39;
	t21 = -t31 * t39 * t30 + t48 * t32;
	t14 = t34 * t43 + (t32 * t51 + t50) * t31;
	t12 = t24 * t38 - t25 * t54;
	t10 = -t22 * t38 - t23 * t54;
	t8 = t25 * t38 + (t24 * t32 + t47) * t34;
	t7 = -t24 * t53 + t25 * t34 - t38 * t47;
	t2 = t17 * t33 + t8 * t37;
	t1 = t17 * t37 - t8 * t33;
	t3 = [-t36 * pkin(1) - t23 * pkin(2) + t6 * pkin(3) + pkin(8) * t55 + t15 * pkin(9) + t61 * r_i_i_C(1) + t62 * r_i_i_C(2) + t60 * t41, (t12 * t37 + t25 * t58) * r_i_i_C(1) + (-t12 * t33 + t25 * t57) * r_i_i_C(2) + t12 * pkin(3) + t24 * pkin(2) + t25 * t59 + t60 * (t24 * t34 + t25 * t53), -t42 * t7 + t60 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t40 * pkin(1) + t25 * pkin(2) + t8 * pkin(3) + pkin(8) * t56 + t17 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t60 * t7, (t10 * t37 + t23 * t58) * r_i_i_C(1) + (-t10 * t33 + t23 * t57) * r_i_i_C(2) + t10 * pkin(3) - t22 * pkin(2) + t23 * t59 + t60 * (-t22 * t34 + t23 * t53), t42 * t41 - t6 * t60, -t62 * r_i_i_C(1) + t61 * r_i_i_C(2), 0; 0, (t42 * (-t32 * t52 + t49) + t39 * pkin(2) + (t33 * r_i_i_C(1) + t37 * r_i_i_C(2) + pkin(9)) * t35 * t30 + t60 * (t32 * t50 + t51)) * t31, t60 * t14 + t42 * (t38 * t43 + (t32 * t49 - t52) * t31), (-t14 * t33 + t21 * t37) * r_i_i_C(1) + (-t14 * t37 - t21 * t33) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:58:31
	% EndTime: 2019-12-31 22:58:32
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (535->104), mult. (1483->180), div. (0->0), fcn. (1939->14), ass. (0->69)
	t51 = sin(qJ(2));
	t54 = cos(qJ(2));
	t55 = cos(qJ(1));
	t72 = cos(pkin(5));
	t65 = t55 * t72;
	t80 = sin(qJ(1));
	t37 = t80 * t51 - t54 * t65;
	t38 = t51 * t65 + t80 * t54;
	t50 = sin(qJ(3));
	t47 = cos(pkin(6));
	t81 = cos(qJ(3));
	t68 = t47 * t81;
	t45 = sin(pkin(6));
	t46 = sin(pkin(5));
	t77 = t46 * t55;
	t70 = t45 * t77;
	t15 = t37 * t68 + t38 * t50 + t81 * t70;
	t76 = t47 * t50;
	t16 = -t37 * t76 + t38 * t81 - t50 * t70;
	t30 = -t37 * t45 + t47 * t77;
	t49 = sin(qJ(4));
	t53 = cos(qJ(4));
	t4 = t16 * t53 - t30 * t49;
	t48 = sin(qJ(5));
	t52 = cos(qJ(5));
	t88 = -t15 * t52 + t4 * t48;
	t87 = -t15 * t48 - t4 * t52;
	t84 = -t16 * t49 - t30 * t53;
	t83 = r_i_i_C(3) + pkin(11);
	t82 = pkin(9) * t45;
	t79 = t45 * t49;
	t78 = t45 * t53;
	t75 = t50 * t51;
	t74 = t50 * t54;
	t73 = t51 * t45;
	t71 = t46 * t73;
	t69 = t46 * t80;
	t67 = t81 * t51;
	t66 = t81 * t54;
	t64 = t72 * t45;
	t63 = t45 * t69;
	t62 = t72 * t80;
	t61 = t52 * r_i_i_C(1) - t48 * r_i_i_C(2) + pkin(4);
	t60 = t48 * r_i_i_C(1) + t52 * r_i_i_C(2) + pkin(10);
	t59 = t55 * t51 + t54 * t62;
	t58 = t59 * t81;
	t57 = -t83 * t49 - t61 * t53 - pkin(3);
	t56 = t59 * t45 + t47 * t69;
	t39 = -t51 * t62 + t55 * t54;
	t36 = -t46 * t54 * t45 + t72 * t47;
	t35 = (-t47 * t75 + t66) * t46;
	t34 = (t47 * t67 + t74) * t46;
	t29 = t50 * t64 + (t47 * t74 + t67) * t46;
	t28 = -t81 * t64 + (-t47 * t66 + t75) * t46;
	t26 = t35 * t53 + t49 * t71;
	t24 = -t39 * t76 - t58;
	t23 = t39 * t68 - t59 * t50;
	t22 = -t37 * t81 - t38 * t76;
	t21 = -t37 * t50 + t38 * t68;
	t20 = t39 * t81 + (-t59 * t47 + t63) * t50;
	t19 = t39 * t50 + t47 * t58 - t81 * t63;
	t14 = t29 * t53 + t36 * t49;
	t12 = t24 * t53 + t39 * t79;
	t10 = t22 * t53 + t38 * t79;
	t8 = t20 * t53 + t56 * t49;
	t7 = t20 * t49 - t56 * t53;
	t2 = t19 * t48 + t8 * t52;
	t1 = t19 * t52 - t8 * t48;
	t3 = [-t80 * pkin(1) - t38 * pkin(2) - t16 * pkin(3) - t4 * pkin(4) + pkin(8) * t77 + t30 * pkin(9) - t15 * pkin(10) + t87 * r_i_i_C(1) + t88 * r_i_i_C(2) + t83 * t84, (t12 * t52 + t23 * t48) * r_i_i_C(1) + (-t12 * t48 + t23 * t52) * r_i_i_C(2) + t12 * pkin(4) + t24 * pkin(3) + t23 * pkin(10) - t59 * pkin(2) + t39 * t82 + t83 * (t24 * t49 - t39 * t78), t57 * t19 + t60 * t20, -t61 * t7 + t83 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t55 * pkin(1) + t39 * pkin(2) + t20 * pkin(3) + t8 * pkin(4) + pkin(8) * t69 + t56 * pkin(9) + t19 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t83 * t7, (t10 * t52 + t21 * t48) * r_i_i_C(1) + (-t10 * t48 + t21 * t52) * r_i_i_C(2) + t10 * pkin(4) + t22 * pkin(3) + t21 * pkin(10) - t37 * pkin(2) + t38 * t82 + t83 * (t22 * t49 - t38 * t78), t57 * t15 + t60 * t16, t83 * t4 + t61 * t84, -t88 * r_i_i_C(1) + t87 * r_i_i_C(2); 0, (t26 * t52 + t34 * t48) * r_i_i_C(1) + (-t26 * t48 + t34 * t52) * r_i_i_C(2) + t26 * pkin(4) + t35 * pkin(3) + t34 * pkin(10) + (t54 * pkin(2) + pkin(9) * t73) * t46 + t83 * (t35 * t49 - t53 * t71), t57 * t28 + t60 * t29, t83 * t14 + t61 * (-t29 * t49 + t36 * t53), (-t14 * t48 + t28 * t52) * r_i_i_C(1) + (-t14 * t52 - t28 * t48) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,5);
end