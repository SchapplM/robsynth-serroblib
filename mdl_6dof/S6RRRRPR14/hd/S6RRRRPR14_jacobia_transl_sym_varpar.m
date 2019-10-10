% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:54
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR14_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR14_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR14_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (19->14), mult. (46->24), div. (0->0), fcn. (54->6), ass. (0->16)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t15 = t8 * t7;
	t9 = cos(qJ(2));
	t14 = t8 * t9;
	t10 = cos(qJ(1));
	t13 = t10 * t7;
	t12 = t10 * t9;
	t5 = sin(pkin(6));
	t11 = (pkin(9) + r_i_i_C(3)) * t5;
	t6 = cos(pkin(6));
	t4 = -t6 * t15 + t12;
	t3 = -t6 * t14 - t13;
	t2 = -t6 * t13 - t14;
	t1 = -t6 * t12 + t15;
	t16 = [-t8 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t10 * t11, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0, 0; t10 * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t11 * t8, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0; 0, (r_i_i_C(1) * t9 - r_i_i_C(2) * t7) * t5, 0, 0, 0, 0;];
	Ja_transl = t16;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (93->46), mult. (252->84), div. (0->0), fcn. (316->10), ass. (0->33)
	t39 = pkin(10) + r_i_i_C(3);
	t12 = cos(pkin(7));
	t13 = sin(qJ(3));
	t16 = cos(qJ(3));
	t10 = sin(pkin(7));
	t11 = sin(pkin(6));
	t18 = cos(qJ(1));
	t35 = t11 * t18;
	t27 = t10 * t35;
	t14 = sin(qJ(2));
	t15 = sin(qJ(1));
	t17 = cos(qJ(2));
	t28 = cos(pkin(6));
	t24 = t18 * t28;
	t3 = t15 * t14 - t17 * t24;
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
	t19 = t13 * t27 - t4 * t16 + t3 * t34;
	t6 = -t14 * t25 + t18 * t17;
	t2 = t22 * t13 + t6 * t16;
	t1 = -t6 * t13 + t22 * t16;
	t7 = [t19 * r_i_i_C(1) + t38 * r_i_i_C(2) - t4 * pkin(2) - t15 * pkin(1) + pkin(9) * t35 + t39 * (-t3 * t10 + t12 * t35), (t5 * t16 - t6 * t34) * r_i_i_C(1) + (-t5 * t13 - t6 * t33) * r_i_i_C(2) + t5 * pkin(2) + t6 * t26, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; t18 * pkin(1) + t6 * pkin(2) + pkin(9) * t36 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t39 * (-t5 * t10 + t12 * t36), (-t3 * t16 - t4 * t34) * r_i_i_C(1) + (t3 * t13 - t4 * t33) * r_i_i_C(2) - t3 * pkin(2) + t4 * t26, -t38 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, ((-t12 * t32 + t29) * r_i_i_C(1) + (-t12 * t30 - t31) * r_i_i_C(2) + t17 * pkin(2) + t14 * t26) * t11, ((t12 * t29 - t32) * r_i_i_C(1) + (-t12 * t31 - t30) * r_i_i_C(2)) * t11 + (t16 * r_i_i_C(1) - t13 * r_i_i_C(2)) * t10 * t28, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:13
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (242->73), mult. (663->129), div. (0->0), fcn. (854->12), ass. (0->49)
	t35 = sin(qJ(2));
	t36 = sin(qJ(1));
	t39 = cos(qJ(2));
	t40 = cos(qJ(1));
	t48 = cos(pkin(6));
	t44 = t40 * t48;
	t22 = t36 * t35 - t39 * t44;
	t30 = sin(pkin(7));
	t32 = cos(pkin(7));
	t31 = sin(pkin(6));
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
	t60 = r_i_i_C(3) + pkin(11);
	t59 = t30 * pkin(10);
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
	t3 = [-t36 * pkin(1) - t23 * pkin(2) + t6 * pkin(3) + pkin(9) * t55 + t15 * pkin(10) + t61 * r_i_i_C(1) + t62 * r_i_i_C(2) + t60 * t41, (t12 * t37 + t25 * t58) * r_i_i_C(1) + (-t12 * t33 + t25 * t57) * r_i_i_C(2) + t12 * pkin(3) + t24 * pkin(2) + t25 * t59 + t60 * (t24 * t34 + t25 * t53), -t42 * t7 + t60 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t40 * pkin(1) + t25 * pkin(2) + t8 * pkin(3) + pkin(9) * t56 + t17 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t60 * t7, (t10 * t37 + t23 * t58) * r_i_i_C(1) + (-t10 * t33 + t23 * t57) * r_i_i_C(2) + t10 * pkin(3) - t22 * pkin(2) + t23 * t59 + t60 * (-t22 * t34 + t23 * t53), t42 * t41 - t6 * t60, -t62 * r_i_i_C(1) + t61 * r_i_i_C(2), 0, 0; 0, (t42 * (-t32 * t52 + t49) + t39 * pkin(2) + (t33 * r_i_i_C(1) + t37 * r_i_i_C(2) + pkin(10)) * t35 * t30 + t60 * (t32 * t50 + t51)) * t31, t60 * t14 + t42 * (t38 * t43 + (t32 * t49 - t52) * t31), (-t14 * t33 + t21 * t37) * r_i_i_C(1) + (-t14 * t37 - t21 * t33) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:13
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (488->100), mult. (1351->170), div. (0->0), fcn. (1766->14), ass. (0->64)
	t48 = sin(qJ(2));
	t51 = cos(qJ(2));
	t52 = cos(qJ(1));
	t67 = cos(pkin(6));
	t63 = t52 * t67;
	t81 = sin(qJ(1));
	t34 = t81 * t48 - t51 * t63;
	t35 = t48 * t63 + t81 * t51;
	t47 = sin(qJ(3));
	t50 = cos(qJ(3));
	t42 = sin(pkin(7));
	t43 = sin(pkin(6));
	t76 = t43 * t52;
	t65 = t42 * t76;
	t45 = cos(pkin(7));
	t75 = t45 * t47;
	t16 = t34 * t75 - t35 * t50 + t47 * t65;
	t29 = -t34 * t42 + t45 * t76;
	t46 = sin(qJ(4));
	t49 = cos(qJ(4));
	t83 = t16 * t46 - t29 * t49;
	t4 = t16 * t49 + t29 * t46;
	t41 = sin(pkin(13));
	t44 = cos(pkin(13));
	t59 = r_i_i_C(1) * t44 - r_i_i_C(2) * t41 + pkin(4);
	t68 = r_i_i_C(3) + qJ(5);
	t54 = t68 * t46 + t59 * t49 + pkin(3);
	t82 = pkin(10) * t42;
	t78 = t42 * t46;
	t77 = t42 * t49;
	t74 = t45 * t50;
	t73 = t47 * t48;
	t72 = t47 * t51;
	t71 = t48 * t42;
	t70 = t48 * t50;
	t69 = t50 * t51;
	t66 = t43 * t71;
	t64 = t43 * t81;
	t62 = t67 * t42;
	t61 = t42 * t64;
	t60 = t67 * t81;
	t58 = t41 * r_i_i_C(1) + t44 * r_i_i_C(2) + pkin(11);
	t57 = -t43 * t51 * t42 + t67 * t45;
	t56 = t52 * t48 + t51 * t60;
	t55 = t56 * t50;
	t15 = -t35 * t47 + (-t34 * t45 - t65) * t50;
	t53 = t56 * t42 + t45 * t64;
	t36 = -t48 * t60 + t52 * t51;
	t32 = (-t45 * t73 + t69) * t43;
	t31 = (t45 * t70 + t72) * t43;
	t28 = t47 * t62 + (t45 * t72 + t70) * t43;
	t24 = t32 * t49 + t46 * t66;
	t22 = -t36 * t75 - t55;
	t21 = t36 * t74 - t56 * t47;
	t20 = -t34 * t50 - t35 * t75;
	t19 = -t34 * t47 + t35 * t74;
	t18 = t36 * t50 + (-t56 * t45 + t61) * t47;
	t17 = t36 * t47 + t45 * t55 - t50 * t61;
	t11 = t28 * t46 - t57 * t49;
	t10 = t22 * t49 + t36 * t78;
	t8 = t20 * t49 + t35 * t78;
	t6 = t18 * t49 + t53 * t46;
	t5 = t18 * t46 - t53 * t49;
	t1 = [(t15 * t41 + t4 * t44) * r_i_i_C(1) + (t15 * t44 - t4 * t41) * r_i_i_C(2) + t4 * pkin(4) + t16 * pkin(3) + t15 * pkin(11) - t35 * pkin(2) - t81 * pkin(1) + pkin(9) * t76 + t68 * t83 + t29 * pkin(10), (t10 * t44 + t21 * t41) * r_i_i_C(1) + (-t10 * t41 + t21 * t44) * r_i_i_C(2) + t10 * pkin(4) + t22 * pkin(3) + t21 * pkin(11) - t56 * pkin(2) + t36 * t82 + t68 * (t22 * t46 - t36 * t77), -t54 * t17 + t58 * t18, -t59 * t5 + t68 * t6, t5, 0; (t17 * t41 + t6 * t44) * r_i_i_C(1) + (t17 * t44 - t6 * t41) * r_i_i_C(2) + t6 * pkin(4) + t18 * pkin(3) + t17 * pkin(11) + t36 * pkin(2) + t52 * pkin(1) + pkin(9) * t64 + t68 * t5 + t53 * pkin(10), (t19 * t41 + t8 * t44) * r_i_i_C(1) + (t19 * t44 - t8 * t41) * r_i_i_C(2) + t8 * pkin(4) + t20 * pkin(3) + t19 * pkin(11) - t34 * pkin(2) + t35 * t82 + t68 * (t20 * t46 - t35 * t77), t54 * t15 - t16 * t58, -t68 * t4 + t59 * t83, -t83, 0; 0, (t24 * t44 + t31 * t41) * r_i_i_C(1) + (-t24 * t41 + t31 * t44) * r_i_i_C(2) + t24 * pkin(4) + t32 * pkin(3) + t31 * pkin(11) + (t51 * pkin(2) + pkin(10) * t71) * t43 + t68 * (t32 * t46 - t49 * t66), t58 * t28 + t54 * (t50 * t62 + (t45 * t69 - t73) * t43), t68 * (t28 * t49 + t57 * t46) - t59 * t11, t11, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:13
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (658->92), mult. (1640->161), div. (0->0), fcn. (2142->16), ass. (0->62)
	t55 = sin(qJ(2));
	t57 = cos(qJ(2));
	t58 = cos(qJ(1));
	t79 = cos(pkin(6));
	t71 = t58 * t79;
	t85 = sin(qJ(1));
	t37 = t85 * t55 - t57 * t71;
	t38 = t55 * t71 + t85 * t57;
	t54 = sin(qJ(3));
	t78 = cos(pkin(7));
	t72 = t54 * t78;
	t50 = sin(pkin(7));
	t51 = sin(pkin(6));
	t81 = t51 * t58;
	t76 = t50 * t81;
	t86 = cos(qJ(3));
	t16 = -t37 * t72 + t38 * t86 - t54 * t76;
	t53 = sin(qJ(4));
	t56 = cos(qJ(4));
	t73 = t51 * t78;
	t65 = t37 * t50 - t58 * t73;
	t4 = t16 * t56 + t65 * t53;
	t3 = t16 * t53 - t65 * t56;
	t68 = t79 * t85;
	t63 = t58 * t55 + t57 * t68;
	t74 = t51 * t85;
	t89 = -t50 * t74 + t63 * t78;
	t88 = pkin(10) * t50;
	t87 = r_i_i_C(3) + pkin(12) + qJ(5);
	t84 = t50 * t53;
	t83 = t50 * t56;
	t82 = t51 * t57;
	t80 = t55 * t50;
	t77 = t51 * t80;
	t75 = sin(pkin(13)) * pkin(5) + pkin(11);
	t70 = t79 * t50;
	t67 = t78 * t86;
	t45 = cos(pkin(13)) * pkin(5) + pkin(4);
	t48 = pkin(13) + qJ(6);
	t46 = sin(t48);
	t47 = cos(t48);
	t66 = t47 * r_i_i_C(1) - t46 * r_i_i_C(2) + t45;
	t64 = t46 * r_i_i_C(1) + t47 * r_i_i_C(2) + t75;
	t62 = -t50 * t82 + t79 * t78;
	t60 = -t87 * t53 - t66 * t56 - pkin(3);
	t15 = t37 * t67 + t38 * t54 + t86 * t76;
	t59 = t63 * t50 + t85 * t73;
	t39 = -t55 * t68 + t58 * t57;
	t35 = (-t55 * t72 + t86 * t57) * t51;
	t30 = t54 * t70 + (t86 * t55 + t57 * t72) * t51;
	t29 = t51 * t55 * t54 - t67 * t82 - t86 * t70;
	t24 = -t39 * t72 - t63 * t86;
	t22 = -t37 * t86 - t38 * t72;
	t20 = t39 * t86 - t89 * t54;
	t19 = t39 * t54 + t89 * t86;
	t14 = t30 * t56 + t62 * t53;
	t13 = t30 * t53 - t62 * t56;
	t8 = t20 * t56 + t59 * t53;
	t7 = t20 * t53 - t59 * t56;
	t2 = t19 * t46 + t8 * t47;
	t1 = t19 * t47 - t8 * t46;
	t5 = [-t85 * pkin(1) - t38 * pkin(2) - t16 * pkin(3) + pkin(9) * t81 - t65 * pkin(10) - t64 * t15 - t87 * t3 - t66 * t4, t24 * pkin(3) - t63 * pkin(2) + t39 * t88 + t66 * (t24 * t56 + t39 * t84) + t64 * (t39 * t67 - t63 * t54) + t87 * (t24 * t53 - t39 * t83), t60 * t19 + t64 * t20, -t66 * t7 + t87 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t58 * pkin(1) + t39 * pkin(2) + t20 * pkin(3) + pkin(9) * t74 + t59 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t75 * t19 + t8 * t45 + t87 * t7, t38 * t88 - t37 * pkin(2) + t22 * pkin(3) + t87 * (t22 * t53 - t38 * t83) + t66 * (t22 * t56 + t38 * t84) + t64 * (-t37 * t54 + t38 * t67), t60 * t15 + t64 * t16, -t66 * t3 + t87 * t4, t3, (t15 * t47 - t4 * t46) * r_i_i_C(1) + (-t15 * t46 - t4 * t47) * r_i_i_C(2); 0, t35 * pkin(3) + t66 * (t35 * t56 + t53 * t77) + t87 * (t35 * t53 - t56 * t77) + (t57 * pkin(2) + pkin(10) * t80 + t64 * (t54 * t57 + t55 * t67)) * t51, t60 * t29 + t64 * t30, -t66 * t13 + t87 * t14, t13, (-t14 * t46 + t29 * t47) * r_i_i_C(1) + (-t14 * t47 - t29 * t46) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end