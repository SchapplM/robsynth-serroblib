% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:15
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPP1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPP1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPP1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
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
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(8) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t2 * t5 + t4 * t8, t6 * t4, 0, 0, 0, 0; t2 * t8 + t4 * t5, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (35->20), mult. (83->34), div. (0->0), fcn. (91->6), ass. (0->18)
	t16 = pkin(9) + r_i_i_C(3);
	t6 = sin(qJ(2));
	t14 = t16 * t6;
	t9 = cos(qJ(2));
	t18 = t9 * pkin(2) + pkin(1) + t14;
	t7 = sin(qJ(1));
	t17 = t7 * t9;
	t10 = cos(qJ(1));
	t15 = t10 * t9;
	t5 = sin(qJ(3));
	t8 = cos(qJ(3));
	t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5 + pkin(2);
	t11 = -t12 * t6 + t16 * t9;
	t4 = t8 * t15 + t7 * t5;
	t3 = -t5 * t15 + t7 * t8;
	t2 = t10 * t5 - t8 * t17;
	t1 = t10 * t8 + t5 * t17;
	t13 = [t10 * pkin(8) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t18 * t7, t11 * t10, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0; t7 * pkin(8) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t18 * t10, t11 * t7, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0, 0; 0, t12 * t9 + t14, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6, 0, 0, 0;];
	Ja_transl = t13;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (115->36), mult. (318->59), div. (0->0), fcn. (385->10), ass. (0->32)
	t10 = sin(qJ(3));
	t13 = cos(qJ(3));
	t6 = sin(pkin(10));
	t8 = cos(pkin(10));
	t22 = r_i_i_C(1) * t6 + r_i_i_C(2) * t8;
	t24 = r_i_i_C(3) + qJ(4);
	t7 = sin(pkin(6));
	t9 = cos(pkin(6));
	t18 = t22 * t9 - t24 * t7;
	t21 = t8 * r_i_i_C(1) - t6 * r_i_i_C(2) + pkin(3);
	t34 = t18 * t10 - t21 * t13 - pkin(2);
	t11 = sin(qJ(2));
	t20 = t22 * t7 + pkin(9);
	t17 = t24 * t9 + t20;
	t33 = t17 * t11;
	t12 = sin(qJ(1));
	t14 = cos(qJ(2));
	t15 = cos(qJ(1));
	t26 = t15 * t10;
	t4 = t12 * t13 - t14 * t26;
	t30 = t4 * t7;
	t29 = t4 * t9;
	t28 = t11 * t9;
	t27 = t12 * t14;
	t25 = t15 * t13;
	t23 = -t14 * pkin(2) - pkin(1);
	t16 = t34 * t11 + t17 * t14;
	t5 = t12 * t10 + t14 * t25;
	t3 = -t13 * t27 + t26;
	t2 = t10 * t27 + t25;
	t1 = t15 * t28 - t30;
	t19 = [t15 * pkin(8) + t21 * t3 + t18 * t2 + (t23 - t33) * t12, t16 * t15, -t18 * t5 + t21 * t4, t1, 0, 0; (t6 * t29 + t5 * t8) * r_i_i_C(1) + (t8 * t29 - t5 * t6) * r_i_i_C(2) + t1 * r_i_i_C(3) + t5 * pkin(3) - qJ(4) * t30 + t12 * pkin(8) + ((t9 * qJ(4) + t20) * t11 - t23) * t15, t16 * t12, t18 * t3 - t21 * t2, t12 * t28 + t2 * t7, 0, 0; 0, -t34 * t14 + t33, (-t21 * t10 - t13 * t18) * t11, t11 * t10 * t7 - t14 * t9, 0, 0;];
	Ja_transl = t19;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (183->62), mult. (513->101), div. (0->0), fcn. (635->10), ass. (0->45)
	t27 = sin(qJ(3));
	t30 = cos(qJ(3));
	t24 = sin(pkin(6));
	t47 = r_i_i_C(1) + qJ(4);
	t65 = t47 * t24;
	t67 = t30 * pkin(3) + t27 * t65 + pkin(2);
	t28 = sin(qJ(2));
	t26 = cos(pkin(6));
	t41 = t47 * t26 + pkin(9);
	t66 = t41 * t28;
	t23 = sin(pkin(10));
	t25 = cos(pkin(10));
	t31 = cos(qJ(2));
	t56 = t27 * t28;
	t37 = t24 * t31 + t26 * t56;
	t54 = t28 * t30;
	t36 = t23 * t54 + t37 * t25;
	t46 = r_i_i_C(3) + qJ(5);
	t62 = r_i_i_C(2) - pkin(4);
	t64 = t41 * t31 - t67 * t28 - t46 * t36 + t62 * (-t37 * t23 + t25 * t54);
	t29 = sin(qJ(1));
	t32 = cos(qJ(1));
	t49 = t32 * t27;
	t51 = t30 * t31;
	t19 = t29 * t51 - t49;
	t48 = t32 * t30;
	t52 = t29 * t27;
	t18 = t31 * t52 + t48;
	t55 = t28 * t29;
	t40 = -t18 * t26 + t24 * t55;
	t63 = -t19 * t23 + t40 * t25;
	t20 = t29 * t30 - t31 * t49;
	t60 = t20 * t24;
	t59 = t23 * t26;
	t58 = t25 * t26;
	t57 = t26 * t30;
	t53 = t28 * t32;
	t50 = t31 * t26;
	t45 = -t31 * pkin(2) - pkin(1);
	t39 = t20 * t26 + t24 * t53;
	t38 = -t24 * t28 + t27 * t50;
	t21 = t31 * t48 + t52;
	t15 = t26 * t53 - t60;
	t3 = t21 * t23 - t39 * t25;
	t1 = [-t19 * pkin(3) + t32 * pkin(8) - t62 * (-t19 * t25 - t40 * t23) - t18 * t65 + t46 * t63 + (t45 - t66) * t29, t64 * t32, t20 * pkin(3) - t62 * (t20 * t25 - t21 * t59) + t46 * (t20 * t23 + t21 * t58) + t21 * t65, t15, t3, 0; -qJ(4) * t60 + t21 * pkin(3) + t29 * pkin(8) + t15 * r_i_i_C(1) - t62 * (t21 * t25 + t39 * t23) + t46 * t3 + ((t26 * qJ(4) + pkin(9)) * t28 - t45) * t32, t64 * t29, -t18 * pkin(3) - t62 * (-t18 * t25 - t19 * t59) + t46 * (-t18 * t23 + t19 * t58) + t19 * t65, t18 * t24 + t26 * t55, -t63, 0; 0, -t62 * (-t38 * t23 + t25 * t51) + t46 * (t23 * t51 + t38 * t25) + t66 + t67 * t31, (t62 * (t23 * t57 + t25 * t27) - t46 * (t23 * t27 - t25 * t57) - pkin(3) * t27 + t30 * t65) * t28, t24 * t56 - t50, t36, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:15:12
	% EndTime: 2019-10-10 11:15:12
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (239->67), mult. (669->105), div. (0->0), fcn. (834->10), ass. (0->52)
	t32 = sin(pkin(10));
	t34 = cos(pkin(10));
	t36 = sin(qJ(3));
	t35 = cos(pkin(6));
	t37 = sin(qJ(2));
	t74 = t37 * t35;
	t33 = sin(pkin(6));
	t40 = cos(qJ(2));
	t77 = t33 * t40;
	t49 = t36 * t74 + t77;
	t39 = cos(qJ(3));
	t72 = t37 * t39;
	t86 = -t32 * t49 + t34 * t72;
	t57 = qJ(4) * t35 + pkin(9);
	t54 = t57 * t37;
	t84 = t40 * pkin(2) + pkin(1) + t54;
	t41 = cos(qJ(1));
	t66 = t41 * t39;
	t38 = sin(qJ(1));
	t70 = t38 * t36;
	t22 = t40 * t70 + t66;
	t67 = t41 * t36;
	t69 = t39 * t40;
	t23 = t38 * t69 - t67;
	t73 = t37 * t38;
	t62 = t33 * t73;
	t83 = (-t22 * t35 + t62) * t34 - t23 * t32;
	t44 = t32 * t72 + t49 * t34;
	t64 = qJ(4) * t33;
	t47 = pkin(3) * t39 + t36 * t64 + pkin(2);
	t68 = t40 * t35;
	t78 = t33 * t37;
	t48 = t36 * t78 - t68;
	t63 = pkin(4) + r_i_i_C(3) + qJ(6);
	t65 = r_i_i_C(2) + qJ(5);
	t81 = pkin(5) + r_i_i_C(1);
	t82 = -t47 * t37 + t57 * t40 - t65 * t44 - t81 * t48 - t63 * t86;
	t79 = t32 * t35;
	t76 = t34 * t35;
	t75 = t35 * t39;
	t71 = t37 * t41;
	t55 = (qJ(4) + t81) * t33;
	t53 = t22 * t33 + t35 * t73;
	t24 = t38 * t39 - t40 * t67;
	t51 = t24 * t35 + t33 * t71;
	t50 = t36 * t68 - t78;
	t46 = t22 * t79 - t23 * t34 - t32 * t62;
	t25 = t40 * t66 + t70;
	t16 = -t24 * t33 + t35 * t71;
	t4 = t25 * t34 + t51 * t32;
	t3 = t25 * t32 - t51 * t34;
	t1 = [-t23 * pkin(3) + t41 * pkin(8) - t22 * t64 - t84 * t38 + t63 * t46 - t81 * t53 + t65 * t83, t82 * t41, t24 * pkin(3) + t65 * (t24 * t32 + t25 * t76) + t63 * (t24 * t34 - t25 * t79) + t25 * t55, t16, t3, t4; t25 * pkin(3) + t38 * pkin(8) + t81 * t16 - t24 * t64 + t65 * t3 + t63 * t4 + t84 * t41, t82 * t38, -t22 * pkin(3) + t65 * (-t22 * t32 + t23 * t76) + t63 * (-t22 * t34 - t23 * t79) + t23 * t55, t53, -t83, -t46; 0, t54 + t81 * (t36 * t77 + t74) + t65 * (t32 * t69 + t50 * t34) + t47 * t40 + t63 * (-t32 * t50 + t34 * t69), (-t65 * (t32 * t36 - t34 * t75) - t63 * (t32 * t75 + t34 * t36) - pkin(3) * t36 + t39 * t55) * t37, t48, t44, t86;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end