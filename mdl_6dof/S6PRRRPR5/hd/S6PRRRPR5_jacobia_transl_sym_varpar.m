% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(12));
	t1 = sin(pkin(12));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (55->24), mult. (156->50), div. (0->0), fcn. (196->10), ass. (0->24)
	t7 = sin(pkin(7));
	t8 = sin(pkin(6));
	t24 = t7 * t8;
	t10 = cos(pkin(7));
	t15 = cos(qJ(2));
	t23 = t10 * t15;
	t11 = cos(pkin(6));
	t13 = sin(qJ(2));
	t22 = t11 * t13;
	t21 = t11 * t15;
	t12 = sin(qJ(3));
	t14 = cos(qJ(3));
	t20 = r_i_i_C(1) * t14 - r_i_i_C(2) * t12;
	t6 = sin(pkin(12));
	t9 = cos(pkin(12));
	t1 = -t6 * t13 + t9 * t21;
	t19 = -t1 * t10 + t9 * t24;
	t3 = -t9 * t13 - t6 * t21;
	t18 = t10 * t3 + t6 * t24;
	t17 = pkin(2) + t20;
	t16 = (pkin(9) + r_i_i_C(3)) * t7 + (-t12 * r_i_i_C(1) - t14 * r_i_i_C(2)) * t10;
	t4 = t15 * t9 - t6 * t22;
	t2 = t15 * t6 + t9 * t22;
	t5 = [0, t16 * t4 + t17 * t3, (-t12 * t4 + t18 * t14) * r_i_i_C(1) + (-t18 * t12 - t14 * t4) * r_i_i_C(2), 0, 0, 0; 0, t17 * t1 + t16 * t2, (-t12 * t2 - t19 * t14) * r_i_i_C(1) + (t19 * t12 - t14 * t2) * r_i_i_C(2), 0, 0, 0; 1, (t16 * t13 + t17 * t15) * t8, t20 * t7 * t11 + ((-t12 * t13 + t14 * t23) * r_i_i_C(1) + (-t12 * t23 - t13 * t14) * r_i_i_C(2)) * t8, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (168->57), mult. (471->110), div. (0->0), fcn. (606->12), ass. (0->43)
	t50 = r_i_i_C(3) + pkin(10);
	t22 = sin(pkin(7));
	t49 = t22 * pkin(9);
	t23 = sin(pkin(6));
	t48 = t22 * t23;
	t26 = cos(pkin(6));
	t47 = t22 * t26;
	t27 = sin(qJ(4));
	t46 = t22 * t27;
	t30 = cos(qJ(4));
	t45 = t22 * t30;
	t25 = cos(pkin(7));
	t44 = t23 * t25;
	t28 = sin(qJ(3));
	t43 = t25 * t28;
	t31 = cos(qJ(3));
	t42 = t25 * t31;
	t29 = sin(qJ(2));
	t41 = t26 * t29;
	t32 = cos(qJ(2));
	t40 = t26 * t32;
	t39 = t28 * t29;
	t38 = t28 * t32;
	t37 = t29 * t31;
	t36 = t31 * t32;
	t35 = t30 * r_i_i_C(1) - t27 * r_i_i_C(2) + pkin(3);
	t21 = sin(pkin(12));
	t24 = cos(pkin(12));
	t16 = -t21 * t29 + t24 * t40;
	t34 = t16 * t25 - t24 * t48;
	t18 = -t21 * t40 - t24 * t29;
	t33 = t18 * t25 + t21 * t48;
	t19 = -t21 * t41 + t24 * t32;
	t17 = t21 * t32 + t24 * t41;
	t15 = t26 * t25 - t32 * t48;
	t12 = -t18 * t22 + t21 * t44;
	t11 = -t16 * t22 - t24 * t44;
	t10 = t28 * t47 + (t25 * t38 + t37) * t23;
	t8 = t18 * t31 - t19 * t43;
	t6 = t16 * t31 - t17 * t43;
	t4 = t19 * t31 + t33 * t28;
	t2 = t17 * t31 + t34 * t28;
	t1 = [0, (t19 * t46 + t8 * t30) * r_i_i_C(1) + (t19 * t45 - t8 * t27) * r_i_i_C(2) + t8 * pkin(3) + t18 * pkin(2) + t19 * t49 + t50 * (t18 * t28 + t19 * t42), t50 * t4 + t35 * (-t19 * t28 + t33 * t31), (t12 * t30 - t4 * t27) * r_i_i_C(1) + (-t12 * t27 - t4 * t30) * r_i_i_C(2), 0, 0; 0, (t17 * t46 + t6 * t30) * r_i_i_C(1) + (t17 * t45 - t6 * t27) * r_i_i_C(2) + t6 * pkin(3) + t16 * pkin(2) + t17 * t49 + t50 * (t16 * t28 + t17 * t42), t50 * t2 + t35 * (-t17 * t28 + t34 * t31), (t11 * t30 - t2 * t27) * r_i_i_C(1) + (-t11 * t27 - t2 * t30) * r_i_i_C(2), 0, 0; 1, (t35 * (-t25 * t39 + t36) + t32 * pkin(2) + (t27 * r_i_i_C(1) + t30 * r_i_i_C(2) + pkin(9)) * t29 * t22 + t50 * (t25 * t37 + t38)) * t23, t50 * t10 + t35 * (t31 * t47 + (t25 * t36 - t39) * t23), (-t10 * t27 + t15 * t30) * r_i_i_C(1) + (-t10 * t30 - t15 * t27) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:08
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (246->60), mult. (582->109), div. (0->0), fcn. (747->14), ass. (0->46)
	t59 = r_i_i_C(3) + qJ(5) + pkin(10);
	t30 = sin(pkin(7));
	t31 = sin(pkin(6));
	t58 = t30 * t31;
	t33 = cos(pkin(7));
	t57 = t31 * t33;
	t36 = sin(qJ(3));
	t56 = t33 * t36;
	t39 = cos(qJ(3));
	t55 = t33 * t39;
	t37 = sin(qJ(2));
	t54 = t36 * t37;
	t40 = cos(qJ(2));
	t53 = t36 * t40;
	t52 = t37 * t39;
	t51 = t39 * t40;
	t50 = cos(pkin(6));
	t49 = sin(pkin(12));
	t32 = cos(pkin(12));
	t48 = t32 * t58;
	t47 = t31 * t49;
	t46 = t32 * t50;
	t45 = t50 * t30;
	t44 = t30 * t47;
	t43 = t50 * t49;
	t29 = qJ(4) + pkin(13);
	t27 = sin(t29);
	t28 = cos(t29);
	t38 = cos(qJ(4));
	t42 = t38 * pkin(4) + t28 * r_i_i_C(1) - t27 * r_i_i_C(2) + pkin(3);
	t35 = sin(qJ(4));
	t41 = (t35 * pkin(4) + t27 * r_i_i_C(1) + t28 * r_i_i_C(2) + pkin(9)) * t30;
	t21 = t32 * t40 - t37 * t43;
	t20 = -t32 * t37 - t40 * t43;
	t19 = t37 * t46 + t49 * t40;
	t18 = -t49 * t37 + t40 * t46;
	t17 = t50 * t33 - t40 * t58;
	t12 = -t20 * t30 + t33 * t47;
	t11 = -t18 * t30 - t32 * t57;
	t10 = t36 * t45 + (t33 * t53 + t52) * t31;
	t9 = t31 * t54 - t39 * t45 - t51 * t57;
	t4 = t21 * t39 + (t20 * t33 + t44) * t36;
	t3 = -t20 * t55 + t21 * t36 - t39 * t44;
	t2 = t19 * t39 + (t18 * t33 - t48) * t36;
	t1 = -t18 * t55 + t19 * t36 + t39 * t48;
	t5 = [0, t20 * pkin(2) + t59 * (t20 * t36 + t21 * t55) + t42 * (t20 * t39 - t21 * t56) + t21 * t41, -t42 * t3 + t59 * t4, (t12 * t28 - t4 * t27) * r_i_i_C(1) + (-t12 * t27 - t4 * t28) * r_i_i_C(2) + (t12 * t38 - t4 * t35) * pkin(4), t3, 0; 0, t18 * pkin(2) + t59 * (t18 * t36 + t19 * t55) + t42 * (t18 * t39 - t19 * t56) + t19 * t41, -t42 * t1 + t59 * t2, (t11 * t28 - t2 * t27) * r_i_i_C(1) + (-t11 * t27 - t2 * t28) * r_i_i_C(2) + (t11 * t38 - t2 * t35) * pkin(4), t1, 0; 1, (t59 * (t33 * t52 + t53) + t42 * (-t33 * t54 + t51) + t40 * pkin(2) + t37 * t41) * t31, t59 * t10 - t42 * t9, (-t10 * t27 + t17 * t28) * r_i_i_C(1) + (-t10 * t28 - t17 * t27) * r_i_i_C(2) + (-t10 * t35 + t17 * t38) * pkin(4), t9, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:08
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (538->94), mult. (1246->168), div. (0->0), fcn. (1622->16), ass. (0->64)
	t40 = sin(pkin(7));
	t47 = sin(qJ(2));
	t50 = cos(qJ(2));
	t42 = cos(pkin(12));
	t70 = cos(pkin(6));
	t64 = t42 * t70;
	t68 = sin(pkin(12));
	t54 = t68 * t47 - t50 * t64;
	t69 = cos(pkin(7));
	t41 = sin(pkin(6));
	t73 = t41 * t42;
	t79 = t40 * t73 + t54 * t69;
	t58 = t70 * t68;
	t55 = t42 * t47 + t50 * t58;
	t65 = t41 * t68;
	t78 = -t40 * t65 + t55 * t69;
	t77 = r_i_i_C(3) + pkin(11);
	t76 = cos(qJ(3));
	t39 = qJ(4) + pkin(13);
	t37 = sin(t39);
	t75 = t37 * t40;
	t38 = cos(t39);
	t74 = t38 * t40;
	t72 = t41 * t47;
	t71 = t41 * t50;
	t67 = t40 * t72;
	t46 = sin(qJ(3));
	t63 = t46 * t69;
	t62 = t70 * t40;
	t45 = sin(qJ(4));
	t60 = (pkin(4) * t45 + pkin(9)) * t40;
	t59 = t69 * t76;
	t44 = sin(qJ(6));
	t48 = cos(qJ(6));
	t57 = t48 * r_i_i_C(1) - t44 * r_i_i_C(2) + pkin(5);
	t43 = -qJ(5) - pkin(10);
	t56 = t44 * r_i_i_C(1) + t48 * r_i_i_C(2) - t43;
	t49 = cos(qJ(4));
	t36 = t49 * pkin(4) + pkin(3);
	t53 = -t77 * t37 - t57 * t38 - t36;
	t31 = t42 * t50 - t47 * t58;
	t30 = t47 * t64 + t68 * t50;
	t29 = -t40 * t71 + t70 * t69;
	t28 = (-t47 * t63 + t76 * t50) * t41;
	t27 = (t46 * t50 + t47 * t59) * t41;
	t24 = t55 * t40 + t69 * t65;
	t23 = t54 * t40 - t69 * t73;
	t22 = t46 * t62 + (t76 * t47 + t50 * t63) * t41;
	t21 = t46 * t72 - t59 * t71 - t76 * t62;
	t20 = t28 * t38 + t37 * t67;
	t18 = -t31 * t63 - t55 * t76;
	t17 = t31 * t59 - t55 * t46;
	t16 = -t30 * t63 - t54 * t76;
	t15 = t30 * t59 - t54 * t46;
	t14 = t31 * t76 - t78 * t46;
	t13 = t31 * t46 + t78 * t76;
	t12 = t30 * t76 - t79 * t46;
	t11 = t30 * t46 + t79 * t76;
	t10 = t22 * t38 + t29 * t37;
	t8 = t18 * t38 + t31 * t75;
	t6 = t16 * t38 + t30 * t75;
	t4 = t14 * t38 + t24 * t37;
	t2 = t12 * t38 + t23 * t37;
	t1 = [0, (t17 * t44 + t8 * t48) * r_i_i_C(1) + (t17 * t48 - t8 * t44) * r_i_i_C(2) + t8 * pkin(5) + t18 * t36 - t17 * t43 - t55 * pkin(2) + t77 * (t18 * t37 - t31 * t74) + t31 * t60, t53 * t13 + t56 * t14, t77 * t4 + (-t14 * t45 + t24 * t49) * pkin(4) + t57 * (-t14 * t37 + t24 * t38), t13, (t13 * t48 - t4 * t44) * r_i_i_C(1) + (-t13 * t44 - t4 * t48) * r_i_i_C(2); 0, (t15 * t44 + t6 * t48) * r_i_i_C(1) + (t15 * t48 - t6 * t44) * r_i_i_C(2) + t6 * pkin(5) + t16 * t36 - t15 * t43 - t54 * pkin(2) + t77 * (t16 * t37 - t30 * t74) + t30 * t60, t53 * t11 + t56 * t12, t77 * t2 + (-t12 * t45 + t23 * t49) * pkin(4) + t57 * (-t12 * t37 + t23 * t38), t11, (t11 * t48 - t2 * t44) * r_i_i_C(1) + (-t11 * t44 - t2 * t48) * r_i_i_C(2); 1, (t20 * t48 + t27 * t44) * r_i_i_C(1) + (-t20 * t44 + t27 * t48) * r_i_i_C(2) + t20 * pkin(5) + t28 * t36 - t27 * t43 + t77 * (t28 * t37 - t38 * t67) + (t50 * pkin(2) + t47 * t60) * t41, t53 * t21 + t56 * t22, t77 * t10 + (-t22 * t45 + t29 * t49) * pkin(4) + t57 * (-t22 * t37 + t29 * t38), t21, (-t10 * t44 + t21 * t48) * r_i_i_C(1) + (-t10 * t48 - t21 * t44) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end