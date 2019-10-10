% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR11_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR11_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR11_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
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
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (13->12), mult. (28->19), div. (0->0), fcn. (35->6), ass. (0->12)
	t1 = sin(pkin(12));
	t5 = sin(qJ(1));
	t11 = t5 * t1;
	t3 = cos(pkin(12));
	t10 = t5 * t3;
	t6 = cos(qJ(1));
	t9 = t6 * t1;
	t8 = t6 * t3;
	t2 = sin(pkin(6));
	t7 = t2 * (r_i_i_C(3) + qJ(2));
	t4 = cos(pkin(6));
	t12 = [(-t4 * t9 - t10) * r_i_i_C(1) + (-t4 * t8 + t11) * r_i_i_C(2) - t5 * pkin(1) + t6 * t7, t5 * t2, 0, 0, 0, 0; (-t4 * t11 + t8) * r_i_i_C(1) + (-t4 * t10 - t9) * r_i_i_C(2) + t6 * pkin(1) + t5 * t7, -t6 * t2, 0, 0, 0, 0; 0, t4, 0, 0, 0, 0;];
	Ja_transl = t12;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (62->32), mult. (166->58), div. (0->0), fcn. (213->10), ass. (0->29)
	t34 = r_i_i_C(3) + pkin(9);
	t13 = cos(pkin(7));
	t15 = sin(qJ(3));
	t17 = cos(qJ(3));
	t10 = sin(pkin(7));
	t11 = sin(pkin(6));
	t18 = cos(qJ(1));
	t26 = t18 * t11;
	t24 = t10 * t26;
	t12 = cos(pkin(12));
	t14 = cos(pkin(6));
	t29 = t14 * t18;
	t16 = sin(qJ(1));
	t9 = sin(pkin(12));
	t32 = t16 * t9;
	t3 = -t12 * t29 + t32;
	t27 = t16 * t12;
	t4 = t9 * t29 + t27;
	t33 = (t13 * t3 + t24) * t17 + t4 * t15;
	t30 = t13 * t15;
	t28 = t16 * t11;
	t25 = t11 * qJ(2);
	t5 = -t14 * t27 - t18 * t9;
	t22 = t10 * t28 + t13 * t5;
	t19 = t15 * t24 - t17 * t4 + t3 * t30;
	t6 = t12 * t18 - t14 * t32;
	t2 = t22 * t15 + t17 * t6;
	t1 = -t15 * t6 + t22 * t17;
	t7 = [t19 * r_i_i_C(1) + t33 * r_i_i_C(2) - t4 * pkin(2) - t16 * pkin(1) + t18 * t25 + t34 * (-t3 * t10 + t13 * t26), t28, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; t18 * pkin(1) + t6 * pkin(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t25 + t34 * (-t10 * t5 + t13 * t28), -t26, -t33 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, t14, (t17 * r_i_i_C(1) - t15 * r_i_i_C(2)) * t14 * t10 + ((t12 * t13 * t17 - t15 * t9) * r_i_i_C(1) + (-t12 * t30 - t17 * t9) * r_i_i_C(2)) * t11, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (176->48), mult. (480->84), div. (0->0), fcn. (625->12), ass. (0->42)
	t27 = cos(pkin(6));
	t25 = cos(pkin(12));
	t33 = cos(qJ(1));
	t41 = t33 * t25;
	t22 = sin(pkin(12));
	t30 = sin(qJ(1));
	t46 = t30 * t22;
	t16 = -t27 * t41 + t46;
	t23 = sin(pkin(7));
	t26 = cos(pkin(7));
	t24 = sin(pkin(6));
	t42 = t33 * t24;
	t11 = -t16 * t23 + t26 * t42;
	t28 = sin(qJ(4));
	t31 = cos(qJ(4));
	t43 = t33 * t22;
	t44 = t30 * t25;
	t17 = t27 * t43 + t44;
	t29 = sin(qJ(3));
	t32 = cos(qJ(3));
	t39 = t23 * t42;
	t47 = t26 * t29;
	t6 = t16 * t47 - t17 * t32 + t29 * t39;
	t52 = t11 * t31 - t6 * t28;
	t51 = t11 * t28 + t6 * t31;
	t36 = t27 * t44 + t43;
	t45 = t30 * t24;
	t50 = -t23 * t45 + t36 * t26;
	t49 = r_i_i_C(3) + pkin(10);
	t48 = t23 * t27;
	t40 = t24 * qJ(2);
	t37 = t31 * r_i_i_C(1) - t28 * r_i_i_C(2) + pkin(3);
	t34 = -t17 * t29 + (-t16 * t26 - t39) * t32;
	t13 = t23 * t36 + t26 * t45;
	t18 = -t27 * t46 + t41;
	t15 = -t24 * t25 * t23 + t27 * t26;
	t10 = t29 * t48 + (t22 * t32 + t25 * t47) * t24;
	t8 = t18 * t32 - t50 * t29;
	t7 = t18 * t29 + t50 * t32;
	t2 = t13 * t28 + t8 * t31;
	t1 = t13 * t31 - t8 * t28;
	t3 = [-t30 * pkin(1) - t17 * pkin(2) + t6 * pkin(3) + t11 * pkin(9) + t51 * r_i_i_C(1) + t52 * r_i_i_C(2) + t33 * t40 + t49 * t34, t45, -t37 * t7 + t49 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t33 * pkin(1) + t18 * pkin(2) + t8 * pkin(3) + t13 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t30 * t40 + t49 * t7, -t42, t37 * t34 - t49 * t6, -t52 * r_i_i_C(1) + t51 * r_i_i_C(2), 0, 0; 0, t27, t49 * t10 + t37 * (t32 * t48 + (t25 * t26 * t32 - t22 * t29) * t24), (-t10 * t28 + t15 * t31) * r_i_i_C(1) + (-t10 * t31 - t15 * t28) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (367->59), mult. (1010->99), div. (0->0), fcn. (1330->14), ass. (0->44)
	t36 = cos(qJ(1));
	t52 = cos(pkin(12));
	t53 = cos(pkin(6));
	t45 = t53 * t52;
	t51 = sin(pkin(12));
	t59 = sin(qJ(1));
	t22 = -t36 * t45 + t59 * t51;
	t29 = sin(pkin(7));
	t32 = cos(pkin(7));
	t30 = sin(pkin(6));
	t55 = t36 * t30;
	t63 = -t22 * t32 - t29 * t55;
	t44 = t53 * t51;
	t23 = t36 * t44 + t59 * t52;
	t34 = sin(qJ(3));
	t60 = cos(qJ(3));
	t10 = t23 * t60 + t63 * t34;
	t19 = -t22 * t29 + t32 * t55;
	t33 = sin(qJ(4));
	t35 = cos(qJ(4));
	t1 = t10 * t33 + t19 * t35;
	t62 = t10 * t35 - t19 * t33;
	t40 = t36 * t51 + t59 * t45;
	t49 = t59 * t30;
	t61 = -t29 * t49 + t40 * t32;
	t28 = sin(pkin(13));
	t31 = cos(pkin(13));
	t43 = r_i_i_C(1) * t31 - r_i_i_C(2) * t28 + pkin(4);
	t54 = r_i_i_C(3) + qJ(5);
	t39 = t54 * t33 + t43 * t35 + pkin(3);
	t48 = t29 * t53;
	t47 = t32 * t52;
	t42 = r_i_i_C(1) * t28 + r_i_i_C(2) * t31 + pkin(10);
	t41 = -t30 * t52 * t29 + t53 * t32;
	t11 = -t23 * t34 + t63 * t60;
	t37 = t40 * t29 + t32 * t49;
	t24 = t36 * t52 - t59 * t44;
	t18 = t34 * t48 + (t34 * t47 + t60 * t51) * t30;
	t14 = t24 * t60 - t61 * t34;
	t13 = t24 * t34 + t61 * t60;
	t7 = t18 * t33 - t41 * t35;
	t6 = t14 * t35 + t37 * t33;
	t5 = t14 * t33 - t37 * t35;
	t2 = [(t11 * t28 - t31 * t62) * r_i_i_C(1) + (t11 * t31 + t28 * t62) * r_i_i_C(2) - t62 * pkin(4) - t10 * pkin(3) + t11 * pkin(10) - t23 * pkin(2) - t59 * pkin(1) + qJ(2) * t55 - t54 * t1 + t19 * pkin(9), t49, -t39 * t13 + t42 * t14, -t43 * t5 + t54 * t6, t5, 0; (t13 * t28 + t31 * t6) * r_i_i_C(1) + (t13 * t31 - t28 * t6) * r_i_i_C(2) + t6 * pkin(4) + t14 * pkin(3) + t13 * pkin(10) + t24 * pkin(2) + t36 * pkin(1) + qJ(2) * t49 + t54 * t5 + t37 * pkin(9), -t55, t42 * t10 + t39 * t11, -t43 * t1 + t54 * t62, t1, 0; 0, t53, t42 * t18 + t39 * (t60 * t48 + (-t51 * t34 + t60 * t47) * t30), t54 * (t18 * t35 + t41 * t33) - t43 * t7, t7, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (509->66), mult. (1267->110), div. (0->0), fcn. (1667->16), ass. (0->53)
	t42 = cos(qJ(1));
	t69 = cos(pkin(12));
	t71 = cos(pkin(6));
	t59 = t71 * t69;
	t67 = sin(pkin(12));
	t73 = sin(qJ(1));
	t53 = -t42 * t59 + t73 * t67;
	t37 = sin(pkin(6));
	t68 = sin(pkin(7));
	t63 = t37 * t68;
	t70 = cos(pkin(7));
	t82 = t42 * t63 + t53 * t70;
	t57 = t71 * t67;
	t25 = t42 * t57 + t73 * t69;
	t40 = sin(qJ(3));
	t74 = cos(qJ(3));
	t14 = -t25 * t74 + t82 * t40;
	t39 = sin(qJ(4));
	t41 = cos(qJ(4));
	t64 = t37 * t70;
	t44 = -t42 * t64 + t53 * t68;
	t81 = t14 * t39 + t44 * t41;
	t80 = t14 * t41 - t44 * t39;
	t11 = t25 * t40 + t82 * t74;
	t48 = t42 * t67 + t73 * t59;
	t76 = t48 * t70 - t73 * t63;
	t75 = r_i_i_C(3) + pkin(11) + qJ(5);
	t72 = t42 * t37;
	t66 = sin(pkin(13)) * pkin(5) + pkin(10);
	t65 = t73 * t37;
	t58 = t71 * t68;
	t56 = t70 * t69;
	t32 = cos(pkin(13)) * pkin(5) + pkin(4);
	t35 = pkin(13) + qJ(6);
	t33 = sin(t35);
	t34 = cos(t35);
	t55 = t34 * r_i_i_C(1) - t33 * r_i_i_C(2) + t32;
	t52 = t33 * r_i_i_C(1) + t34 * r_i_i_C(2) + t66;
	t47 = -t75 * t39 - t55 * t41 - pkin(3);
	t46 = -t69 * t63 + t71 * t70;
	t43 = t48 * t68 + t73 * t64;
	t26 = t42 * t69 - t73 * t57;
	t20 = t40 * t58 + (t40 * t56 + t67 * t74) * t37;
	t19 = -t74 * t58 + (t40 * t67 - t56 * t74) * t37;
	t16 = t26 * t74 - t76 * t40;
	t15 = t26 * t40 + t76 * t74;
	t10 = t20 * t41 + t46 * t39;
	t9 = t20 * t39 - t46 * t41;
	t8 = t16 * t41 + t43 * t39;
	t7 = t16 * t39 - t43 * t41;
	t2 = t15 * t33 + t8 * t34;
	t1 = t15 * t34 - t8 * t33;
	t3 = [-t73 * pkin(1) - t25 * pkin(2) + t14 * pkin(3) - t44 * pkin(9) + qJ(2) * t72 - t52 * t11 + t55 * t80 + t75 * t81, t65, t47 * t15 + t52 * t16, -t55 * t7 + t75 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t42 * pkin(1) + t26 * pkin(2) + t16 * pkin(3) + t43 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t65 + t66 * t15 + t8 * t32 + t75 * t7, -t72, t47 * t11 - t14 * t52, t55 * t81 - t75 * t80, -t81, (t11 * t34 + t33 * t80) * r_i_i_C(1) + (-t11 * t33 + t34 * t80) * r_i_i_C(2); 0, t71, t47 * t19 + t52 * t20, t75 * t10 - t55 * t9, t9, (-t10 * t33 + t19 * t34) * r_i_i_C(1) + (-t10 * t34 - t19 * t33) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end