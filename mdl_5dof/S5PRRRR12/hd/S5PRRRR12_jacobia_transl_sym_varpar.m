% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR12
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRRR12_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR12_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR12_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(5));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(11));
	t1 = sin(pkin(11));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(5)), 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (87->19), mult. (64->29), div. (0->0), fcn. (60->12), ass. (0->20)
	t17 = qJ(2) + qJ(3);
	t15 = sin(t17);
	t16 = cos(t17);
	t18 = sin(pkin(11));
	t19 = cos(pkin(11));
	t13 = pkin(5) + t17;
	t11 = sin(t13);
	t14 = pkin(5) - t17;
	t9 = sin(t14) / 0.2e1;
	t7 = t9 - t11 / 0.2e1;
	t10 = cos(t13) / 0.2e1;
	t12 = cos(t14);
	t8 = t12 / 0.2e1 + t10;
	t26 = (-t18 * t15 + t19 * t8) * r_i_i_C(1) + (-t18 * t16 + t19 * t7) * r_i_i_C(2);
	t25 = (-t19 * t15 - t18 * t8) * r_i_i_C(1) + (-t19 * t16 - t18 * t7) * r_i_i_C(2);
	t24 = (t11 / 0.2e1 + t9) * r_i_i_C(1) + (t10 - t12 / 0.2e1) * r_i_i_C(2);
	t22 = cos(qJ(2));
	t23 = cos(pkin(5)) * t22;
	t21 = sin(qJ(2));
	t1 = [0, (-t18 * t23 - t19 * t21) * pkin(2) + t25, t25, 0, 0; 0, (-t18 * t21 + t19 * t23) * pkin(2) + t26, t26, 0, 0; 1, sin(pkin(5)) * t22 * pkin(2) + t24, t24, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (196->29), mult. (128->40), div. (0->0), fcn. (122->15), ass. (0->30)
	t33 = qJ(2) + qJ(3);
	t22 = qJ(4) + t33;
	t16 = pkin(5) + t22;
	t13 = cos(t16) / 0.2e1;
	t17 = pkin(5) - t22;
	t15 = cos(t17);
	t10 = t15 / 0.2e1 + t13;
	t18 = sin(t22);
	t19 = cos(t22);
	t23 = sin(pkin(11));
	t25 = cos(pkin(11));
	t12 = sin(t17) / 0.2e1;
	t14 = sin(t16);
	t9 = t12 - t14 / 0.2e1;
	t38 = (t10 * t25 - t18 * t23) * r_i_i_C(1) + (-t19 * t23 + t25 * t9) * r_i_i_C(2);
	t37 = (-t10 * t23 - t18 * t25) * r_i_i_C(1) + (-t19 * t25 - t23 * t9) * r_i_i_C(2);
	t36 = (t14 / 0.2e1 + t12) * r_i_i_C(1) + (t13 - t15 / 0.2e1) * r_i_i_C(2);
	t35 = pkin(3) * sin(t33);
	t28 = sin(qJ(2));
	t34 = sin(qJ(3)) * t28;
	t29 = cos(qJ(3));
	t30 = cos(qJ(2));
	t32 = -pkin(3) * t34 + (pkin(3) * t29 + pkin(2)) * t30;
	t31 = pkin(3) * (t29 * t30 - t34);
	t26 = cos(pkin(5));
	t24 = sin(pkin(5));
	t11 = -pkin(2) * t28 - t35;
	t6 = t26 * t31;
	t5 = t32 * t26;
	t1 = [0, t11 * t25 - t23 * t5 + t37, -t23 * t6 - t25 * t35 + t37, t37, 0; 0, t11 * t23 + t25 * t5 + t38, -t23 * t35 + t25 * t6 + t38, t38, 0; 1, t32 * t24 + t36, t24 * t31 + t36, t36, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (360->87), mult. (651->153), div. (0->0), fcn. (719->16), ass. (0->68)
	t54 = sin(qJ(3));
	t58 = cos(qJ(3));
	t46 = sin(pkin(11));
	t51 = cos(pkin(5));
	t82 = t46 * t51;
	t49 = cos(pkin(11));
	t47 = sin(pkin(6));
	t91 = pkin(10) * t47;
	t72 = t46 * t91;
	t32 = pkin(4) * t49 + t51 * t72;
	t71 = t49 * t91;
	t34 = pkin(4) * t82 - t71;
	t53 = sin(qJ(4));
	t57 = cos(qJ(4));
	t92 = t32 * t53 + t34 * t57;
	t7 = -pkin(3) * t82 - t92;
	t13 = -t32 * t57 + t34 * t53;
	t9 = pkin(3) * t49 - t13;
	t96 = -t54 * t9 + t7 * t58;
	t33 = -pkin(4) * t46 + t51 * t71;
	t79 = t49 * t51;
	t35 = pkin(4) * t79 + t72;
	t14 = t33 * t57 - t35 * t53;
	t12 = pkin(3) * t46 - t14;
	t64 = t33 * t53 + t35 * t57;
	t8 = pkin(3) * t79 + t64;
	t95 = -t12 * t54 + t8 * t58;
	t94 = -t7 * t54 - t58 * t9;
	t93 = -t12 * t58 - t8 * t54;
	t90 = r_i_i_C(3) * t47;
	t40 = -pkin(4) * t53 + t57 * t91;
	t83 = t40 * t54;
	t48 = sin(pkin(5));
	t81 = t47 * t48;
	t59 = cos(qJ(2));
	t80 = t48 * t59;
	t50 = cos(pkin(6));
	t52 = sin(qJ(5));
	t78 = t50 * t52;
	t56 = cos(qJ(5));
	t77 = t50 * t56;
	t76 = t51 * t52;
	t75 = t51 * t56;
	t68 = t50 * t76;
	t20 = t46 * t68 - t49 * t56;
	t67 = t50 * t75;
	t21 = -t46 * t67 - t49 * t52;
	t26 = t46 * t75 + t49 * t78;
	t27 = t46 * t76 - t49 * t77;
	t45 = qJ(2) + qJ(3) + qJ(4);
	t42 = sin(t45);
	t43 = cos(t45);
	t74 = (t20 * t42 - t26 * t43) * r_i_i_C(1) + (-t42 * t82 + t43 * t49) * t90 + (-t21 * t42 + t27 * t43) * r_i_i_C(2);
	t22 = t46 * t56 + t49 * t68;
	t23 = -t46 * t52 + t49 * t67;
	t24 = t46 * t78 - t49 * t75;
	t25 = -t46 * t77 - t49 * t76;
	t73 = (t42 * t79 + t43 * t46) * t90 + (-t22 * t42 - t24 * t43) * r_i_i_C(1) + (-t23 * t42 + t25 * t43) * r_i_i_C(2);
	t70 = t52 * t81;
	t69 = t56 * t81;
	t66 = t42 * r_i_i_C(3) * t81 + ((-t42 * t78 + t43 * t56) * r_i_i_C(1) + (-t42 * t77 - t43 * t52) * r_i_i_C(2)) * t48;
	t36 = t40 * t58;
	t39 = -pkin(4) * t57 - t53 * t91;
	t38 = pkin(3) - t39;
	t55 = sin(qJ(2));
	t63 = (-t38 * t54 + t36) * t48 * t55 + t66;
	t60 = -t38 * t58 - t83;
	t1 = [0, -(pkin(2) * t49 - t94) * t55 + (-pkin(2) * t82 + t96) * t59 + t74, t94 * t55 + t96 * t59 + t74, (t13 * t54 - t58 * t92) * t59 + (t13 * t58 + t54 * t92) * t55 + t74, (t21 * t43 + t27 * t42 + t46 * t69) * r_i_i_C(1) + (t20 * t43 + t26 * t42 - t46 * t70) * r_i_i_C(2); 0, -(pkin(2) * t46 - t93) * t55 + (pkin(2) * t79 + t95) * t59 + t73, t93 * t55 + t95 * t59 + t73, (t14 * t54 + t58 * t64) * t59 + (t14 * t58 - t54 * t64) * t55 + t73, (t23 * t43 + t25 * t42 - t49 * t69) * r_i_i_C(1) + (-t22 * t43 + t24 * t42 + t49 * t70) * r_i_i_C(2); 1, (pkin(2) - t60) * t80 + t63, -t60 * t80 + t63, ((t39 * t54 + t36) * t55 - (t39 * t58 - t83) * t59) * t48 + t66, (r_i_i_C(1) * t56 - r_i_i_C(2) * t52) * t51 * t47 + ((-t42 * t52 + t43 * t77) * r_i_i_C(1) + (-t42 * t56 - t43 * t78) * r_i_i_C(2)) * t48;];
	Ja_transl = t1;
end