% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR15
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRR15_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR15_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR15_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
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
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (129->29), mult. (88->42), div. (0->0), fcn. (96->12), ass. (0->28)
	t27 = cos(qJ(1));
	t35 = -t27 / 0.2e1;
	t21 = qJ(2) + qJ(3);
	t17 = pkin(5) + t21;
	t12 = sin(t17);
	t18 = pkin(5) - t21;
	t13 = sin(t18);
	t10 = t12 - t13;
	t20 = cos(t21);
	t25 = sin(qJ(1));
	t14 = cos(t17);
	t15 = cos(t18);
	t11 = t15 + t14;
	t19 = sin(t21);
	t5 = t11 * t35 + t25 * t19;
	t34 = -t5 * r_i_i_C(1) + (t10 * t35 - t25 * t20) * r_i_i_C(2);
	t6 = -t27 * t19 - t25 * t11 / 0.2e1;
	t33 = t6 * r_i_i_C(1) + (-t27 * t20 + t25 * t10 / 0.2e1) * r_i_i_C(2);
	t32 = (t12 / 0.2e1 + t13 / 0.2e1) * r_i_i_C(1) + (t14 / 0.2e1 - t15 / 0.2e1) * r_i_i_C(2);
	t26 = cos(qJ(2));
	t31 = t26 * pkin(2);
	t23 = cos(pkin(5));
	t30 = t23 * t26;
	t29 = r_i_i_C(1) * t20 + pkin(1) + t31;
	t22 = sin(pkin(5));
	t24 = sin(qJ(2));
	t28 = -t10 * r_i_i_C(1) / 0.2e1 - t23 * t24 * pkin(2) + (r_i_i_C(3) + pkin(8) + pkin(9)) * t22;
	t1 = [t5 * r_i_i_C(2) - t29 * t25 + t28 * t27, (-t24 * t27 - t25 * t30) * pkin(2) + t33, t33, 0, 0; t6 * r_i_i_C(2) + t28 * t25 + t29 * t27, (-t24 * t25 + t27 * t30) * pkin(2) + t34, t34, 0, 0; 0, t22 * t31 + t32, t32, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (260->42), mult. (174->55), div. (0->0), fcn. (168->16), ass. (0->36)
	t26 = qJ(2) + qJ(3);
	t25 = qJ(4) + t26;
	t39 = pkin(5) - t25;
	t36 = sin(t39);
	t20 = pkin(5) + t25;
	t41 = sin(t20) / 0.2e1;
	t12 = t41 - t36 / 0.2e1;
	t22 = cos(t25);
	t31 = sin(qJ(1));
	t34 = cos(qJ(1));
	t18 = cos(t20) / 0.2e1;
	t19 = cos(t39);
	t13 = t19 / 0.2e1 + t18;
	t21 = sin(t25);
	t5 = -t34 * t13 + t31 * t21;
	t46 = -t5 * r_i_i_C(1) + (-t34 * t12 - t31 * t22) * r_i_i_C(2);
	t6 = -t31 * t13 - t34 * t21;
	t45 = t6 * r_i_i_C(1) + (t31 * t12 - t34 * t22) * r_i_i_C(2);
	t44 = pkin(3) * sin(t26);
	t29 = sin(qJ(3));
	t30 = sin(qJ(2));
	t43 = t29 * t30;
	t42 = (t41 + t36 / 0.2e1) * r_i_i_C(1) + (t18 - t19 / 0.2e1) * r_i_i_C(2);
	t33 = cos(qJ(2));
	t40 = t22 * r_i_i_C(1) + pkin(1) + pkin(3) * cos(t26) + t33 * pkin(2);
	t32 = cos(qJ(3));
	t23 = t32 * pkin(3) + pkin(2);
	t27 = sin(pkin(5));
	t28 = cos(pkin(5));
	t38 = -t12 * r_i_i_C(1) - (t33 * t29 * pkin(3) + t30 * t23) * t28 + (r_i_i_C(3) + pkin(8) + pkin(9) + pkin(10)) * t27;
	t37 = -pkin(3) * t43 + t23 * t33;
	t35 = pkin(3) * (t32 * t33 - t43);
	t15 = -t30 * pkin(2) - t44;
	t9 = t28 * t35;
	t8 = t37 * t28;
	t1 = [t5 * r_i_i_C(2) - t40 * t31 + t38 * t34, t34 * t15 - t31 * t8 + t45, -t31 * t9 - t34 * t44 + t45, t45, 0; t6 * r_i_i_C(2) + t38 * t31 + t40 * t34, t31 * t15 + t34 * t8 + t46, -t31 * t44 + t34 * t9 + t46, t46, 0; 0, t37 * t27 + t42, t27 * t35 + t42, t42, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (462->83), mult. (741->134), div. (0->0), fcn. (829->16), ass. (0->76)
	t48 = sin(qJ(4));
	t53 = cos(qJ(4));
	t43 = sin(pkin(6));
	t89 = pkin(11) * t43;
	t38 = -t48 * pkin(4) + t53 * t89;
	t54 = cos(qJ(3));
	t34 = t38 * t54;
	t37 = -pkin(4) * t53 - t48 * t89;
	t36 = pkin(3) - t37;
	t49 = sin(qJ(3));
	t19 = -t49 * t36 + t34;
	t55 = cos(qJ(2));
	t16 = t19 * t55;
	t85 = t38 * t49;
	t18 = -t36 * t54 - t85;
	t17 = pkin(2) - t18;
	t50 = sin(qJ(2));
	t6 = -t17 * t50 + t16;
	t88 = r_i_i_C(3) * t43;
	t86 = t19 * t50;
	t42 = qJ(2) + qJ(3) + qJ(4);
	t40 = sin(t42);
	t46 = cos(pkin(5));
	t84 = t40 * t46;
	t44 = sin(pkin(5));
	t83 = t43 * t44;
	t82 = t43 * t46;
	t81 = t44 * t55;
	t45 = cos(pkin(6));
	t47 = sin(qJ(5));
	t80 = t45 * t47;
	t52 = cos(qJ(5));
	t79 = t45 * t52;
	t51 = sin(qJ(1));
	t78 = t51 * t47;
	t77 = t52 * t51;
	t56 = cos(qJ(1));
	t76 = t56 * t47;
	t75 = t56 * t52;
	t74 = t40 * t88;
	t73 = t47 * t83;
	t72 = t52 * t83;
	t71 = t56 * t83;
	t70 = t46 * t76;
	t69 = t46 * t78;
	t68 = t46 * t77;
	t67 = t46 * t75;
	t26 = t45 * t69 - t75;
	t27 = -t45 * t68 - t76;
	t32 = t45 * t76 + t68;
	t33 = -t45 * t75 + t69;
	t41 = cos(t42);
	t66 = (t26 * t40 - t32 * t41) * r_i_i_C(1) + (-t27 * t40 + t33 * t41) * r_i_i_C(2) + (t41 * t56 - t51 * t84) * t88;
	t28 = t45 * t70 + t77;
	t29 = t45 * t67 - t78;
	t30 = t45 * t78 - t67;
	t31 = t45 * t77 + t70;
	t65 = (-t28 * t40 - t30 * t41) * r_i_i_C(1) + (-t29 * t40 - t31 * t41) * r_i_i_C(2) + (t41 * t51 + t56 * t84) * t88;
	t64 = ((-t40 * t79 - t41 * t47) * r_i_i_C(2) + (-t40 * t80 + t41 * t52) * r_i_i_C(1) + t74) * t44;
	t63 = t44 * t86 + t64;
	t61 = t17 * t55 + t86;
	t62 = pkin(1) + t61 + t74;
	t20 = t37 * t54 - t85;
	t21 = t37 * t49 + t34;
	t60 = t20 * t55 - t21 * t50;
	t59 = -t26 * t41 - t32 * t40;
	t58 = -t29 * t41 + t31 * t40;
	t57 = t44 * (t45 * pkin(11) + pkin(8) + pkin(9) + pkin(10)) + t6 * t46 + (t41 * t82 + t44 * t45) * r_i_i_C(3);
	t10 = -t28 * t41 + t30 * t40 + t47 * t71;
	t9 = t27 * t41 + t33 * t40 + t51 * t72;
	t8 = t20 * t50 + t21 * t55;
	t7 = t18 * t50 + t16;
	t4 = t60 * t46;
	t3 = (t18 * t55 - t86) * t46;
	t2 = t61 * t46;
	t1 = [t10 * r_i_i_C(1) + t58 * r_i_i_C(2) - t62 * t51 + (r_i_i_C(2) * t72 + t57) * t56, -t2 * t51 + t6 * t56 + t66, t3 * t51 + t7 * t56 + t66, t4 * t51 + t8 * t56 + t66, t9 * r_i_i_C(1) + (-t51 * t73 - t59) * r_i_i_C(2); t59 * r_i_i_C(1) + t9 * r_i_i_C(2) + t62 * t56 + (r_i_i_C(1) * t73 + t57) * t51, t2 * t56 + t6 * t51 + t65, -t3 * t56 + t7 * t51 + t65, -t4 * t56 + t8 * t51 + t65, (-t52 * t71 - t58) * r_i_i_C(1) + t10 * r_i_i_C(2); 0, t17 * t81 + t63, -t18 * t81 + t63, -t60 * t44 + t64, (r_i_i_C(1) * t52 - r_i_i_C(2) * t47) * t82 + ((-t40 * t47 + t41 * t79) * r_i_i_C(1) + (-t40 * t52 - t41 * t80) * r_i_i_C(2)) * t44;];
	Ja_transl = t1;
end