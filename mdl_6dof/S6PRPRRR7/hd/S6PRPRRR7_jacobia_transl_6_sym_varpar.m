% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobia_transl_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:29
% EndTime: 2019-02-26 19:57:30
% DurationCPUTime: 0.50s
% Computational Cost: add. (1087->120), mult. (3155->222), div. (0->0), fcn. (4179->18), ass. (0->89)
t50 = sin(pkin(8));
t54 = cos(pkin(8));
t53 = cos(pkin(13));
t58 = sin(qJ(2));
t61 = cos(qJ(2));
t93 = sin(pkin(13));
t96 = cos(pkin(6));
t82 = t96 * t93;
t47 = t53 * t61 - t58 * t82;
t49 = sin(pkin(14));
t76 = t53 * t58 + t61 * t82;
t94 = cos(pkin(14));
t73 = t76 * t94;
t51 = sin(pkin(7));
t52 = sin(pkin(6));
t87 = t52 * t93;
t83 = t51 * t87;
t95 = cos(pkin(7));
t63 = t47 * t49 + t95 * t73 - t94 * t83;
t69 = t76 * t51 + t95 * t87;
t107 = -t69 * t50 + t63 * t54;
t86 = t53 * t96;
t46 = t58 * t86 + t93 * t61;
t77 = t93 * t58 - t61 * t86;
t72 = t77 * t94;
t88 = t52 * t94;
t65 = t53 * t51 * t88 + t46 * t49 + t95 * t72;
t89 = t52 * t95;
t71 = t77 * t51 - t53 * t89;
t106 = -t71 * t50 + t65 * t54;
t81 = t95 * t94;
t85 = t96 * t51;
t68 = t94 * t85 + (-t49 * t58 + t61 * t81) * t52;
t100 = t51 * t52;
t78 = -t61 * t100 + t96 * t95;
t105 = t78 * t50 + t68 * t54;
t104 = r_i_i_C(3) + pkin(12);
t103 = cos(qJ(4));
t44 = (-t49 * t61 - t58 * t81) * t52;
t102 = t44 * t50;
t101 = t50 * t51;
t99 = t51 * t54;
t98 = t58 * t51;
t97 = t51 * qJ(3);
t92 = t52 * t98;
t91 = t54 * t103;
t90 = t49 * t95;
t84 = t103 * t101;
t55 = sin(qJ(6));
t59 = cos(qJ(6));
t80 = t59 * r_i_i_C(1) - t55 * r_i_i_C(2) + pkin(5);
t79 = t55 * r_i_i_C(1) + t59 * r_i_i_C(2) + pkin(11);
t34 = -t46 * t81 + t77 * t49;
t25 = -t34 * t50 + t46 * t99;
t36 = -t47 * t81 + t76 * t49;
t26 = -t36 * t50 + t47 * t99;
t56 = sin(qJ(5));
t60 = cos(qJ(5));
t74 = -t104 * t56 - t80 * t60 - pkin(4);
t57 = sin(qJ(4));
t45 = (-t58 * t90 + t94 * t61) * t52;
t42 = t58 * t88 + (t61 * t89 + t85) * t49;
t38 = t54 * t92 - t102;
t37 = -t47 * t90 - t73;
t35 = -t46 * t90 - t72;
t33 = t47 * t94 + (-t76 * t95 + t83) * t49;
t32 = t46 * t94 + (-t53 * t100 - t77 * t95) * t49;
t31 = -t68 * t50 + t78 * t54;
t28 = t45 * t103 + (t44 * t54 + t50 * t92) * t57;
t27 = -t52 * t58 * t84 - t44 * t91 + t45 * t57;
t24 = t63 * t50 + t69 * t54;
t23 = t65 * t50 + t71 * t54;
t22 = t42 * t103 + t105 * t57;
t21 = -t105 * t103 + t42 * t57;
t20 = t28 * t60 + t38 * t56;
t18 = t37 * t103 + (t47 * t101 + t36 * t54) * t57;
t17 = -t36 * t91 + t37 * t57 - t47 * t84;
t16 = t35 * t103 + (t46 * t101 + t34 * t54) * t57;
t15 = -t34 * t91 + t35 * t57 - t46 * t84;
t14 = t33 * t103 - t107 * t57;
t13 = t107 * t103 + t33 * t57;
t12 = t32 * t103 - t106 * t57;
t11 = t106 * t103 + t32 * t57;
t10 = t22 * t60 + t31 * t56;
t8 = t18 * t60 + t26 * t56;
t6 = t16 * t60 + t25 * t56;
t4 = t14 * t60 + t24 * t56;
t2 = t12 * t60 + t23 * t56;
t1 = [0 (t17 * t55 + t8 * t59) * r_i_i_C(1) + (t17 * t59 - t8 * t55) * r_i_i_C(2) + t8 * pkin(5) + t18 * pkin(4) + t17 * pkin(11) + t37 * pkin(3) - t76 * pkin(2) + t47 * t97 + t104 * (t18 * t56 - t26 * t60) + t26 * pkin(10), t69, t74 * t13 + t79 * t14, t104 * t4 + t80 * (-t14 * t56 + t24 * t60) (t13 * t59 - t4 * t55) * r_i_i_C(1) + (-t13 * t55 - t4 * t59) * r_i_i_C(2); 0 (t15 * t55 + t6 * t59) * r_i_i_C(1) + (t15 * t59 - t6 * t55) * r_i_i_C(2) + t6 * pkin(5) + t16 * pkin(4) + t15 * pkin(11) + t35 * pkin(3) - t77 * pkin(2) + t46 * t97 + t104 * (t16 * t56 - t25 * t60) + t25 * pkin(10), t71, t74 * t11 + t79 * t12, t104 * t2 + t80 * (-t12 * t56 + t23 * t60) (t11 * t59 - t2 * t55) * r_i_i_C(1) + (-t11 * t55 - t2 * t59) * r_i_i_C(2); 1 (t20 * t59 + t27 * t55) * r_i_i_C(1) + (-t20 * t55 + t27 * t59) * r_i_i_C(2) + t20 * pkin(5) + t28 * pkin(4) + t27 * pkin(11) + t45 * pkin(3) - pkin(10) * t102 + t104 * (t28 * t56 - t38 * t60) + (t61 * pkin(2) + (pkin(10) * t54 + qJ(3)) * t98) * t52, t78, t74 * t21 + t79 * t22, t104 * t10 + t80 * (-t22 * t56 + t31 * t60) (-t10 * t55 + t21 * t59) * r_i_i_C(1) + (-t10 * t59 - t21 * t55) * r_i_i_C(2);];
Ja_transl  = t1;
