% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:44
% EndTime: 2019-02-26 22:55:45
% DurationCPUTime: 0.64s
% Computational Cost: add. (666->109), mult. (1904->195), div. (0->0), fcn. (2507->16), ass. (0->79)
t58 = sin(qJ(2));
t59 = sin(qJ(1));
t62 = cos(qJ(2));
t63 = cos(qJ(1));
t79 = cos(pkin(6));
t76 = t63 * t79;
t44 = t58 * t76 + t59 * t62;
t49 = sin(pkin(14));
t53 = cos(pkin(14));
t43 = t59 * t58 - t62 * t76;
t51 = sin(pkin(7));
t55 = cos(pkin(7));
t52 = sin(pkin(6));
t84 = t52 * t63;
t68 = t43 * t55 + t51 * t84;
t26 = t44 * t49 + t68 * t53;
t37 = -t43 * t51 + t55 * t84;
t50 = sin(pkin(8));
t54 = cos(pkin(8));
t15 = t26 * t50 - t37 * t54;
t56 = sin(qJ(5));
t27 = -t44 * t53 + t68 * t49;
t57 = sin(qJ(4));
t61 = cos(qJ(4));
t73 = t26 * t54 + t37 * t50;
t6 = t27 * t61 + t73 * t57;
t60 = cos(qJ(5));
t104 = t15 * t60 + t6 * t56;
t103 = -t15 * t56 + t6 * t60;
t100 = t27 * t57 - t73 * t61;
t77 = t59 * t79;
t45 = -t63 * t58 - t62 * t77;
t85 = t52 * t59;
t39 = -t45 * t51 + t55 * t85;
t46 = -t58 * t77 + t63 * t62;
t67 = t45 * t55 + t51 * t85;
t65 = t46 * t49 - t67 * t53;
t95 = -t39 * t50 + t65 * t54;
t94 = r_i_i_C(3) + pkin(12);
t82 = t55 * t58;
t40 = (-t49 * t62 - t53 * t82) * t52;
t92 = t40 * t50;
t89 = t49 * t55;
t88 = t50 * t51;
t87 = t51 * t54;
t86 = t52 * t58;
t83 = t53 * t55;
t81 = t55 * t62;
t80 = t51 * qJ(3);
t78 = t51 * t86;
t75 = t79 * t51;
t35 = t53 * t75 + (-t49 * t58 + t53 * t81) * t52;
t42 = -t52 * t62 * t51 + t79 * t55;
t72 = t35 * t54 + t42 * t50;
t71 = t60 * r_i_i_C(1) - t56 * r_i_i_C(2) + pkin(4);
t29 = t43 * t49 - t44 * t83;
t20 = -t29 * t50 + t44 * t87;
t70 = t29 * t54 + t44 * t88;
t31 = -t45 * t49 - t46 * t83;
t21 = -t31 * t50 + t46 * t87;
t69 = t31 * t54 + t46 * t88;
t66 = t40 * t54 + t50 * t78;
t17 = t39 * t54 + t65 * t50;
t41 = (-t49 * t82 + t53 * t62) * t52;
t36 = t53 * t86 + (t52 * t81 + t75) * t49;
t33 = t54 * t78 - t92;
t32 = t45 * t53 - t46 * t89;
t30 = -t43 * t53 - t44 * t89;
t28 = t46 * t53 + t67 * t49;
t23 = -t35 * t50 + t42 * t54;
t19 = t41 * t61 + t66 * t57;
t14 = t36 * t61 + t72 * t57;
t12 = t32 * t61 + t69 * t57;
t10 = t30 * t61 + t70 * t57;
t8 = t28 * t61 - t95 * t57;
t7 = t28 * t57 + t95 * t61;
t2 = t17 * t56 + t8 * t60;
t1 = t17 * t60 - t8 * t56;
t3 = [-t59 * pkin(1) - t44 * pkin(2) + t27 * pkin(3) + t6 * pkin(4) + pkin(10) * t84 - t15 * pkin(11) + t103 * r_i_i_C(1) - t104 * r_i_i_C(2) + t37 * qJ(3) + t94 * t100 (t12 * t60 + t21 * t56) * r_i_i_C(1) + (-t12 * t56 + t21 * t60) * r_i_i_C(2) + t12 * pkin(4) + t32 * pkin(3) + t45 * pkin(2) + t46 * t80 + t94 * (t32 * t57 - t69 * t61) + t21 * pkin(11), t39, -t71 * t7 + t94 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t63 * pkin(1) + t46 * pkin(2) + t28 * pkin(3) + t8 * pkin(4) + pkin(10) * t85 + t17 * pkin(11) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t39 * qJ(3) + t94 * t7 (t10 * t60 + t20 * t56) * r_i_i_C(1) + (-t10 * t56 + t20 * t60) * r_i_i_C(2) + t10 * pkin(4) + t30 * pkin(3) - t43 * pkin(2) + t44 * t80 + t94 * (t30 * t57 - t70 * t61) + t20 * pkin(11), -t37, t71 * t100 - t6 * t94, t104 * r_i_i_C(1) + t103 * r_i_i_C(2), 0; 0 (t19 * t60 + t33 * t56) * r_i_i_C(1) + (-t19 * t56 + t33 * t60) * r_i_i_C(2) + t19 * pkin(4) + t41 * pkin(3) - pkin(11) * t92 + t94 * (t41 * t57 - t66 * t61) + (t62 * pkin(2) + (pkin(11) * t54 + qJ(3)) * t58 * t51) * t52, t42, t94 * t14 + t71 * (-t36 * t57 + t72 * t61) (-t14 * t56 + t23 * t60) * r_i_i_C(1) + (-t14 * t60 - t23 * t56) * r_i_i_C(2), 0;];
Ja_transl  = t3;
