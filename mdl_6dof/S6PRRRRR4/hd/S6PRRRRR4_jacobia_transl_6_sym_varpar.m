% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR4
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:20:36
% EndTime: 2019-02-26 20:20:36
% DurationCPUTime: 0.35s
% Computational Cost: add. (660->94), mult. (1496->168), div. (0->0), fcn. (1949->16), ass. (0->67)
t90 = r_i_i_C(3) + pkin(12);
t52 = sin(qJ(6));
t56 = cos(qJ(6));
t93 = t56 * r_i_i_C(1) - t52 * r_i_i_C(2) + pkin(5);
t49 = sin(pkin(7));
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t51 = cos(pkin(13));
t81 = cos(pkin(6));
t75 = t51 * t81;
t79 = sin(pkin(13));
t63 = t79 * t55 - t58 * t75;
t80 = cos(pkin(7));
t50 = sin(pkin(6));
t84 = t50 * t51;
t92 = t49 * t84 + t63 * t80;
t69 = t81 * t79;
t64 = t51 * t55 + t58 * t69;
t76 = t50 * t79;
t91 = -t49 * t76 + t64 * t80;
t87 = cos(qJ(3));
t48 = qJ(4) + qJ(5);
t46 = sin(t48);
t86 = t46 * t49;
t47 = cos(t48);
t85 = t47 * t49;
t83 = t50 * t55;
t82 = t50 * t58;
t78 = t49 * t83;
t54 = sin(qJ(3));
t74 = t54 * t80;
t73 = t81 * t49;
t53 = sin(qJ(4));
t71 = (pkin(4) * t53 + pkin(9)) * t49;
t70 = t80 * t87;
t59 = -pkin(11) - pkin(10);
t68 = t52 * r_i_i_C(1) + t56 * r_i_i_C(2) - t59;
t39 = t55 * t75 + t79 * t58;
t21 = t39 * t87 - t92 * t54;
t32 = t63 * t49 - t80 * t84;
t8 = t21 * t47 + t32 * t46;
t67 = t90 * t8 + t93 * (-t21 * t46 + t32 * t47);
t40 = t51 * t58 - t55 * t69;
t23 = t40 * t87 - t91 * t54;
t33 = t64 * t49 + t80 * t76;
t10 = t23 * t47 + t33 * t46;
t66 = t90 * t10 + t93 * (-t23 * t46 + t33 * t47);
t31 = t54 * t73 + (t87 * t55 + t58 * t74) * t50;
t38 = -t49 * t82 + t81 * t80;
t19 = t31 * t47 + t38 * t46;
t65 = t90 * t19 + t93 * (-t31 * t46 + t38 * t47);
t57 = cos(qJ(4));
t45 = t57 * pkin(4) + pkin(3);
t62 = -t90 * t46 - t47 * t93 - t45;
t37 = (-t55 * t74 + t87 * t58) * t50;
t36 = (t54 * t58 + t55 * t70) * t50;
t30 = t54 * t83 - t70 * t82 - t87 * t73;
t29 = t37 * t47 + t46 * t78;
t27 = -t40 * t74 - t64 * t87;
t26 = t40 * t70 - t64 * t54;
t25 = -t39 * t74 - t63 * t87;
t24 = t39 * t70 - t63 * t54;
t22 = t40 * t54 + t91 * t87;
t20 = t39 * t54 + t92 * t87;
t14 = t27 * t47 + t40 * t86;
t12 = t25 * t47 + t39 * t86;
t1 = [0 (t14 * t56 + t26 * t52) * r_i_i_C(1) + (-t14 * t52 + t26 * t56) * r_i_i_C(2) + t14 * pkin(5) + t27 * t45 - t26 * t59 - t64 * pkin(2) + t40 * t71 + t90 * (t27 * t46 - t40 * t85) t62 * t22 + t68 * t23 (-t23 * t53 + t33 * t57) * pkin(4) + t66, t66 (-t10 * t52 + t22 * t56) * r_i_i_C(1) + (-t10 * t56 - t22 * t52) * r_i_i_C(2); 0 (t12 * t56 + t24 * t52) * r_i_i_C(1) + (-t12 * t52 + t24 * t56) * r_i_i_C(2) + t12 * pkin(5) + t25 * t45 - t24 * t59 - t63 * pkin(2) + t39 * t71 + t90 * (t25 * t46 - t39 * t85) t62 * t20 + t68 * t21 (-t21 * t53 + t32 * t57) * pkin(4) + t67, t67 (t20 * t56 - t8 * t52) * r_i_i_C(1) + (-t20 * t52 - t8 * t56) * r_i_i_C(2); 1 (t29 * t56 + t36 * t52) * r_i_i_C(1) + (-t29 * t52 + t36 * t56) * r_i_i_C(2) + t29 * pkin(5) + t37 * t45 - t36 * t59 + t90 * (t37 * t46 - t47 * t78) + (t58 * pkin(2) + t55 * t71) * t50, t62 * t30 + t68 * t31 (-t31 * t53 + t38 * t57) * pkin(4) + t65, t65 (-t19 * t52 + t30 * t56) * r_i_i_C(1) + (-t19 * t56 - t30 * t52) * r_i_i_C(2);];
Ja_transl  = t1;
