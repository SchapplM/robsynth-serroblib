% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR10_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:46
% EndTime: 2019-02-26 22:52:47
% DurationCPUTime: 0.82s
% Computational Cost: add. (840->137), mult. (2410->242), div. (0->0), fcn. (3166->16), ass. (0->89)
t64 = sin(qJ(2));
t65 = sin(qJ(1));
t69 = cos(qJ(2));
t70 = cos(qJ(1));
t87 = cos(pkin(6));
t83 = t70 * t87;
t49 = t65 * t64 - t69 * t83;
t50 = t64 * t83 + t65 * t69;
t60 = cos(pkin(7));
t63 = sin(qJ(3));
t68 = cos(qJ(3));
t57 = sin(pkin(7));
t58 = sin(pkin(6));
t96 = t58 * t70;
t85 = t57 * t96;
t32 = t68 * (t49 * t60 + t85) + t50 * t63;
t44 = -t49 * t57 + t60 * t96;
t56 = sin(pkin(8));
t59 = cos(pkin(8));
t19 = t32 * t56 - t44 * t59;
t93 = t60 * t63;
t33 = t49 * t93 - t50 * t68 + t63 * t85;
t62 = sin(qJ(4));
t67 = cos(qJ(4));
t80 = t32 * t59 + t44 * t56;
t6 = t33 * t67 + t62 * t80;
t61 = sin(qJ(5));
t66 = cos(qJ(5));
t116 = t19 * t66 + t6 * t61;
t115 = -t19 * t61 + t6 * t66;
t112 = t33 * t62 - t67 * t80;
t107 = r_i_i_C(3) + pkin(13);
t106 = t56 * pkin(12);
t105 = t57 * pkin(11);
t89 = t64 * t68;
t90 = t63 * t69;
t46 = (-t60 * t89 - t90) * t58;
t104 = t46 * t56;
t102 = t56 * t57;
t101 = t56 * t61;
t100 = t56 * t66;
t99 = t57 * t59;
t98 = t57 * t64;
t97 = t58 * t65;
t95 = t59 * t62;
t94 = t59 * t67;
t92 = t60 * t68;
t91 = t63 * t64;
t88 = t68 * t69;
t86 = t58 * t98;
t84 = t65 * t87;
t82 = t87 * t57;
t42 = t68 * t82 + (t60 * t88 - t91) * t58;
t48 = -t58 * t69 * t57 + t60 * t87;
t79 = t42 * t59 + t48 * t56;
t78 = t66 * r_i_i_C(1) - t61 * r_i_i_C(2) + pkin(4);
t36 = t49 * t63 - t50 * t92;
t26 = -t36 * t56 + t50 * t99;
t77 = t102 * t50 + t36 * t59;
t51 = -t70 * t64 - t69 * t84;
t52 = -t64 * t84 + t70 * t69;
t38 = -t51 * t63 - t52 * t92;
t27 = -t38 * t56 + t52 * t99;
t76 = t102 * t52 + t38 * t59;
t74 = -t51 * t57 + t60 * t97;
t73 = t51 * t60 + t57 * t97;
t72 = t46 * t59 + t56 * t86;
t71 = t74 * t56;
t34 = -t52 * t63 + t68 * t73;
t21 = -t34 * t56 + t59 * t74;
t47 = (-t60 * t91 + t88) * t58;
t43 = t63 * t82 + (t60 * t90 + t89) * t58;
t40 = t59 * t86 - t104;
t39 = t51 * t68 - t52 * t93;
t37 = -t49 * t68 - t50 * t93;
t35 = t52 * t68 + t63 * t73;
t29 = -t42 * t56 + t48 * t59;
t25 = t47 * t67 + t62 * t72;
t23 = t42 * t67 - t43 * t95;
t18 = t43 * t67 + t62 * t79;
t16 = t34 * t67 - t35 * t95;
t14 = -t32 * t67 + t33 * t95;
t12 = t39 * t67 + t62 * t76;
t10 = t37 * t67 + t62 * t77;
t8 = t35 * t67 + (t34 * t59 + t71) * t62;
t7 = -t34 * t94 + t35 * t62 - t67 * t71;
t2 = t21 * t61 + t8 * t66;
t1 = t21 * t66 - t8 * t61;
t3 = [-t65 * pkin(1) - t50 * pkin(2) + t33 * pkin(3) + t6 * pkin(4) + pkin(10) * t96 + t44 * pkin(11) - t19 * pkin(12) + t115 * r_i_i_C(1) - t116 * r_i_i_C(2) + t107 * t112 (t12 * t66 + t27 * t61) * r_i_i_C(1) + (-t12 * t61 + t27 * t66) * r_i_i_C(2) + t12 * pkin(4) + t39 * pkin(3) + t51 * pkin(2) + t52 * t105 + t107 * (t39 * t62 - t67 * t76) + t27 * pkin(12) (t101 * t35 + t16 * t66) * r_i_i_C(1) + (t100 * t35 - t16 * t61) * r_i_i_C(2) + t16 * pkin(4) + t34 * pkin(3) + t35 * t106 + t107 * (t34 * t62 + t35 * t94) t107 * t8 - t78 * t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t70 * pkin(1) + t52 * pkin(2) + t35 * pkin(3) + t8 * pkin(4) + pkin(10) * t97 + t74 * pkin(11) + t21 * pkin(12) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t107 * t7 (t10 * t66 + t26 * t61) * r_i_i_C(1) + (-t10 * t61 + t26 * t66) * r_i_i_C(2) + t10 * pkin(4) + t37 * pkin(3) - t49 * pkin(2) + t50 * t105 + t107 * (t37 * t62 - t67 * t77) + t26 * pkin(12) (-t101 * t33 + t14 * t66) * r_i_i_C(1) + (-t100 * t33 - t14 * t61) * r_i_i_C(2) + t14 * pkin(4) - t32 * pkin(3) - t33 * t106 + t107 * (-t32 * t62 - t33 * t94) -t107 * t6 + t78 * t112, t116 * r_i_i_C(1) + t115 * r_i_i_C(2), 0; 0 (t25 * t66 + t40 * t61) * r_i_i_C(1) + (-t25 * t61 + t40 * t66) * r_i_i_C(2) + t25 * pkin(4) + t47 * pkin(3) - pkin(12) * t104 + t107 * (t47 * t62 - t67 * t72) + (t69 * pkin(2) + (pkin(12) * t59 + pkin(11)) * t98) * t58 (t101 * t43 + t23 * t66) * r_i_i_C(1) + (t100 * t43 - t23 * t61) * r_i_i_C(2) + t23 * pkin(4) + t42 * pkin(3) + t43 * t106 + t107 * (t42 * t62 + t43 * t94) t107 * t18 + t78 * (-t43 * t62 + t67 * t79) (-t18 * t61 + t29 * t66) * r_i_i_C(1) + (-t18 * t66 - t29 * t61) * r_i_i_C(2), 0;];
Ja_transl  = t3;
