% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:05
% EndTime: 2019-02-26 22:52:05
% DurationCPUTime: 0.38s
% Computational Cost: add. (751->99), mult. (1867->170), div. (0->0), fcn. (2438->16), ass. (0->69)
t61 = sin(qJ(2));
t64 = cos(qJ(2));
t65 = cos(qJ(1));
t84 = cos(pkin(6));
t76 = t65 * t84;
t92 = sin(qJ(1));
t43 = t92 * t61 - t64 * t76;
t44 = t61 * t76 + t92 * t64;
t60 = sin(qJ(3));
t55 = sin(pkin(7));
t56 = sin(pkin(6));
t89 = t56 * t65;
t82 = t55 * t89;
t57 = cos(pkin(7));
t88 = t57 * t60;
t93 = cos(qJ(3));
t22 = -t43 * t88 + t44 * t93 - t60 * t82;
t36 = -t43 * t55 + t57 * t89;
t59 = sin(qJ(4));
t63 = cos(qJ(4));
t10 = t22 * t63 - t36 * t59;
t99 = -t22 * t59 - t36 * t63;
t79 = t57 * t93;
t21 = t43 * t79 + t44 * t60 + t93 * t82;
t54 = qJ(5) + qJ(6);
t52 = sin(t54);
t53 = cos(t54);
t98 = (-t10 * t52 + t21 * t53) * r_i_i_C(1) + (-t10 * t53 - t21 * t52) * r_i_i_C(2);
t73 = t84 * t92;
t45 = -t61 * t73 + t65 * t64;
t70 = t65 * t61 + t64 * t73;
t80 = t56 * t92;
t74 = t55 * t80;
t26 = t45 * t93 + (-t70 * t57 + t74) * t60;
t67 = t70 * t55 + t57 * t80;
t14 = t26 * t63 + t67 * t59;
t69 = t70 * t93;
t25 = t45 * t60 + t57 * t69 - t93 * t74;
t5 = -t14 * t52 + t25 * t53;
t6 = t14 * t53 + t25 * t52;
t97 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t75 = t84 * t55;
t78 = t93 * t61;
t86 = t60 * t64;
t35 = t60 * t75 + (t57 * t86 + t78) * t56;
t42 = -t56 * t64 * t55 + t84 * t57;
t20 = t35 * t63 + t42 * t59;
t77 = t93 * t64;
t87 = t60 * t61;
t34 = -t93 * t75 + (-t57 * t77 + t87) * t56;
t96 = (-t20 * t52 + t34 * t53) * r_i_i_C(1) + (-t20 * t53 - t34 * t52) * r_i_i_C(2);
t95 = pkin(10) * t55;
t94 = r_i_i_C(3) + pkin(13) + pkin(12);
t91 = t55 * t59;
t90 = t55 * t63;
t85 = t61 * t55;
t83 = t56 * t85;
t58 = sin(qJ(5));
t81 = t58 * pkin(5) + pkin(11);
t62 = cos(qJ(5));
t51 = t62 * pkin(5) + pkin(4);
t72 = r_i_i_C(1) * t53 - r_i_i_C(2) * t52 + t51;
t71 = t52 * r_i_i_C(1) + t53 * r_i_i_C(2) + t81;
t68 = -t94 * t59 - t72 * t63 - pkin(3);
t41 = (-t57 * t87 + t77) * t56;
t30 = -t45 * t88 - t69;
t28 = -t43 * t93 - t44 * t88;
t13 = t26 * t59 - t67 * t63;
t1 = [-t92 * pkin(1) - t44 * pkin(2) - t22 * pkin(3) + pkin(9) * t89 + t36 * pkin(10) - t72 * t10 - t71 * t21 + t94 * t99, t30 * pkin(3) - t70 * pkin(2) + t45 * t95 + t72 * (t30 * t63 + t45 * t91) + t71 * (t45 * t79 - t70 * t60) + t94 * (t30 * t59 - t45 * t90) t68 * t25 + t71 * t26, -t72 * t13 + t94 * t14 (-t14 * t58 + t25 * t62) * pkin(5) + t97, t97; t65 * pkin(1) + t45 * pkin(2) + t26 * pkin(3) + pkin(9) * t80 + t67 * pkin(10) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t94 * t13 + t14 * t51 + t81 * t25, t44 * t95 - t43 * pkin(2) + t28 * pkin(3) + t72 * (t28 * t63 + t44 * t91) + t71 * (-t43 * t60 + t44 * t79) + t94 * (t28 * t59 - t44 * t90) t68 * t21 + t71 * t22, t94 * t10 + t72 * t99 (-t10 * t58 + t21 * t62) * pkin(5) + t98, t98; 0, t41 * pkin(3) + t72 * (t41 * t63 + t59 * t83) + t94 * (t41 * t59 - t63 * t83) + (t64 * pkin(2) + pkin(10) * t85 + t71 * (t57 * t78 + t86)) * t56, t68 * t34 + t71 * t35, t94 * t20 + t72 * (-t35 * t59 + t42 * t63) (-t20 * t58 + t34 * t62) * pkin(5) + t96, t96;];
Ja_transl  = t1;
