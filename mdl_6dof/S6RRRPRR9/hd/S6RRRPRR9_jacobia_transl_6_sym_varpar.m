% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:37
% EndTime: 2019-02-26 22:20:37
% DurationCPUTime: 0.56s
% Computational Cost: add. (842->121), mult. (2275->205), div. (0->0), fcn. (3013->16), ass. (0->72)
t57 = sin(pkin(7));
t56 = sin(pkin(13));
t62 = sin(qJ(3));
t67 = cos(qJ(3));
t81 = cos(pkin(13));
t71 = -t62 * t56 + t67 * t81;
t39 = t71 * t57;
t59 = cos(pkin(7));
t41 = t71 * t59;
t63 = sin(qJ(2));
t64 = sin(qJ(1));
t68 = cos(qJ(2));
t69 = cos(qJ(1));
t82 = cos(pkin(6));
t77 = t69 * t82;
t46 = t64 * t63 - t68 * t77;
t47 = t63 * t77 + t64 * t68;
t51 = -t67 * t56 - t62 * t81;
t58 = sin(pkin(6));
t84 = t58 * t69;
t18 = t39 * t84 + t46 * t41 - t47 * t51;
t40 = t51 * t57;
t42 = t51 * t59;
t19 = -t40 * t84 - t46 * t42 - t47 * t71;
t35 = -t46 * t57 + t59 * t84;
t61 = sin(qJ(5));
t66 = cos(qJ(5));
t6 = t19 * t66 + t35 * t61;
t60 = sin(qJ(6));
t65 = cos(qJ(6));
t96 = t18 * t65 + t6 * t60;
t95 = -t18 * t60 + t6 * t65;
t92 = t19 * t61 - t35 * t66;
t25 = (-t42 * t68 + t63 * t71) * t58 - t82 * t40;
t91 = r_i_i_C(3) + pkin(12);
t90 = pkin(3) * t62;
t89 = t57 * t61;
t88 = t57 * t66;
t87 = t58 * t63;
t86 = t58 * t64;
t85 = t58 * t68;
t83 = pkin(10) + qJ(4);
t80 = t57 * t87;
t79 = t58 * (t57 * t90 + t83 * t59 + pkin(9));
t78 = t64 * t82;
t74 = t65 * r_i_i_C(1) - t60 * r_i_i_C(2) + pkin(5);
t73 = -t60 * r_i_i_C(1) - t65 * r_i_i_C(2) - pkin(11);
t48 = -t69 * t63 - t68 * t78;
t72 = -t48 * t57 + t59 * t86;
t49 = -t63 * t78 + t69 * t68;
t22 = -t40 * t86 - t48 * t42 + t49 * t71;
t70 = t91 * t61 + t74 * t66 + pkin(4);
t55 = t67 * pkin(3) + pkin(2);
t45 = -t57 * t85 + t82 * t59;
t44 = -t83 * t57 + t59 * t90;
t33 = (t42 * t63 + t68 * t71) * t58;
t32 = -t41 * t87 + t51 * t85;
t31 = t33 * t66 + t61 * t80;
t29 = t49 * t42 + t48 * t71;
t28 = -t49 * t41 + t48 * t51;
t27 = t47 * t42 - t46 * t71;
t26 = -t47 * t41 - t46 * t51;
t24 = t82 * t39 + (t41 * t68 + t51 * t63) * t58;
t21 = t39 * t86 + t48 * t41 + t49 * t51;
t14 = t29 * t66 + t49 * t89;
t12 = t27 * t66 + t47 * t89;
t10 = t25 * t66 + t45 * t61;
t8 = t22 * t66 + t72 * t61;
t7 = t22 * t61 - t72 * t66;
t2 = -t21 * t60 + t8 * t65;
t1 = -t21 * t65 - t8 * t60;
t3 = [-t64 * pkin(1) + t19 * pkin(4) + t6 * pkin(5) - t18 * pkin(11) + t95 * r_i_i_C(1) - r_i_i_C(2) * t96 + t46 * t44 - t47 * t55 + t69 * t79 + t91 * t92 (t14 * t65 - t28 * t60) * r_i_i_C(1) + (-t14 * t60 - t28 * t65) * r_i_i_C(2) + t14 * pkin(5) + t29 * pkin(4) - t28 * pkin(11) + t48 * t55 - t49 * t44 + t91 * (t29 * t61 - t49 * t88) -t73 * t22 + (-t49 * t62 + (t48 * t59 + t57 * t86) * t67) * pkin(3) + t70 * t21, t72, -t74 * t7 + t91 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t69 * pkin(1) + t22 * pkin(4) + t8 * pkin(5) - t21 * pkin(11) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t48 * t44 + t49 * t55 + t64 * t79 + t91 * t7 (t12 * t65 - t26 * t60) * r_i_i_C(1) + (-t12 * t60 - t26 * t65) * r_i_i_C(2) + t12 * pkin(5) + t27 * pkin(4) - t26 * pkin(11) - t46 * t55 - t47 * t44 + t91 * (t27 * t61 - t47 * t88) t73 * t19 + (-t47 * t62 + (-t46 * t59 - t57 * t84) * t67) * pkin(3) - t70 * t18, -t35, -t6 * t91 + t74 * t92, r_i_i_C(1) * t96 + t95 * r_i_i_C(2); 0 (t31 * t65 - t32 * t60) * r_i_i_C(1) + (-t31 * t60 - t32 * t65) * r_i_i_C(2) + t31 * pkin(5) + t33 * pkin(4) - t32 * pkin(11) + (-t63 * t44 + t68 * t55) * t58 + t91 * (t33 * t61 - t66 * t80) -t73 * t25 + (t82 * t57 * t67 + (t59 * t67 * t68 - t62 * t63) * t58) * pkin(3) + t70 * t24, t45, t91 * t10 + t74 * (-t25 * t61 + t45 * t66) (-t10 * t60 - t24 * t65) * r_i_i_C(1) + (-t10 * t65 + t24 * t60) * r_i_i_C(2);];
Ja_transl  = t3;
