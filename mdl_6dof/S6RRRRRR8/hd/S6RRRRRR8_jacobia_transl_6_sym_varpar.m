% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR8
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
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:22
% EndTime: 2019-02-26 22:51:23
% DurationCPUTime: 0.54s
% Computational Cost: add. (822->114), mult. (1862->190), div. (0->0), fcn. (2429->16), ass. (0->78)
t64 = sin(qJ(2));
t67 = cos(qJ(2));
t68 = cos(qJ(1));
t89 = cos(pkin(6));
t82 = t68 * t89;
t97 = sin(qJ(1));
t46 = t64 * t97 - t67 * t82;
t47 = t64 * t82 + t67 * t97;
t63 = sin(qJ(3));
t58 = sin(pkin(7));
t59 = sin(pkin(6));
t93 = t59 * t68;
t87 = t58 * t93;
t60 = cos(pkin(7));
t92 = t60 * t63;
t98 = cos(qJ(3));
t25 = -t46 * t92 + t47 * t98 - t63 * t87;
t39 = -t46 * t58 + t60 * t93;
t57 = qJ(4) + qJ(5);
t55 = sin(t57);
t56 = cos(t57);
t10 = t25 * t56 - t39 * t55;
t85 = t60 * t98;
t24 = t46 * t85 + t47 * t63 + t87 * t98;
t61 = sin(qJ(6));
t65 = cos(qJ(6));
t109 = t10 * t61 - t24 * t65;
t108 = -t10 * t65 - t24 * t61;
t102 = r_i_i_C(3) + pkin(13);
t105 = r_i_i_C(1) * t65 - r_i_i_C(2) * t61 + pkin(5);
t62 = sin(qJ(4));
t104 = pkin(4) * t62 + pkin(10);
t66 = cos(qJ(4));
t54 = pkin(4) * t66 + pkin(3);
t103 = t102 * t55 + t54;
t96 = t55 * t58;
t95 = t56 * t58;
t94 = t58 * t59;
t91 = t63 * t64;
t90 = t63 * t67;
t88 = t64 * t94;
t86 = t59 * t97;
t84 = t98 * t64;
t83 = t98 * t67;
t81 = t89 * t58;
t80 = t58 * t86;
t79 = t104 * t58;
t78 = t89 * t97;
t69 = -pkin(12) - pkin(11);
t77 = r_i_i_C(1) * t61 + r_i_i_C(2) * t65 - t69;
t76 = t102 * t10 + t105 * (-t25 * t55 - t39 * t56);
t48 = -t64 * t78 + t67 * t68;
t73 = t64 * t68 + t67 * t78;
t29 = t48 * t98 + (-t60 * t73 + t80) * t63;
t70 = t58 * t73 + t60 * t86;
t13 = t29 * t55 - t56 * t70;
t14 = t29 * t56 + t55 * t70;
t75 = t102 * t14 - t105 * t13;
t38 = t63 * t81 + (t60 * t90 + t84) * t59;
t45 = t60 * t89 - t67 * t94;
t23 = t38 * t56 + t45 * t55;
t74 = t102 * t23 + t105 * (-t38 * t55 + t45 * t56);
t72 = t73 * t98;
t71 = -t105 * t56 - t103;
t44 = (-t60 * t91 + t83) * t59;
t43 = (t60 * t84 + t90) * t59;
t37 = -t98 * t81 + (-t60 * t83 + t91) * t59;
t35 = t44 * t56 + t55 * t88;
t33 = -t48 * t92 - t72;
t32 = t48 * t85 - t63 * t73;
t31 = -t46 * t98 - t47 * t92;
t30 = -t46 * t63 + t47 * t85;
t28 = t48 * t63 + t60 * t72 - t80 * t98;
t18 = t33 * t56 + t48 * t96;
t16 = t31 * t56 + t47 * t96;
t2 = t14 * t65 + t28 * t61;
t1 = -t14 * t61 + t28 * t65;
t3 = [t108 * r_i_i_C(1) + t109 * r_i_i_C(2) - t10 * pkin(5) + t24 * t69 - t47 * pkin(2) - t97 * pkin(1) + pkin(9) * t93 - t103 * t25 + (-t102 * t56 + t104) * t39 (t18 * t65 + t32 * t61) * r_i_i_C(1) + (-t18 * t61 + t32 * t65) * r_i_i_C(2) + t18 * pkin(5) + t33 * t54 - t32 * t69 - t73 * pkin(2) + t48 * t79 + t102 * (t33 * t55 - t48 * t95) t28 * t71 + t29 * t77 (-t29 * t62 + t66 * t70) * pkin(4) + t75, t75, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t68 * pkin(1) + t48 * pkin(2) + t14 * pkin(5) + pkin(9) * t86 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t102 * t13 + t104 * t70 - t28 * t69 + t29 * t54 (t16 * t65 + t30 * t61) * r_i_i_C(1) + (-t16 * t61 + t30 * t65) * r_i_i_C(2) + t16 * pkin(5) + t31 * t54 - t30 * t69 - t46 * pkin(2) + t47 * t79 + t102 * (t31 * t55 - t47 * t95) t24 * t71 + t25 * t77 (-t25 * t62 - t39 * t66) * pkin(4) + t76, t76, -r_i_i_C(1) * t109 + r_i_i_C(2) * t108; 0 (t35 * t65 + t43 * t61) * r_i_i_C(1) + (-t35 * t61 + t43 * t65) * r_i_i_C(2) + t35 * pkin(5) + t44 * t54 - t43 * t69 + t102 * (t44 * t55 - t56 * t88) + (t67 * pkin(2) + t64 * t79) * t59, t37 * t71 + t38 * t77 (-t38 * t62 + t45 * t66) * pkin(4) + t74, t74 (-t23 * t61 + t37 * t65) * r_i_i_C(1) + (-t23 * t65 - t37 * t61) * r_i_i_C(2);];
Ja_transl  = t3;
