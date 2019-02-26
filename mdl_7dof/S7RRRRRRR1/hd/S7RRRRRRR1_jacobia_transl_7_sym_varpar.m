% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 7 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Ja_transl [3x7]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S7RRRRRRR1_jacobia_transl_7_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobia_transl_7_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobia_transl_7_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobia_transl_7_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_7_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:30
% EndTime: 2019-02-26 22:54:30
% DurationCPUTime: 0.89s
% Computational Cost: add. (562->129), mult. (1558->252), div. (0->0), fcn. (2050->14), ass. (0->95)
t78 = sin(qJ(2));
t79 = sin(qJ(1));
t103 = t78 * t79;
t84 = cos(qJ(3));
t77 = sin(qJ(3));
t86 = cos(qJ(1));
t94 = t86 * t77;
t85 = cos(qJ(2));
t99 = t79 * t85;
t64 = t84 * t99 + t94;
t76 = sin(qJ(4));
t83 = cos(qJ(4));
t46 = t103 * t76 + t64 * t83;
t93 = t86 * t84;
t63 = t77 * t99 - t93;
t75 = sin(qJ(5));
t82 = cos(qJ(5));
t23 = t46 * t75 + t63 * t82;
t24 = t46 * t82 - t63 * t75;
t102 = t78 * t83;
t45 = -t79 * t102 + t64 * t76;
t74 = sin(qJ(6));
t81 = cos(qJ(6));
t4 = t24 * t81 + t45 * t74;
t73 = sin(qJ(7));
t80 = cos(qJ(7));
t123 = t23 * t80 + t4 * t73;
t122 = t23 * t73 - t4 * t80;
t119 = t24 * t74 - t45 * t81;
t116 = pkin(4) + r_i_i_C(3);
t115 = pkin(3) * t76;
t112 = t73 * t75;
t111 = t74 * t76;
t110 = t74 * t82;
t109 = t75 * t80;
t108 = t75 * t83;
t107 = t76 * t81;
t106 = t77 * t78;
t105 = t77 * t85;
t104 = t78 * t76;
t101 = t78 * t84;
t100 = t78 * t86;
t98 = t81 * t82;
t97 = t82 * t83;
t96 = t85 * t76;
t95 = t85 * t83;
t92 = t77 * t104;
t91 = t75 * t106;
t90 = t82 * t106;
t89 = t80 * r_i_i_C(1) - t73 * r_i_i_C(2);
t88 = -t73 * r_i_i_C(1) - t80 * r_i_i_C(2);
t62 = t101 * t83 - t96;
t61 = t101 * t76 + t95;
t87 = t116 * t74 - t89 * t81;
t68 = -t79 * t77 + t85 * t93;
t67 = -t79 * t84 - t85 * t94;
t66 = t84 * t95 + t104;
t65 = t84 * t96 - t102;
t58 = t62 * t86;
t57 = t61 * t86;
t56 = t62 * t79;
t55 = t61 * t79;
t54 = (-t75 * t84 - t77 * t97) * t78;
t53 = (-t108 * t77 + t82 * t84) * t78;
t52 = t100 * t76 + t68 * t83;
t51 = -t100 * t83 + t68 * t76;
t50 = -t105 * t75 + t66 * t82;
t49 = t105 * t82 + t66 * t75;
t44 = t62 * t82 - t91;
t43 = t62 * t75 + t90;
t42 = -t58 * t82 + t86 * t91;
t41 = -t58 * t75 - t86 * t90;
t40 = -t56 * t82 + t79 * t91;
t39 = -t56 * t75 - t79 * t90;
t38 = t54 * t81 - t74 * t92;
t36 = t67 * t97 - t68 * t75;
t35 = t108 * t67 + t68 * t82;
t34 = -t63 * t97 - t64 * t75;
t33 = -t108 * t63 + t64 * t82;
t32 = -t61 * t98 + t62 * t74;
t30 = t52 * t82 + t67 * t75;
t29 = t52 * t75 - t67 * t82;
t28 = t50 * t81 + t65 * t74;
t22 = t44 * t81 + t61 * t74;
t20 = t42 * t81 - t57 * t74;
t18 = t40 * t81 - t55 * t74;
t16 = t111 * t67 + t36 * t81;
t14 = -t111 * t63 + t34 * t81;
t12 = -t51 * t98 + t52 * t74;
t10 = -t45 * t98 + t46 * t74;
t8 = t30 * t81 + t51 * t74;
t7 = -t30 * t74 + t51 * t81;
t2 = -t29 * t73 + t8 * t80;
t1 = -t29 * t80 - t8 * t73;
t3 = [pkin(2) * t103 - t45 * pkin(3) + t122 * r_i_i_C(1) + t123 * r_i_i_C(2) + t116 * t119 (t20 * t80 - t41 * t73) * r_i_i_C(1) + (-t20 * t73 - t41 * t80) * r_i_i_C(2) - t57 * pkin(3) - t86 * t85 * pkin(2) + t116 * (-t42 * t74 - t57 * t81) (t16 * t80 - t35 * t73) * r_i_i_C(1) + (-t16 * t73 - t35 * t80) * r_i_i_C(2) + t67 * t115 + t116 * (t107 * t67 - t36 * t74) (t112 * t51 + t12 * t80) * r_i_i_C(1) + (t109 * t51 - t12 * t73) * r_i_i_C(2) + t52 * pkin(3) + t116 * (t110 * t51 + t52 * t81) t29 * t87 + t30 * t88, -t116 * t8 + t89 * t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -pkin(2) * t100 + t51 * pkin(3) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t116 * t7 (t18 * t80 - t39 * t73) * r_i_i_C(1) + (-t18 * t73 - t39 * t80) * r_i_i_C(2) - t55 * pkin(3) - pkin(2) * t99 + t116 * (-t40 * t74 - t55 * t81) (t14 * t80 - t33 * t73) * r_i_i_C(1) + (-t14 * t73 - t33 * t80) * r_i_i_C(2) - t63 * t115 + t116 * (-t107 * t63 - t34 * t74) (t10 * t80 + t112 * t45) * r_i_i_C(1) + (-t10 * t73 + t109 * t45) * r_i_i_C(2) + t46 * pkin(3) + t116 * (t110 * t45 + t46 * t81) t23 * t87 + t24 * t88, -t116 * t4 - t89 * t119, -t123 * r_i_i_C(1) + t122 * r_i_i_C(2); 0 (t28 * t80 - t49 * t73) * r_i_i_C(1) + (-t28 * t73 - t49 * t80) * r_i_i_C(2) + t65 * pkin(3) - t78 * pkin(2) + t116 * (-t50 * t74 + t65 * t81) (t38 * t80 - t53 * t73) * r_i_i_C(1) + (-t38 * t73 - t53 * t80) * r_i_i_C(2) - pkin(3) * t92 + t116 * (-t54 * t74 - t81 * t92) (t112 * t61 + t32 * t80) * r_i_i_C(1) + (t109 * t61 - t32 * t73) * r_i_i_C(2) + t62 * pkin(3) + t116 * (t110 * t61 + t62 * t81) t43 * t87 + t44 * t88, -t116 * t22 + t89 * (-t44 * t74 + t61 * t81) (-t22 * t73 - t43 * t80) * r_i_i_C(1) + (-t22 * t80 + t43 * t73) * r_i_i_C(2);];
Ja_transl  = t3;
