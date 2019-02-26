% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR12_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR12_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobia_transl_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:07
% EndTime: 2019-02-26 21:21:08
% DurationCPUTime: 0.85s
% Computational Cost: add. (1425->125), mult. (4066->209), div. (0->0), fcn. (5408->18), ass. (0->86)
t107 = cos(qJ(4));
t106 = sin(qJ(1));
t50 = sin(pkin(14));
t61 = cos(qJ(1));
t97 = cos(pkin(14));
t99 = cos(pkin(6));
t88 = t99 * t97;
t47 = t106 * t50 - t61 * t88;
t52 = sin(pkin(6));
t96 = sin(pkin(7));
t91 = t52 * t96;
t98 = cos(pkin(7));
t113 = t98 * t47 + t61 * t91;
t93 = t50 * t99;
t48 = t106 * t97 + t61 * t93;
t57 = sin(qJ(3));
t60 = cos(qJ(3));
t110 = t113 * t60 + t48 * t57;
t51 = sin(pkin(8));
t92 = t52 * t98;
t82 = t47 * t96 - t61 * t92;
t114 = t82 * t51;
t37 = -t113 * t57 + t48 * t60;
t53 = cos(pkin(8));
t56 = sin(qJ(4));
t16 = t37 * t107 + (-t110 * t53 + t114) * t56;
t30 = t110 * t51 + t82 * t53;
t55 = sin(qJ(5));
t59 = cos(qJ(5));
t4 = t16 * t59 + t30 * t55;
t54 = sin(qJ(6));
t122 = t4 * t54;
t58 = cos(qJ(6));
t121 = t4 * t58;
t120 = -t16 * t55 + t30 * t59;
t116 = t107 * t114 - t37 * t56;
t109 = r_i_i_C(3) + pkin(13);
t80 = -t106 * t88 - t61 * t50;
t70 = t106 * t91 + t80 * t98;
t79 = -t106 * t93 + t61 * t97;
t65 = t79 * t57 - t70 * t60;
t69 = t106 * t92 - t80 * t96;
t62 = t65 * t51 + t69 * t53;
t108 = pkin(11) * t51;
t103 = t51 * t55;
t102 = t51 * t59;
t101 = t53 * t56;
t100 = t61 * t52;
t95 = t53 * t107;
t94 = t106 * t52;
t87 = t99 * t96;
t86 = t98 * t97;
t85 = t58 * r_i_i_C(1) - t54 * r_i_i_C(2) + pkin(5);
t84 = t54 * r_i_i_C(1) + t58 * r_i_i_C(2) + pkin(12);
t77 = -t109 * t55 - t85 * t59 - pkin(4);
t76 = -t97 * t91 + t99 * t98;
t73 = t76 * t51;
t72 = t110 * t107;
t71 = t60 * t87 + (-t50 * t57 + t60 * t86) * t52;
t68 = t71 * t107;
t67 = t69 * t51;
t63 = t65 * t107;
t44 = t57 * t87 + (t50 * t60 + t57 * t86) * t52;
t40 = t70 * t57 + t79 * t60;
t36 = -t71 * t51 + t76 * t53;
t33 = -t44 * t101 + t68;
t32 = t44 * t95 + t71 * t56;
t28 = t44 * t107 + (t71 * t53 + t73) * t56;
t27 = -t107 * t73 + t44 * t56 - t53 * t68;
t26 = t44 * t103 + t33 * t59;
t24 = -t40 * t101 - t63;
t23 = t40 * t95 - t65 * t56;
t22 = -t37 * t101 - t72;
t21 = -t110 * t56 + t37 * t95;
t20 = t40 * t107 + (-t65 * t53 + t67) * t56;
t19 = -t107 * t67 + t40 * t56 + t53 * t63;
t17 = -t110 * t95 + t116;
t15 = t53 * t72 - t116;
t14 = t28 * t59 + t36 * t55;
t12 = t40 * t103 + t24 * t59;
t10 = t37 * t103 + t22 * t59;
t8 = t20 * t59 + t62 * t55;
t7 = t20 * t55 - t62 * t59;
t2 = t19 * t54 + t8 * t58;
t1 = t19 * t58 - t8 * t54;
t3 = [(t17 * t54 - t121) * r_i_i_C(1) + (t17 * t58 + t122) * r_i_i_C(2) - t4 * pkin(5) - t16 * pkin(4) + t17 * pkin(12) - t37 * pkin(3) - t48 * pkin(2) - t106 * pkin(1) + qJ(2) * t100 + t109 * t120 - t30 * pkin(11) - t82 * pkin(10), t94 (t12 * t58 + t23 * t54) * r_i_i_C(1) + (-t12 * t54 + t23 * t58) * r_i_i_C(2) + t12 * pkin(5) + t24 * pkin(4) + t23 * pkin(12) - t65 * pkin(3) + t40 * t108 + t109 * (-t40 * t102 + t24 * t55) t77 * t19 + t84 * t20, t109 * t8 - t85 * t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t61 * pkin(1) + t79 * pkin(2) + t40 * pkin(3) + t20 * pkin(4) + t8 * pkin(5) + t69 * pkin(10) + t62 * pkin(11) + t19 * pkin(12) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t94 + t109 * t7, -t100 (t10 * t58 + t21 * t54) * r_i_i_C(1) + (-t10 * t54 + t21 * t58) * r_i_i_C(2) + t10 * pkin(5) + t22 * pkin(4) + t21 * pkin(12) - t110 * pkin(3) + t37 * t108 + t109 * (-t37 * t102 + t22 * t55) t77 * t15 + t84 * t16, t109 * t4 + t85 * t120 (t15 * t58 - t122) * r_i_i_C(1) + (-t15 * t54 - t121) * r_i_i_C(2); 0, t99 (t26 * t58 + t32 * t54) * r_i_i_C(1) + (-t26 * t54 + t32 * t58) * r_i_i_C(2) + t26 * pkin(5) + t33 * pkin(4) + t32 * pkin(12) + t71 * pkin(3) + t44 * t108 + t109 * (-t44 * t102 + t33 * t55) t77 * t27 + t84 * t28, t109 * t14 + t85 * (-t28 * t55 + t36 * t59) (-t14 * t54 + t27 * t58) * r_i_i_C(1) + (-t14 * t58 - t27 * t54) * r_i_i_C(2);];
Ja_transl  = t3;
