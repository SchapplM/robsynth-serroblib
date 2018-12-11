% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14_jacobia_transl_4_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_transl_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobia_transl_4_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_transl_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:20
% EndTime: 2018-12-10 18:38:21
% DurationCPUTime: 0.37s
% Computational Cost: add. (987->106), mult. (1013->155), div. (0->0), fcn. (1007->26), ass. (0->76)
t60 = sin(qJ(2));
t61 = sin(qJ(1));
t64 = cos(qJ(1));
t82 = pkin(6) + qJ(2);
t73 = cos(t82) / 0.2e1;
t83 = pkin(6) - qJ(2);
t79 = cos(t83);
t65 = t79 / 0.2e1 + t73;
t20 = t61 * t60 - t64 * t65;
t53 = sin(pkin(7));
t57 = cos(pkin(7));
t54 = sin(pkin(6));
t85 = t54 * t64;
t16 = -t20 * t53 + t57 * t85;
t80 = pkin(8) + qJ(4);
t70 = sin(t80) / 0.2e1;
t81 = pkin(8) - qJ(4);
t76 = sin(t81);
t31 = t70 + t76 / 0.2e1;
t72 = cos(t80) / 0.2e1;
t78 = cos(t81);
t36 = t78 / 0.2e1 + t72;
t71 = sin(t82) / 0.2e1;
t77 = sin(t83);
t34 = t71 - t77 / 0.2e1;
t63 = cos(qJ(2));
t21 = t64 * t34 + t61 * t63;
t49 = pkin(7) + pkin(14);
t38 = sin(t49) / 0.2e1;
t50 = pkin(7) - pkin(14);
t47 = sin(t50);
t27 = t38 + t47 / 0.2e1;
t39 = cos(t50) / 0.2e1;
t48 = cos(t49);
t29 = t39 + t48 / 0.2e1;
t51 = sin(pkin(14));
t4 = t20 * t29 + t21 * t51 + t27 * t85;
t28 = t38 - t47 / 0.2e1;
t30 = t39 - t48 / 0.2e1;
t55 = cos(pkin(14));
t5 = t20 * t28 - t21 * t55 + t30 * t85;
t59 = sin(qJ(4));
t93 = t16 * t31 + t4 * t36 - t5 * t59;
t32 = t70 - t76 / 0.2e1;
t35 = t72 - t78 / 0.2e1;
t62 = cos(qJ(4));
t92 = -t16 * t35 + t4 * t32 + t5 * t62;
t91 = r_i_i_C(3) + pkin(11);
t90 = t31 * t53;
t89 = t35 * t53;
t37 = t73 - t79 / 0.2e1;
t88 = t37 * t53;
t56 = cos(pkin(8));
t87 = t53 * t56;
t86 = t54 * t61;
t84 = t53 * qJ(3);
t25 = t61 * t34 - t64 * t63;
t23 = -t64 * t60 - t61 * t65;
t18 = -t23 * t53 + t57 * t86;
t6 = t23 * t29 + t25 * t51 + t27 * t86;
t7 = t23 * t28 - t25 * t55 + t30 * t86;
t66 = t18 * t35 - t6 * t32 - t7 * t62;
t58 = cos(pkin(6));
t52 = sin(pkin(8));
t33 = t71 + t77 / 0.2e1;
t19 = -t33 * t53 + t58 * t57;
t15 = t37 * t28 + t33 * t55;
t14 = t37 * t29 - t33 * t51;
t13 = t33 * t28 + t58 * t30 - t37 * t55;
t12 = t58 * t27 + t33 * t29 + t37 * t51;
t11 = t23 * t55 + t25 * t28;
t10 = -t23 * t51 + t25 * t29;
t9 = -t20 * t55 - t21 * t28;
t8 = t20 * t51 - t21 * t29;
t1 = t18 * t31 + t6 * t36 - t7 * t59;
t2 = [t92 * r_i_i_C(1) + t93 * r_i_i_C(2) + t5 * pkin(3) - t21 * pkin(2) - t61 * pkin(1) + pkin(10) * t85 + t16 * qJ(3) + t91 * (t16 * t56 - t4 * t52) (t10 * t32 + t11 * t62 + t25 * t89) * r_i_i_C(1) + (t10 * t36 - t11 * t59 - t25 * t90) * r_i_i_C(2) + t11 * pkin(3) + t23 * pkin(2) - t25 * t84 + t91 * (-t10 * t52 - t25 * t87) t18, t1 * r_i_i_C(1) + t66 * r_i_i_C(2), 0, 0; t64 * pkin(1) - t25 * pkin(2) + t7 * pkin(3) + pkin(10) * t86 - t66 * r_i_i_C(1) + t1 * r_i_i_C(2) + t18 * qJ(3) + t91 * (t18 * t56 - t6 * t52) (-t21 * t89 + t8 * t32 + t9 * t62) * r_i_i_C(1) + (t21 * t90 + t8 * t36 - t9 * t59) * r_i_i_C(2) + t9 * pkin(3) - t20 * pkin(2) + t21 * t84 + t91 * (t21 * t87 - t8 * t52) -t16, -t93 * r_i_i_C(1) + t92 * r_i_i_C(2), 0, 0; 0 (t14 * t32 + t15 * t62 + t35 * t88) * r_i_i_C(1) + (t14 * t36 - t15 * t59 - t31 * t88) * r_i_i_C(2) + t15 * pkin(3) + t33 * pkin(2) - t37 * t84 + t91 * (-t14 * t52 - t37 * t87) t19 (t12 * t36 - t13 * t59 + t19 * t31) * r_i_i_C(1) + (-t12 * t32 - t13 * t62 + t19 * t35) * r_i_i_C(2), 0, 0;];
Ja_transl  = t2;
