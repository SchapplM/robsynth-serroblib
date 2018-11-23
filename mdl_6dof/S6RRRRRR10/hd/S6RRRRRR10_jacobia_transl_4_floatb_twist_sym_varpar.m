% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR10_jacobia_transl_4_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_transl_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobia_transl_4_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_transl_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:16
% EndTime: 2018-11-23 11:27:17
% DurationCPUTime: 0.41s
% Computational Cost: add. (1331->116), mult. (1341->172), div. (0->0), fcn. (1316->26), ass. (0->79)
t60 = sin(qJ(2));
t61 = sin(qJ(1));
t65 = cos(qJ(1));
t92 = pkin(6) + qJ(2);
t78 = cos(t92) / 0.2e1;
t93 = pkin(6) - qJ(2);
t86 = cos(t93);
t66 = t86 / 0.2e1 + t78;
t23 = t61 * t60 - t65 * t66;
t53 = sin(pkin(7));
t56 = cos(pkin(7));
t54 = sin(pkin(6));
t94 = t54 * t65;
t19 = -t23 * t53 + t56 * t94;
t75 = sin(t92) / 0.2e1;
t83 = sin(t93);
t35 = t75 - t83 / 0.2e1;
t64 = cos(qJ(2));
t24 = t65 * t35 + t61 * t64;
t90 = pkin(7) + qJ(3);
t74 = sin(t90) / 0.2e1;
t91 = pkin(7) - qJ(3);
t82 = sin(t91);
t33 = t74 - t82 / 0.2e1;
t77 = cos(t90) / 0.2e1;
t85 = cos(t91);
t38 = t77 - t85 / 0.2e1;
t63 = cos(qJ(3));
t3 = -t23 * t33 + t24 * t63 + t38 * t94;
t88 = pkin(8) + qJ(4);
t73 = sin(t88) / 0.2e1;
t89 = pkin(8) - qJ(4);
t81 = sin(t89);
t30 = t73 + t81 / 0.2e1;
t76 = cos(t88) / 0.2e1;
t84 = cos(t89);
t37 = t84 / 0.2e1 + t76;
t32 = t74 + t82 / 0.2e1;
t39 = t85 / 0.2e1 + t77;
t59 = sin(qJ(3));
t5 = t23 * t39 + t24 * t59 + t32 * t94;
t58 = sin(qJ(4));
t103 = t19 * t30 + t3 * t58 + t5 * t37;
t31 = t73 - t81 / 0.2e1;
t36 = t76 - t84 / 0.2e1;
t62 = cos(qJ(4));
t102 = -t19 * t36 - t3 * t62 + t5 * t31;
t101 = pkin(12) + r_i_i_C(3);
t100 = t53 * pkin(11);
t99 = t30 * t53;
t98 = t36 * t53;
t40 = t78 - t86 / 0.2e1;
t97 = t40 * t53;
t55 = cos(pkin(8));
t96 = t53 * t55;
t95 = t54 * t61;
t52 = sin(pkin(8));
t87 = t101 * t52;
t28 = t61 * t35 - t65 * t64;
t72 = t62 * r_i_i_C(1) - t58 * r_i_i_C(2) + pkin(3);
t26 = -t65 * t60 - t61 * t66;
t21 = -t26 * t53 + t56 * t95;
t7 = t26 * t39 + t28 * t59 + t32 * t95;
t9 = -t26 * t33 + t28 * t63 + t38 * t95;
t68 = t21 * t36 - t7 * t31 + t62 * t9;
t34 = t75 + t83 / 0.2e1;
t57 = cos(pkin(6));
t15 = t34 * t33 - t57 * t38 - t40 * t63;
t67 = t31 * r_i_i_C(1) + t37 * r_i_i_C(2) - t87;
t22 = -t34 * t53 + t57 * t56;
t18 = t40 * t33 + t34 * t63;
t17 = -t34 * t59 + t40 * t39;
t14 = t57 * t32 + t34 * t39 + t40 * t59;
t13 = t26 * t63 + t28 * t33;
t12 = -t26 * t59 + t28 * t39;
t11 = -t23 * t63 - t24 * t33;
t10 = t23 * t59 - t24 * t39;
t1 = t21 * t30 + t7 * t37 + t58 * t9;
t2 = [t102 * r_i_i_C(1) + t103 * r_i_i_C(2) - t3 * pkin(3) - t24 * pkin(2) - t61 * pkin(1) + pkin(10) * t94 + t19 * pkin(11) + t101 * (t19 * t55 - t5 * t52) (t12 * t31 + t13 * t62 + t28 * t98) * r_i_i_C(1) + (t12 * t37 - t13 * t58 - t28 * t99) * r_i_i_C(2) + t13 * pkin(3) + t26 * pkin(2) - t28 * t100 + t101 * (-t12 * t52 - t28 * t96) (t9 * t31 + t7 * t62) * r_i_i_C(1) + (t9 * t37 - t7 * t58) * r_i_i_C(2) + t7 * pkin(3) - t9 * t87, t1 * r_i_i_C(1) + t68 * r_i_i_C(2), 0, 0; t65 * pkin(1) - t28 * pkin(2) - t9 * pkin(3) + pkin(10) * t95 + t21 * pkin(11) - t68 * r_i_i_C(1) + t1 * r_i_i_C(2) + t101 * (t21 * t55 - t7 * t52) (t10 * t31 + t11 * t62 - t24 * t98) * r_i_i_C(1) + (t10 * t37 - t11 * t58 + t24 * t99) * r_i_i_C(2) + t11 * pkin(3) - t23 * pkin(2) + t24 * t100 + t101 * (-t10 * t52 + t24 * t96) -t67 * t3 - t5 * t72, -t103 * r_i_i_C(1) + t102 * r_i_i_C(2), 0, 0; 0 (t17 * t31 + t18 * t62 + t36 * t97) * r_i_i_C(1) + (t17 * t37 - t18 * t58 - t30 * t97) * r_i_i_C(2) + t18 * pkin(3) + t34 * pkin(2) - pkin(11) * t97 + t101 * (-t17 * t52 - t40 * t96) t72 * t14 - t67 * t15 (t14 * t37 - t15 * t58 + t22 * t30) * r_i_i_C(1) + (-t14 * t31 - t15 * t62 + t22 * t36) * r_i_i_C(2), 0, 0;];
Ja_transl  = t2;
