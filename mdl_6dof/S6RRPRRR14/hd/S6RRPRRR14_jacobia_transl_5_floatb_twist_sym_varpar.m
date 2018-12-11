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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14_jacobia_transl_5_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_transl_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobia_transl_5_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_transl_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:20
% EndTime: 2018-12-10 18:38:21
% DurationCPUTime: 0.64s
% Computational Cost: add. (2472->134), mult. (2506->199), div. (0->0), fcn. (2507->28), ass. (0->93)
t81 = sin(pkin(6));
t93 = cos(qJ(1));
t109 = t81 * t93;
t88 = sin(qJ(2));
t89 = sin(qJ(1));
t107 = pkin(6) - qJ(2);
t103 = cos(t107);
t106 = pkin(6) + qJ(2);
t99 = cos(t106) / 0.2e1;
t94 = t103 / 0.2e1 + t99;
t47 = t89 * t88 - t93 * t94;
t101 = sin(t107);
t97 = sin(t106) / 0.2e1;
t61 = t97 - t101 / 0.2e1;
t92 = cos(qJ(2));
t48 = t93 * t61 + t89 * t92;
t76 = pkin(7) + pkin(14);
t65 = sin(t76) / 0.2e1;
t77 = pkin(7) - pkin(14);
t74 = sin(t77);
t54 = t65 + t74 / 0.2e1;
t66 = cos(t77) / 0.2e1;
t75 = cos(t76);
t56 = t66 + t75 / 0.2e1;
t78 = sin(pkin(14));
t29 = t54 * t109 + t47 * t56 + t48 * t78;
t80 = sin(pkin(7));
t84 = cos(pkin(7));
t43 = t84 * t109 - t47 * t80;
t79 = sin(pkin(8));
t83 = cos(pkin(8));
t20 = t29 * t79 - t43 * t83;
t55 = t65 - t74 / 0.2e1;
t57 = t66 - t75 / 0.2e1;
t82 = cos(pkin(14));
t30 = t57 * t109 + t47 * t55 - t48 * t82;
t105 = pkin(8) - qJ(4);
t100 = sin(t105);
t104 = pkin(8) + qJ(4);
t96 = sin(t104) / 0.2e1;
t59 = t96 - t100 / 0.2e1;
t102 = cos(t105);
t98 = cos(t104) / 0.2e1;
t62 = t98 - t102 / 0.2e1;
t91 = cos(qJ(4));
t7 = t29 * t59 + t30 * t91 - t43 * t62;
t86 = sin(qJ(5));
t90 = cos(qJ(5));
t118 = t20 * t90 + t7 * t86;
t117 = -t20 * t86 + t7 * t90;
t58 = t96 + t100 / 0.2e1;
t63 = t102 / 0.2e1 + t98;
t87 = sin(qJ(4));
t116 = -t29 * t63 + t30 * t87 - t43 * t58;
t115 = r_i_i_C(3) + pkin(12);
t114 = t58 * t80;
t113 = t62 * t80;
t64 = t99 - t103 / 0.2e1;
t112 = t64 * t80;
t111 = t80 * t83;
t110 = t81 * t89;
t108 = t80 * qJ(3);
t50 = -t93 * t88 - t89 * t94;
t52 = t89 * t61 - t93 * t92;
t31 = t54 * t110 + t50 * t56 + t52 * t78;
t45 = t84 * t110 - t50 * t80;
t22 = -t31 * t79 + t45 * t83;
t95 = t90 * r_i_i_C(1) - t86 * r_i_i_C(2) + pkin(4);
t34 = t47 * t78 - t48 * t56;
t23 = t111 * t48 - t34 * t79;
t36 = -t50 * t78 + t52 * t56;
t24 = -t52 * t111 - t36 * t79;
t60 = t97 + t101 / 0.2e1;
t41 = t64 * t56 - t60 * t78;
t33 = -t64 * t111 - t41 * t79;
t32 = t57 * t110 + t50 * t55 - t52 * t82;
t9 = t31 * t59 + t32 * t91 - t45 * t62;
t85 = cos(pkin(6));
t38 = t85 * t54 + t60 * t56 + t64 * t78;
t39 = t60 * t55 + t85 * t57 - t64 * t82;
t46 = -t60 * t80 + t85 * t84;
t16 = t38 * t59 + t39 * t91 - t46 * t62;
t42 = t64 * t55 + t60 * t82;
t37 = t50 * t82 + t52 * t55;
t35 = -t47 * t82 - t48 * t55;
t26 = -t38 * t79 + t46 * t83;
t19 = t62 * t112 + t41 * t59 + t42 * t91;
t14 = t52 * t113 + t36 * t59 + t37 * t91;
t12 = -t113 * t48 + t34 * t59 + t35 * t91;
t8 = -t31 * t63 + t32 * t87 - t45 * t58;
t2 = t22 * t86 + t9 * t90;
t1 = t22 * t90 - t9 * t86;
t3 = [-t89 * pkin(1) - t48 * pkin(2) + t30 * pkin(3) + t7 * pkin(4) + pkin(10) * t109 - t20 * pkin(11) + t117 * r_i_i_C(1) - t118 * r_i_i_C(2) + t43 * qJ(3) + t115 * t116 (t14 * t90 + t24 * t86) * r_i_i_C(1) + (-t14 * t86 + t24 * t90) * r_i_i_C(2) + t14 * pkin(4) + t37 * pkin(3) + t50 * pkin(2) - t52 * t108 + t115 * (t114 * t52 - t36 * t63 + t37 * t87) + t24 * pkin(11), t45, t115 * t9 - t95 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t93 * pkin(1) - pkin(2) * t52 + t32 * pkin(3) + t9 * pkin(4) + pkin(10) * t110 + t22 * pkin(11) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t45 * qJ(3) + t115 * t8 (t12 * t90 + t23 * t86) * r_i_i_C(1) + (-t12 * t86 + t23 * t90) * r_i_i_C(2) + t12 * pkin(4) + t35 * pkin(3) - t47 * pkin(2) + t48 * t108 + t115 * (-t114 * t48 - t34 * t63 + t35 * t87) + t23 * pkin(11), -t43, -t115 * t7 + t95 * t116, t118 * r_i_i_C(1) + t117 * r_i_i_C(2), 0; 0 (t19 * t90 + t33 * t86) * r_i_i_C(1) + (-t19 * t86 + t33 * t90) * r_i_i_C(2) + t19 * pkin(4) + t42 * pkin(3) + t60 * pkin(2) - t64 * t108 + t115 * (t112 * t58 - t41 * t63 + t42 * t87) + t33 * pkin(11), t46, t115 * t16 + t95 * (t38 * t63 - t39 * t87 + t46 * t58) (-t16 * t86 + t26 * t90) * r_i_i_C(1) + (-t16 * t90 - t26 * t86) * r_i_i_C(2), 0;];
Ja_transl  = t3;
