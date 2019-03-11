% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t65 = sin(qJ(3));
t102 = t65 * pkin(2);
t51 = pkin(9) + t102;
t64 = sin(qJ(4));
t62 = t64 ^ 2;
t67 = cos(qJ(4));
t63 = t67 ^ 2;
t87 = t62 + t63;
t89 = t87 * t51;
t55 = t64 * qJ(5);
t111 = -t67 * pkin(4) - t55;
t68 = pkin(4) + pkin(5);
t86 = t67 * qJ(5);
t38 = t68 * t64 - t86;
t100 = cos(qJ(2));
t66 = sin(qJ(2));
t99 = cos(qJ(3));
t36 = t65 * t100 + t99 * t66;
t110 = 0.2e1 * t36;
t53 = -t100 * pkin(2) - pkin(1);
t109 = 0.2e1 * t53;
t108 = -0.2e1 * t64;
t107 = 0.2e1 * t64;
t106 = -0.2e1 * t67;
t105 = 0.2e1 * t67;
t35 = -t99 * t100 + t65 * t66;
t104 = pkin(9) * t35;
t33 = t35 * pkin(4);
t103 = t64 * pkin(9);
t61 = t99 * pkin(2);
t52 = -t61 - pkin(3);
t101 = pkin(3) - t52;
t43 = (-pkin(8) - pkin(7)) * t66;
t80 = t100 * pkin(7);
t44 = t100 * pkin(8) + t80;
t19 = -t99 * t43 + t65 * t44;
t76 = pkin(4) * t64 - t86;
t11 = t36 * t76 + t19;
t98 = t11 * t64;
t97 = t11 * t67;
t96 = t19 * t67;
t95 = t35 * t51;
t94 = t64 * t36;
t93 = t64 * t51;
t92 = t64 * t67;
t26 = t67 * t36;
t14 = t35 * pkin(3) - t36 * pkin(9) + t53;
t20 = t65 * t43 + t99 * t44;
t9 = t67 * t14 - t64 * t20;
t10 = t64 * t14 + t67 * t20;
t82 = pkin(3) - t111;
t28 = -t61 - t82;
t59 = t67 * pkin(5);
t21 = t59 - t28;
t29 = t59 + t82;
t91 = t21 + t29;
t90 = -t28 + t82;
t88 = t87 * pkin(9);
t54 = t64 * qJ(6);
t85 = t67 * qJ(6);
t84 = 0.2e1 * t100;
t83 = -0.2e1 * t36 * t35;
t5 = -t33 - t9;
t30 = t35 * qJ(5);
t4 = t30 + t10;
t81 = 0.2e1 * t30 + t10;
t79 = -pkin(3) * t36 - t104;
t73 = t36 * t85 - t5;
t2 = -t35 * pkin(5) - t73;
t23 = t36 * t54;
t3 = t23 + t4;
t78 = -t2 * t64 - t3 * t67;
t1 = t4 * t67 + t5 * t64;
t77 = t36 * t82 + t104;
t75 = t28 * t36 - t95;
t74 = t36 * t52 - t95;
t71 = qJ(5) ^ 2;
t70 = 0.2e1 * qJ(5);
t58 = t67 * pkin(9);
t48 = 0.2e1 * t92;
t47 = t67 * t51;
t42 = t58 - t85;
t40 = -t54 + t103;
t34 = t36 ^ 2;
t32 = t47 - t85;
t31 = -t54 + t93;
t25 = t67 * t35;
t24 = t64 * t35;
t22 = t64 * t26;
t18 = t19 * t64;
t15 = (-t62 + t63) * t36;
t8 = -t38 * t36 - t19;
t7 = t8 * t67;
t6 = t8 * t64;
t12 = [1, 0, 0, t66 ^ 2, t66 * t84, 0, 0, 0, pkin(1) * t84, -0.2e1 * pkin(1) * t66, t34, t83, 0, 0, 0, t35 * t109, t36 * t109, t63 * t34, -0.2e1 * t34 * t92, 0.2e1 * t35 * t26, t64 * t83, t35 ^ 2, 0.2e1 * t19 * t94 + 0.2e1 * t9 * t35, -0.2e1 * t10 * t35 + 0.2e1 * t19 * t26, 0.2e1 * t11 * t94 - 0.2e1 * t5 * t35 (-t4 * t64 + t5 * t67) * t110, -0.2e1 * t11 * t26 + 0.2e1 * t4 * t35, t11 ^ 2 + t4 ^ 2 + t5 ^ 2, -0.2e1 * t2 * t35 - 0.2e1 * t8 * t94, 0.2e1 * t8 * t26 + 0.2e1 * t3 * t35 (-t2 * t67 + t3 * t64) * t110, t2 ^ 2 + t3 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, t66, t100, 0, -t66 * pkin(7), -t80, 0, 0, t36, -t35, 0, -t19, -t20, t22, t15, t24, t25, 0, t64 * t74 - t96, t67 * t74 + t18, t64 * t75 - t97, t1, -t67 * t75 - t98, t1 * t51 + t11 * t28, -t21 * t94 - t31 * t35 + t7, t21 * t26 + t32 * t35 + t6 (-t31 * t67 + t32 * t64) * t36 + t78, t2 * t31 + t8 * t21 + t3 * t32; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t61, -0.2e1 * t102, t62, t48, 0, 0, 0, t52 * t106, t52 * t107, t28 * t106, 0.2e1 * t89, t28 * t108, t87 * t51 ^ 2 + t28 ^ 2, t21 * t105, t21 * t107, -0.2e1 * t31 * t64 - 0.2e1 * t32 * t67, t21 ^ 2 + t31 ^ 2 + t32 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, 0, -t19, -t20, t22, t15, t24, t25, 0, t64 * t79 - t96, t67 * t79 + t18, -t64 * t77 - t97, t1, t67 * t77 - t98, pkin(9) * t1 - t11 * t82, -t29 * t94 - t40 * t35 + t7, t29 * t26 + t42 * t35 + t6 (-t40 * t67 + t42 * t64) * t36 + t78, t2 * t40 + t8 * t29 + t3 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t61, -t102, t62, t48, 0, 0, 0, t101 * t67, -t101 * t64, t90 * t67, t88 + t89, t90 * t64, pkin(9) * t89 - t28 * t82, t91 * t67, t91 * t64 (-t32 - t42) * t67 + (-t31 - t40) * t64, t21 * t29 + t31 * t40 + t32 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t62, t48, 0, 0, 0, pkin(3) * t105, pkin(3) * t108, -t82 * t106, 0.2e1 * t88, -t82 * t108, t87 * pkin(9) ^ 2 + t82 ^ 2, t29 * t105, t29 * t107, -0.2e1 * t40 * t64 - 0.2e1 * t42 * t67, t29 ^ 2 + t40 ^ 2 + t42 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t94, t35, t9, -t10, -t5 + t33, t111 * t36, t81, -t5 * pkin(4) + t4 * qJ(5) (pkin(5) + t68) * t35 + t73, t23 + t81 (t67 * t68 + t55) * t36, t3 * qJ(5) - t2 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t67, 0, -t93, -t47, -t93, -t76, t47, -t76 * t51, -t31, t32, t38, t32 * qJ(5) - t31 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t67, 0, -t103, -t58, -t103, -t76, t58, -t76 * pkin(9), -t40, t42, t38, t42 * qJ(5) - t40 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, t70, pkin(4) ^ 2 + t71, 0.2e1 * t68, t70, 0, t68 ^ 2 + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t26, 0, t5, -t35, 0, -t26, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, t93, 0, 0, -t64, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, t103, 0, 0, -t64, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), -1, 0, 0, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, t26, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t64, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t64, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t12;
