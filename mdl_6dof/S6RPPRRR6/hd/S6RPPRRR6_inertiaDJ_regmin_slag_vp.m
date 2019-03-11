% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPPRRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:31:29
% EndTime: 2019-03-09 02:31:32
% DurationCPUTime: 1.10s
% Computational Cost: add. (755->143), mult. (1816->268), div. (0->0), fcn. (1533->6), ass. (0->103)
t64 = cos(qJ(4));
t51 = t64 * qJD(4);
t57 = -pkin(7) + qJ(2);
t61 = sin(qJ(4));
t69 = -t61 * qJD(2) - t57 * t51;
t59 = sin(qJ(6));
t60 = sin(qJ(5));
t62 = cos(qJ(6));
t63 = cos(qJ(5));
t33 = t59 * t63 + t62 * t60;
t121 = t33 * t64;
t102 = t61 * qJD(4);
t88 = t63 * t102;
t106 = qJD(5) * t64;
t94 = t60 * t106;
t25 = -t88 - t94;
t54 = t61 ^ 2;
t56 = t64 ^ 2;
t81 = (t54 - t56) * qJD(4);
t55 = t63 ^ 2;
t110 = t60 ^ 2 - t55;
t82 = t110 * qJD(5);
t100 = qJD(5) + qJD(6);
t68 = -t64 * qJD(2) + t57 * t102;
t79 = pkin(4) * t64 + pkin(8) * t61;
t120 = t79 * qJD(5) + t68;
t119 = pkin(8) + pkin(9);
t90 = t60 * t102;
t92 = t63 * t106;
t27 = t90 - t92;
t31 = t79 * qJD(4) + qJD(3);
t58 = pkin(1) + qJ(3);
t78 = t61 * pkin(4) - t64 * pkin(8);
t36 = t78 + t58;
t52 = qJD(5) * t63;
t107 = qJD(5) * t60;
t91 = t57 * t107;
t6 = -t60 * t31 - t36 * t52 + t61 * t91 + t69 * t63;
t5 = t27 * pkin(9) - t6;
t118 = t62 * t5;
t17 = t100 * t33;
t117 = t17 * t61;
t114 = t60 * t64;
t28 = t60 * t36;
t113 = t61 * t63;
t44 = t57 * t113;
t15 = -pkin(9) * t114 + t28 + t44;
t116 = t59 * t15;
t115 = t59 * t60;
t112 = t62 * t15;
t111 = t63 * t64;
t108 = t54 + t56;
t105 = qJD(6) * t59;
t104 = qJD(6) * t62;
t101 = qJ(2) * qJD(2);
t99 = -0.2e1 * pkin(4) * qJD(5);
t98 = pkin(5) * t107;
t97 = pkin(5) * t51;
t96 = pkin(5) * t105;
t95 = pkin(5) * t104;
t93 = t61 * t52;
t89 = t60 * t52;
t87 = t61 * t51;
t73 = t63 * t31 - t57 * t93;
t4 = (pkin(5) * t64 + pkin(9) * t113) * qJD(4) + ((pkin(9) * t64 - t36) * qJD(5) + t69) * t60 + t73;
t85 = t62 * t4 - t59 * t5;
t29 = t63 * t36;
t14 = -pkin(9) * t111 + t29 + (-t57 * t60 + pkin(5)) * t61;
t84 = -t61 * pkin(5) - t14;
t83 = qJD(5) * t119;
t80 = t60 * t88;
t77 = t62 * t14 - t116;
t76 = t59 * t14 + t112;
t42 = t119 * t60;
t43 = t119 * t63;
t75 = -t62 * t42 - t59 * t43;
t74 = -t59 * t42 + t62 * t43;
t32 = -t62 * t63 + t115;
t16 = t100 * t115 - t63 * t104 - t62 * t52;
t72 = -t16 * t61 + t33 * t51;
t71 = t32 * t51 + t117;
t20 = t32 * t64;
t67 = qJD(4) * t32;
t66 = t78 * qJD(4) - t57 * t106;
t65 = 0.2e1 * qJD(2);
t50 = -t63 * pkin(5) - pkin(4);
t45 = 0.2e1 * t87;
t35 = t63 * t83;
t34 = t60 * t83;
t30 = (pkin(5) * t60 - t57) * t64;
t26 = t60 * t51 + t93;
t24 = t61 * t107 - t63 * t51;
t18 = -t27 * pkin(5) + t68;
t13 = -t74 * qJD(6) + t59 * t34 - t62 * t35;
t12 = -t75 * qJD(6) + t62 * t34 + t59 * t35;
t11 = -t105 * t114 + (t100 * t111 - t90) * t62 + t25 * t59;
t10 = t100 * t61 * t32 - qJD(4) * t121;
t9 = -t100 * t121 + t61 * t67;
t8 = t64 * t67 + t117;
t7 = (-qJD(5) * t36 + t69) * t60 + t73;
t2 = -t76 * qJD(6) + t85;
t1 = -t77 * qJD(6) - t59 * t4 - t118;
t3 = [0, 0, 0, 0, t65, 0.2e1 * t101, t65, 0.2e1 * qJD(3), 0.2e1 * t58 * qJD(3) + 0.2e1 * t101, -0.2e1 * t87, 0.2e1 * t81, 0, 0, 0, 0.2e1 * qJD(3) * t61 + 0.2e1 * t58 * t51, 0.2e1 * qJD(3) * t64 - 0.2e1 * t58 * t102, -0.2e1 * t55 * t87 - 0.2e1 * t56 * t89, 0.2e1 * t56 * t82 + 0.4e1 * t64 * t80, -0.2e1 * t61 * t94 - 0.2e1 * t63 * t81, 0.2e1 * t60 * t81 - 0.2e1 * t61 * t92, t45, -0.2e1 * t56 * qJD(2) * t60 + 0.2e1 * t29 * t51 + 0.2e1 * t7 * t61 + 0.2e1 * (-t56 * t52 + t60 * t87) * t57, 0.2e1 * t6 * t61 + 0.2e1 * (-qJD(2) * t63 + t91) * t56 + 0.2e1 * (-t28 + t44) * t51, -0.2e1 * t20 * t9, 0.2e1 * t20 * t11 - 0.2e1 * t121 * t9, -0.2e1 * t20 * t51 + 0.2e1 * t9 * t61, -0.2e1 * t11 * t61 - 0.2e1 * t121 * t51, t45, 0.2e1 * t30 * t11 + 0.2e1 * t121 * t18 + 0.2e1 * t2 * t61 + 0.2e1 * t77 * t51, 0.2e1 * t1 * t61 - 0.2e1 * t18 * t20 + 0.2e1 * t30 * t9 - 0.2e1 * t76 * t51; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, 0, 0, 0, 0, -t51, t102, 0, 0, 0, 0, 0, t24, t26, 0, 0, 0, 0, 0, t71, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108 * t52, t108 * t107, 0, 0, 0, 0, 0, t10 * t61 - t64 * t11, t8 * t61 - t64 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t51, 0, -t68, t69, -t64 * t82 - t80, t110 * t102 - 0.4e1 * t64 * t89, t26, -t24, 0, -t120 * t63 + t66 * t60, t120 * t60 + t66 * t63, t20 * t16 + t9 * t33, -t33 * t11 + t121 * t16 + t20 * t17 - t9 * t32, t72, -t71, 0, t50 * t11 + t121 * t98 + t13 * t61 + t30 * t17 + t18 * t32 + t75 * t51, t12 * t61 - t30 * t16 + t18 * t33 - t20 * t98 + t50 * t9 - t74 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t51, 0, 0, 0, 0, 0, t25, t27, 0, 0, 0, 0, 0, t32 * t102 - t64 * t17, t33 * t102 + t64 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t89, -0.2e1 * t82, 0, 0, 0, t60 * t99, t63 * t99, -0.2e1 * t33 * t16, 0.2e1 * t16 * t32 - 0.2e1 * t33 * t17, 0, 0, 0, 0.2e1 * t50 * t17 + 0.2e1 * t32 * t98, -0.2e1 * t50 * t16 + 0.2e1 * t33 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t27, t51, t7, t6, 0, 0, t9, -t11, t51, t62 * t97 + (t84 * t59 - t112) * qJD(6) + t85, -t118 + (-t4 - t97) * t59 + (t84 * t62 + t116) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t52, 0, 0, 0, 0, 0, t17, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t24, 0, 0, 0, 0, 0, t10, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t107, 0, -pkin(8) * t52, pkin(8) * t107, 0, 0, -t16, -t17, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t96, -0.2e1 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t11, t51, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
