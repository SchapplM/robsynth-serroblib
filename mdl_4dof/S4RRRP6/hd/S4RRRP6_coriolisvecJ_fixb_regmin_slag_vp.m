% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauc_reg [4x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:10
% EndTime: 2019-12-31 17:19:13
% DurationCPUTime: 0.74s
% Computational Cost: add. (694->166), mult. (1832->266), div. (0->0), fcn. (1072->4), ass. (0->96)
t84 = (qJD(1) * qJD(2));
t121 = -2 * t84;
t58 = sin(qJ(3));
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t83 = qJD(2) * qJD(3);
t88 = t58 * qJD(2);
t60 = cos(qJ(3));
t90 = qJD(3) * t60;
t17 = (t59 * t90 + t61 * t88) * qJD(1) + t58 * t83;
t93 = qJD(1) * t59;
t36 = t60 * t93 + t88;
t120 = t36 ^ 2;
t85 = t61 * qJD(1);
t50 = -qJD(3) + t85;
t40 = -t61 * pkin(2) - t59 * pkin(6) - pkin(1);
t28 = t40 * qJD(1);
t54 = pkin(5) * t85;
t45 = qJD(2) * pkin(6) + t54;
t12 = t60 * t28 - t58 * t45;
t6 = -t36 * qJ(4) + t12;
t5 = -t50 * pkin(3) + t6;
t119 = t5 - t6;
t94 = t60 * qJ(4);
t67 = pkin(3) * t59 - t61 * t94;
t103 = -qJ(4) - pkin(6);
t75 = qJD(3) * t103;
t70 = pkin(2) * t59 - pkin(6) * t61;
t38 = t70 * qJD(1);
t80 = t58 * t93;
t98 = pkin(5) * t80 + t60 * t38;
t118 = -t67 * qJD(1) - t58 * qJD(4) + t60 * t75 - t98;
t78 = t61 * t84;
t16 = qJD(3) * t80 + (-t78 - t83) * t60;
t117 = t16 * t58;
t87 = t60 * qJD(2);
t34 = t80 - t87;
t116 = t34 * t50;
t115 = t36 * t50;
t44 = -qJD(2) * pkin(2) + pkin(5) * t93;
t114 = t44 * t58;
t113 = t44 * t60;
t112 = t50 * t58;
t111 = t50 * t60;
t110 = t50 * t61;
t109 = t58 * t28;
t108 = t59 * t60;
t107 = t60 * t61;
t63 = qJD(1) ^ 2;
t106 = t61 * t63;
t62 = qJD(2) ^ 2;
t105 = t62 * t59;
t104 = t62 * t61;
t39 = t70 * qJD(2);
t29 = qJD(1) * t39;
t76 = t59 * t84;
t71 = pkin(5) * t76;
t102 = -t60 * t29 - t58 * t71;
t24 = t58 * t38;
t86 = t60 * qJD(4);
t95 = t58 * qJ(4);
t101 = t58 * t75 + t86 - t24 - (-pkin(5) * t108 - t61 * t95) * qJD(1);
t100 = t58 * t39 + t40 * t90;
t99 = t59 * pkin(5) * t88 + t60 * t39;
t51 = pkin(5) * t107;
t97 = t58 * t40 + t51;
t56 = t59 ^ 2;
t96 = -t61 ^ 2 + t56;
t92 = qJD(2) * t61;
t91 = qJD(3) * t58;
t89 = t44 * qJD(3);
t82 = t59 * t91;
t81 = t50 * t90;
t79 = pkin(3) * t58 + pkin(5);
t9 = t17 * pkin(3) + pkin(5) * t78;
t74 = t34 + t87;
t73 = -t36 + t88;
t72 = pkin(1) * t121;
t13 = t60 * t45 + t109;
t69 = qJD(1) * t56 - t110;
t68 = t28 * t90 + t58 * t29 - t45 * t91;
t65 = -t13 * qJD(3) - t102;
t64 = -t60 * t71 + t68;
t43 = t103 * t60;
t42 = t103 * t58;
t33 = t60 * t40;
t31 = t34 ^ 2;
t15 = t34 * pkin(3) + qJD(4) + t44;
t14 = -t59 * t95 + t97;
t11 = -t59 * t94 + t33 + (-pkin(5) * t58 - pkin(3)) * t61;
t7 = -t34 * qJ(4) + t13;
t4 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t108 + (-qJD(4) * t59 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t61) * t58 + t100;
t3 = -t59 * t86 + t67 * qJD(2) + (-t51 + (qJ(4) * t59 - t40) * t58) * qJD(3) + t99;
t2 = -t17 * qJ(4) - t34 * qJD(4) + t64;
t1 = pkin(3) * t76 + t16 * qJ(4) - t36 * qJD(4) + t65;
t8 = [0, 0, 0, 0.2e1 * t61 * t76, t96 * t121, t104, -t105, 0, -pkin(5) * t104 + t59 * t72, pkin(5) * t105 + t61 * t72, -t16 * t108 + (t61 * t87 - t82) * t36, (-t34 * t60 - t36 * t58) * t92 + (t117 - t17 * t60 + (t34 * t58 - t36 * t60) * qJD(3)) * t59, t50 * t82 + t16 * t61 + (t36 * t59 + t69 * t60) * qJD(2), t59 * t81 + t17 * t61 + (-t34 * t59 - t69 * t58) * qJD(2), (-t50 - t85) * t59 * qJD(2), -(-t40 * t91 + t99) * t50 + (t60 * t89 + pkin(5) * t17 + (qJD(1) * t33 + t12) * qJD(2)) * t59 + ((pkin(5) * t34 + t114) * qJD(2) + (t109 + (pkin(5) * t50 + t45) * t60) * qJD(3) + t102) * t61, (-t61 * pkin(5) * t91 + t100) * t50 + t68 * t61 + (-pkin(5) * t16 - t58 * t89) * t59 + ((pkin(5) * t36 + t113) * t61 + (-pkin(5) * t111 - t97 * qJD(1) - t13) * t59) * qJD(2), t11 * t16 - t14 * t17 - t3 * t36 - t4 * t34 + (-t5 * t60 - t58 * t7) * t92 + (-t1 * t60 - t2 * t58 + (t5 * t58 - t60 * t7) * qJD(3)) * t59, t1 * t11 + t2 * t14 + t5 * t3 + t7 * t4 + t15 * t79 * t92 + (t15 * pkin(3) * t90 + t9 * t79) * t59; 0, 0, 0, -t59 * t106, t96 * t63, 0, 0, 0, t63 * pkin(1) * t59, pkin(1) * t106, -t36 * t111 - t117, (-t16 + t116) * t60 + (-t17 + t115) * t58, -t81 + (t50 * t107 + t73 * t59) * qJD(1), t50 * t91 + (-t58 * t110 + t74 * t59) * qJD(1), t50 * t93, -pkin(2) * t17 + t98 * t50 + (pkin(6) * t111 + t114) * qJD(3) + ((-pkin(6) * t88 - t12) * t59 + (-t74 * pkin(5) - t114) * t61) * qJD(1), pkin(2) * t16 - t24 * t50 + (-pkin(6) * t112 + t113) * qJD(3) + (-t44 * t107 + (-pkin(6) * t87 + t13) * t59 + (t50 * t108 + t73 * t61) * pkin(5)) * qJD(1), t42 * t16 + t43 * t17 - t118 * t36 - t101 * t34 + (t50 * t5 + t2) * t60 + (t50 * t7 - t1) * t58, -t2 * t43 + t1 * t42 + t9 * (-t60 * pkin(3) - pkin(2)) + t101 * t7 + t118 * t5 + (-pkin(3) * t112 - t54) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * t34, -t31 + t120, -t16 - t116, -t115 - t17, t76, -t13 * t50 - t44 * t36 + t65, -t12 * t50 + t44 * t34 - t64, pkin(3) * t16 - t119 * t34, t119 * t7 + (-t15 * t36 + t1) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 - t120, t7 * t34 + t5 * t36 + t9;];
tauc_reg = t8;
