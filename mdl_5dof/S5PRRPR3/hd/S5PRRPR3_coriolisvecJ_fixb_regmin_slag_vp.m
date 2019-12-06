% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:40
% EndTime: 2019-12-05 16:19:44
% DurationCPUTime: 0.65s
% Computational Cost: add. (870->143), mult. (2285->214), div. (0->0), fcn. (1653->6), ass. (0->99)
t90 = sin(qJ(5));
t111 = qJD(5) * t90;
t88 = sin(pkin(9));
t89 = cos(pkin(9));
t91 = sin(qJ(3));
t93 = cos(qJ(3));
t73 = -t88 * t91 + t89 * t93;
t68 = t73 * qJD(2);
t92 = cos(qJ(5));
t59 = t92 * t68;
t74 = t88 * t93 + t89 * t91;
t69 = t74 * qJD(3);
t62 = qJD(2) * t69;
t110 = qJD(2) * qJD(3);
t104 = t93 * t110;
t105 = t91 * t110;
t63 = t104 * t89 - t105 * t88;
t70 = t74 * qJD(2);
t1 = qJD(5) * t59 - t111 * t70 - t62 * t90 + t63 * t92;
t29 = -t70 * t90 + t59;
t85 = qJD(3) + qJD(5);
t120 = t29 * t85;
t134 = t1 - t120;
t98 = t68 * t90 + t70 * t92;
t133 = t98 * t29;
t121 = t98 * t85;
t2 = qJD(5) * t98 + t62 * t92 + t90 * t63;
t132 = -t2 + t121;
t131 = -t29 ^ 2 + t98 ^ 2;
t124 = t68 * pkin(7);
t115 = -qJ(4) - pkin(6);
t80 = t115 * t93;
t66 = t91 * qJD(1) - qJD(2) * t80;
t119 = t89 * t66;
t113 = qJD(3) * pkin(3);
t79 = t115 * t91;
t64 = qJD(1) * t93 + qJD(2) * t79;
t58 = t64 + t113;
t19 = t58 * t88 + t119;
t11 = t19 + t124;
t107 = -pkin(3) * t93 - pkin(2);
t100 = t107 * qJD(2);
t78 = qJD(4) + t100;
t36 = -t68 * pkin(4) + t78;
t130 = t11 * t111 - t36 * t29;
t128 = -0.2e1 * t110;
t127 = qJD(5) - t85;
t109 = qJD(3) * qJD(1);
t102 = qJD(3) * t115;
t65 = t93 * qJD(4) + t102 * t91;
t38 = qJD(2) * t65 + t109 * t93;
t67 = -t91 * qJD(4) + t102 * t93;
t39 = qJD(2) * t67 - t109 * t91;
t12 = -t38 * t88 + t39 * t89;
t4 = -pkin(7) * t63 + t12;
t13 = t38 * t89 + t39 * t88;
t5 = -pkin(7) * t62 + t13;
t126 = -t36 * t98 + t4 * t92 - t90 * t5;
t125 = pkin(3) * t88;
t72 = t73 * qJD(3);
t97 = t73 * t92 - t74 * t90;
t7 = qJD(5) * t97 - t90 * t69 + t92 * t72;
t123 = t7 * t85;
t122 = t70 * pkin(7);
t51 = t88 * t66;
t95 = qJD(2) ^ 2;
t118 = t93 * t95;
t94 = qJD(3) ^ 2;
t117 = t94 * t91;
t116 = t94 * t93;
t22 = t64 * t89 - t51;
t23 = t65 * t89 + t67 * t88;
t41 = t79 * t88 - t80 * t89;
t114 = t91 ^ 2 - t93 ^ 2;
t112 = qJD(2) * t91;
t108 = t91 * t113;
t18 = t58 * t89 - t51;
t20 = -t64 * t88 - t119;
t21 = -t65 * t88 + t67 * t89;
t40 = t79 * t89 + t80 * t88;
t101 = pkin(2) * t128;
t10 = qJD(3) * pkin(4) - t122 + t18;
t99 = -t90 * t10 - t92 * t11;
t32 = t73 * t90 + t74 * t92;
t83 = pkin(3) * t89 + pkin(4);
t82 = pkin(3) * t105;
t50 = -pkin(4) * t73 + t107;
t43 = pkin(4) * t69 + t108;
t42 = pkin(3) * t112 + pkin(4) * t70;
t37 = pkin(4) * t62 + t82;
t25 = pkin(7) * t73 + t41;
t24 = -pkin(7) * t74 + t40;
t17 = -pkin(7) * t69 + t23;
t16 = t22 - t122;
t15 = -pkin(7) * t72 + t21;
t14 = t20 - t124;
t8 = qJD(5) * t32 + t69 * t92 + t90 * t72;
t6 = t8 * t85;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, -t116, -t62 * t74 - t63 * t73 + t68 * t72 + t69 * t70, t12 * t73 + t13 * t74 - t18 * t69 + t19 * t72, 0, 0, 0, 0, 0, -t6, -t123; 0, 0, 0, 0, 0.2e1 * t91 * t104, t114 * t128, t116, -t117, 0, -pkin(6) * t116 + t101 * t91, pkin(6) * t117 + t101 * t93, -t12 * t74 + t13 * t73 - t18 * t72 - t19 * t69 - t21 * t70 + t23 * t68 - t40 * t63 - t41 * t62, t12 * t40 + t13 * t41 + t18 * t21 + t19 * t23 + (t78 + t100) * t108, t1 * t32 + t7 * t98, t1 * t97 - t2 * t32 + t29 * t7 - t8 * t98, t123, -t6, 0, -t43 * t29 + t50 * t2 - t37 * t97 + t36 * t8 + (t92 * t15 - t90 * t17 + (-t24 * t90 - t25 * t92) * qJD(5)) * t85, t43 * t98 + t50 * t1 + t37 * t32 + t36 * t7 - (t90 * t15 + t92 * t17 + (t24 * t92 - t25 * t90) * qJD(5)) * t85; 0, 0, 0, 0, -t91 * t118, t114 * t95, 0, 0, 0, t95 * pkin(2) * t91, pkin(2) * t118, (t19 + t20) * t70 + (t18 - t22) * t68 + (-t62 * t88 - t63 * t89) * pkin(3), -t18 * t20 - t19 * t22 + (-t112 * t78 + t12 * t89 + t13 * t88) * pkin(3), -t133, t131, t134, t132, 0, t42 * t29 - (t92 * t14 - t90 * t16) * t85 + ((-t125 * t92 - t83 * t90) * t85 + t99) * qJD(5) + t126, -t92 * t5 - t90 * t4 - t42 * t98 + (t90 * t14 + t92 * t16) * t85 + (-(-t125 * t90 + t83 * t92) * t85 - t92 * t10) * qJD(5) + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68 ^ 2 - t70 ^ 2, t18 * t70 - t19 * t68 + t82, 0, 0, 0, 0, 0, t2 + t121, t1 + t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, t131, t134, t132, 0, t127 * t99 + t126, (-t11 * t85 - t4) * t90 + (-t10 * t127 - t5) * t92 + t130;];
tauc_reg = t3;
