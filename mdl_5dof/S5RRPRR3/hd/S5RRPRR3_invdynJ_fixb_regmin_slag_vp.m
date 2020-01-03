% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:37
% EndTime: 2020-01-03 12:00:39
% DurationCPUTime: 0.55s
% Computational Cost: add. (1140->128), mult. (1902->176), div. (0->0), fcn. (1145->14), ass. (0->97)
t71 = sin(qJ(2));
t128 = pkin(1) * t71;
t105 = qJD(1) * t128;
t75 = cos(qJ(2));
t109 = qJD(1) * t75;
t63 = qJD(1) + qJD(2);
t43 = pkin(1) * t109 + t63 * pkin(2);
t67 = sin(pkin(9));
t68 = cos(pkin(9));
t23 = -t67 * t105 + t68 * t43;
t21 = t63 * pkin(3) + t23;
t24 = t68 * t105 + t67 * t43;
t70 = sin(qJ(4));
t74 = cos(qJ(4));
t11 = t70 * t21 + t74 * t24;
t106 = qJDD(1) * t71;
t130 = pkin(1) * (qJD(2) * t109 + t106);
t120 = t75 * pkin(1);
t56 = qJDD(1) * t120;
t62 = qJDD(1) + qJDD(2);
t31 = t62 * pkin(2) - qJD(2) * t105 + t56;
t18 = -t67 * t130 + t68 * t31;
t13 = t62 * pkin(3) + t18;
t19 = t68 * t130 + t67 * t31;
t66 = qJ(1) + qJ(2);
t54 = pkin(9) + qJ(4) + t66;
t49 = sin(t54);
t50 = cos(t54);
t78 = -g(2) * t50 - g(3) * t49 - t11 * qJD(4) + t74 * t13 - t70 * t19;
t58 = qJDD(4) + t62;
t122 = t58 * pkin(4);
t91 = t122 + t78;
t127 = pkin(2) * t67;
t52 = t68 * pkin(2) + pkin(3);
t111 = t74 * t127 + t70 * t52;
t117 = t68 * t71;
t88 = pkin(1) * (-t67 * t75 - t117);
t35 = qJD(1) * t88;
t118 = t67 * t71;
t87 = pkin(1) * (t68 * t75 - t118);
t37 = qJD(1) * t87;
t59 = qJD(4) + t63;
t103 = (-t111 * qJD(4) - t74 * t35 + t70 * t37) * t59;
t89 = -t70 * t127 + t74 * t52;
t33 = -pkin(4) - t89;
t34 = pkin(8) + t111;
t77 = qJD(5) ^ 2;
t131 = t33 * t58 + t34 * t77 - t103;
t55 = pkin(2) + t120;
t96 = -pkin(1) * t118 + t68 * t55;
t32 = pkin(3) + t96;
t39 = pkin(1) * t117 + t67 * t55;
t114 = t70 * t32 + t74 * t39;
t116 = t70 * t24;
t79 = g(2) * t49 - (qJD(4) * t21 + t19) * t74 + qJD(4) * t116 - t70 * t13 - g(3) * t50;
t36 = qJD(2) * t88;
t38 = qJD(2) * t87;
t123 = (t114 * qJD(4) - t74 * t36 + t70 * t38) * t59;
t121 = t59 * pkin(4);
t119 = t11 * t59;
t73 = cos(qJ(5));
t115 = t73 * t58;
t113 = -t89 * qJD(4) + t70 * t35 + t74 * t37;
t69 = sin(qJ(5));
t64 = t69 ^ 2;
t110 = -t73 ^ 2 + t64;
t108 = t73 * qJD(5);
t107 = qJDD(3) - g(1);
t60 = sin(t66);
t61 = cos(t66);
t104 = g(2) * t60 - g(3) * t61;
t10 = t74 * t21 - t116;
t8 = -t10 - t121;
t102 = t8 * t108 - t69 * t91;
t98 = qJD(1) * (-qJD(2) + t63);
t97 = qJD(2) * (-qJD(1) - t63);
t94 = -g(2) * t61 - g(3) * t60;
t93 = t74 * t32 - t70 * t39;
t90 = t56 + t94;
t85 = pkin(8) * t77 - t119 - t122;
t14 = -pkin(4) - t93;
t15 = pkin(8) + t114;
t84 = t14 * t58 + t15 * t77 + t123;
t83 = -t58 * pkin(8) - t8 * t59 + t79;
t82 = -pkin(8) * qJDD(5) + (t10 - t121) * qJD(5);
t4 = t93 * qJD(4) + t70 * t36 + t74 * t38;
t81 = -qJDD(5) * t15 + (t14 * t59 - t4) * qJD(5);
t80 = -qJDD(5) * t34 + (t33 * t59 + t113) * qJD(5);
t76 = cos(qJ(1));
t72 = sin(qJ(1));
t57 = t59 ^ 2;
t45 = qJDD(5) * t73 - t77 * t69;
t44 = qJDD(5) * t69 + t77 * t73;
t26 = 0.2e1 * t69 * t59 * t108 + t64 * t58;
t20 = -0.2e1 * t110 * t59 * qJD(5) + 0.2e1 * t69 * t115;
t6 = t8 * qJD(5) * t69;
t1 = [qJDD(1), -g(2) * t76 - g(3) * t72, g(2) * t72 - g(3) * t76, t62, (t62 * t75 + t71 * t97) * pkin(1) + t90, ((-qJDD(1) - t62) * t71 + t75 * t97) * pkin(1) + t104, t19 * t39 + t24 * t38 + t18 * t96 + t23 * t36 - g(2) * (t76 * pkin(1) + pkin(2) * t61) - g(3) * (t72 * pkin(1) + pkin(2) * t60), t58, t93 * t58 - t123 + t78, -t114 * t58 - t4 * t59 + t79, t26, t20, t44, t45, 0, t6 + t81 * t69 + (-t84 + t91) * t73, t84 * t69 + t81 * t73 + t102; 0, 0, 0, t62, t98 * t128 + t90, (t75 * t98 - t106) * pkin(1) + t104, -t23 * t35 - t24 * t37 + (t18 * t68 + t19 * t67 + t94) * pkin(2), t58, t89 * t58 + t103 + t78, -t111 * t58 + t113 * t59 + t79, t26, t20, t44, t45, 0, t6 + t80 * t69 + (t91 - t131) * t73, t131 * t69 + t80 * t73 + t102; 0, 0, 0, 0, 0, 0, t107, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t44; 0, 0, 0, 0, 0, 0, 0, t58, t78 + t119, t10 * t59 + t79, t26, t20, t44, t45, 0, t6 + t82 * t69 + (-t85 + t91) * t73, t85 * t69 + t82 * t73 + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 * t57 * t73, t110 * t57, t69 * t58, t115, qJDD(5), t107 * t73 + t83 * t69, -t107 * t69 + t83 * t73;];
tau_reg = t1;
