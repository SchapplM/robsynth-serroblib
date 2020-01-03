% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPR8
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:44
% EndTime: 2019-12-31 17:42:46
% DurationCPUTime: 0.60s
% Computational Cost: add. (815->131), mult. (1518->188), div. (0->0), fcn. (1159->14), ass. (0->96)
t75 = sin(pkin(8));
t77 = cos(pkin(8));
t100 = g(1) * t77 + g(2) * t75;
t73 = qJ(2) + qJ(3);
t66 = sin(t73);
t67 = cos(t73);
t129 = -g(3) * t67 + t100 * t66;
t106 = qJD(1) * qJD(2);
t83 = cos(qJ(2));
t64 = t83 * qJDD(1);
t80 = sin(qJ(2));
t40 = qJDD(2) * pkin(2) - t80 * t106 + t64;
t79 = sin(qJ(3));
t82 = cos(qJ(3));
t56 = qJD(2) * pkin(2) + t83 * qJD(1);
t95 = -t80 * qJDD(1) - t83 * t106;
t89 = qJD(3) * t56 - t95;
t127 = t79 * t40 + t89 * t82;
t109 = qJD(1) * t80;
t102 = qJD(3) * t109;
t55 = t79 * t102;
t15 = t127 - t55;
t74 = sin(pkin(9));
t76 = cos(pkin(9));
t69 = qJDD(2) + qJDD(3);
t35 = t82 * t40;
t86 = -t82 * t102 - t89 * t79 + t35;
t9 = t69 * pkin(3) + t86;
t4 = -t74 * t15 + t76 * t9;
t65 = pkin(9) + t73;
t59 = cos(t65);
t130 = -t69 * pkin(4) + g(3) * t59 - t4;
t58 = sin(t65);
t96 = t100 * t58;
t70 = qJD(2) + qJD(3);
t110 = pkin(2) * qJD(3);
t117 = t76 * t79;
t44 = t79 * t83 + t82 * t80;
t38 = t44 * qJD(1);
t43 = -t79 * t80 + t82 * t83;
t39 = t43 * qJD(1);
t114 = t76 * t38 + t74 * t39 - (t74 * t82 + t117) * t110;
t119 = t74 * t79;
t63 = t82 * pkin(2) + pkin(3);
t98 = -pkin(2) * t119 + t76 * t63;
t33 = -pkin(4) - t98;
t112 = pkin(2) * t117 + t74 * t63;
t34 = pkin(7) + t112;
t84 = qJD(5) ^ 2;
t128 = t114 * t70 - t33 * t69 - t34 * t84;
t126 = pkin(2) * t69;
t5 = t76 * t15 + t74 * t9;
t31 = t82 * t109 + t79 * t56;
t120 = t74 * t31;
t118 = t76 * t31;
t81 = cos(qJ(5));
t115 = t81 * t69;
t113 = t74 * t38 - t76 * t39 + (t76 * t82 - t119) * t110;
t78 = sin(qJ(5));
t71 = t78 ^ 2;
t111 = -t81 ^ 2 + t71;
t108 = qJD(5) * t81;
t107 = qJDD(1) - g(3);
t30 = -t79 * t109 + t82 * t56;
t28 = t70 * pkin(3) + t30;
t16 = t76 * t28 - t120;
t13 = -t70 * pkin(4) - t16;
t105 = t13 * t108 + t130 * t78;
t104 = t13 * qJD(5) * t78 + t81 * t96;
t101 = g(3) * t66 + t100 * t67 + t55;
t97 = g(1) * t75 - g(2) * t77 - qJDD(4);
t22 = -t76 * t43 + t74 * t44;
t23 = t74 * t43 + t76 * t44;
t24 = t70 * t43;
t25 = t70 * t44;
t6 = t74 * t24 + t76 * t25;
t94 = t22 * t69 + t23 * t84 + t6 * t70;
t18 = t74 * t30 + t118;
t60 = t74 * pkin(3) + pkin(7);
t61 = -t76 * pkin(3) - pkin(4);
t93 = -t18 * t70 + t60 * t84 + t61 * t69;
t7 = t76 * t24 - t74 * t25;
t92 = -qJDD(5) * t23 + (t22 * t70 - t7) * qJD(5);
t19 = t76 * t30 - t120;
t91 = -qJDD(5) * t60 + (t61 * t70 + t19) * qJD(5);
t90 = -qJDD(5) * t34 + (t33 * t70 - t113) * qJD(5);
t88 = -t69 * pkin(7) + g(3) * t58 + t100 * t59 - t13 * t70 - t5;
t87 = (-pkin(2) * t70 - t56) * qJD(3) + t95;
t85 = qJD(2) ^ 2;
t68 = t70 ^ 2;
t48 = qJDD(5) * t81 - t84 * t78;
t47 = qJDD(5) * t78 + t84 * t81;
t32 = 0.2e1 * t78 * t70 * t108 + t71 * t69;
t26 = -0.2e1 * t111 * t70 * qJD(5) + 0.2e1 * t78 * t115;
t17 = t74 * t28 + t118;
t1 = [t107, 0, t83 * qJDD(2) - t85 * t80, -qJDD(2) * t80 - t85 * t83, 0, -t25 * t70 + t43 * t69, -t24 * t70 - t44 * t69, -t16 * t6 + t17 * t7 - t4 * t22 + t5 * t23 - g(3), 0, 0, 0, 0, 0, t92 * t78 - t94 * t81, t94 * t78 + t92 * t81; 0, qJDD(2), -g(3) * t83 + t100 * t80 + t64, t100 * t83 - t107 * t80, t69, t38 * t70 + t35 + (-t102 + t126) * t82 + t87 * t79 + t129, t39 * t70 + (-t40 - t126) * t79 + t87 * t82 + t101, t5 * t112 + t4 * t98 - g(3) * (t83 * pkin(2) + pkin(3) * t67) - t100 * (-t80 * pkin(2) - pkin(3) * t66) + t113 * t17 + t114 * t16, t32, t26, t47, t48, 0, t90 * t78 + (-t130 + t128) * t81 + t104, t90 * t81 + (-t128 - t96) * t78 + t105; 0, 0, 0, 0, t69, t31 * t70 + t129 + t86, t30 * t70 + t101 - t127, t16 * t18 - t17 * t19 + (t4 * t76 + t5 * t74 + t129) * pkin(3), t32, t26, t47, t48, 0, t91 * t78 + (-t130 - t93) * t81 + t104, t91 * t81 + (-t96 + t93) * t78 + t105; 0, 0, 0, 0, 0, 0, 0, -t97, 0, 0, 0, 0, 0, t48, -t47; 0, 0, 0, 0, 0, 0, 0, 0, -t78 * t68 * t81, t111 * t68, t78 * t69, t115, qJDD(5), t88 * t78 - t97 * t81, t97 * t78 + t88 * t81;];
tau_reg = t1;
