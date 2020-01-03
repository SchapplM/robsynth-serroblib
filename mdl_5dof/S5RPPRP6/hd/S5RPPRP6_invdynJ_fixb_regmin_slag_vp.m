% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:19
% EndTime: 2019-12-31 17:55:21
% DurationCPUTime: 0.74s
% Computational Cost: add. (1085->191), mult. (1992->213), div. (0->0), fcn. (1300->8), ass. (0->108)
t81 = sin(qJ(1));
t83 = cos(qJ(1));
t146 = g(1) * t81 - g(2) * t83;
t76 = sin(pkin(7));
t77 = cos(pkin(7));
t80 = sin(qJ(4));
t82 = cos(qJ(4));
t42 = t82 * t76 + t80 * t77;
t142 = t42 * qJD(1);
t147 = t142 * qJD(4);
t128 = t76 ^ 2 + t77 ^ 2;
t73 = qJDD(1) * qJ(2);
t74 = qJD(1) * qJD(2);
t145 = t73 + t74;
t79 = -pkin(1) - qJ(3);
t50 = t79 * qJD(1) + qJD(2);
t144 = t128 * t50;
t105 = g(1) * t83 + g(2) * t81;
t47 = qJDD(3) + t145;
t143 = t47 - t105;
t123 = qJD(4) * t82;
t124 = qJD(4) * t80;
t37 = -t76 * t123 - t77 * t124;
t132 = t82 * t77;
t41 = t80 * t76 - t132;
t131 = t37 * qJD(4) - t41 * qJDD(4);
t141 = -qJD(1) * t142 + t131;
t138 = -qJD(1) * qJD(3) + qJDD(1) * t79;
t43 = qJDD(2) + t138;
t107 = -pkin(6) * qJDD(1) + t43;
t27 = t107 * t76;
t28 = t107 * t77;
t109 = -pkin(6) * qJD(1) + t50;
t30 = t109 * t77;
t115 = t30 * t123 + t82 * t27 + t80 * t28;
t72 = pkin(7) + qJ(4);
t62 = sin(t72);
t63 = cos(t72);
t140 = -g(3) * t63 - t146 * t62 + t115;
t135 = -pkin(6) + t79;
t44 = t135 * t76;
t45 = t135 * t77;
t19 = t82 * t44 + t80 * t45;
t18 = t80 * t44 - t82 * t45;
t8 = -t42 * qJD(3) - t18 * qJD(4);
t139 = -t8 * qJD(4) - t19 * qJDD(4) - t105 * t63;
t125 = qJD(1) * t76;
t114 = t80 * t125;
t93 = -qJD(4) * t114 + t42 * qJDD(1);
t113 = qJD(1) * t132;
t35 = t113 - t114;
t137 = t35 ^ 2;
t136 = 0.2e1 * t74;
t64 = t76 * pkin(3);
t134 = t35 * t142;
t29 = t109 * t76;
t133 = t80 * t29;
t130 = t83 * pkin(1) + t81 * qJ(2);
t56 = qJ(2) + t64;
t127 = pkin(1) * qJDD(1);
t122 = qJDD(4) * pkin(4);
t13 = t82 * t29 + t80 * t30;
t121 = t13 * qJD(4);
t61 = qJD(1) * qJ(2) + qJD(3);
t12 = t82 * t30 - t133;
t120 = qJD(5) - t12;
t119 = t76 * qJDD(1);
t118 = t77 * qJDD(1);
t116 = qJDD(4) * qJ(5);
t46 = pkin(3) * t125 + t61;
t112 = g(2) * t130;
t111 = t128 * t43;
t110 = t12 + t133;
t108 = t29 * t123 + t30 * t124 + t80 * t27 - t82 * t28;
t40 = pkin(3) * t119 + t47;
t106 = qJDD(2) - t127;
t103 = -t82 * t118 + t80 * t119;
t101 = -t62 * pkin(4) + t63 * qJ(5);
t16 = t103 + t147;
t100 = t41 * t16 + t35 * t37;
t38 = t77 * t123 - t76 * t124;
t97 = -t38 * qJD(4) - t42 * qJDD(4);
t94 = -t101 + t64;
t92 = qJD(1) * t35 - t97;
t91 = g(3) * t62 - t146 * t63 - t108;
t17 = qJD(4) * t113 + t93;
t90 = t17 * pkin(4) + t16 * qJ(5) + t40;
t1 = t116 + (qJD(5) - t133) * qJD(4) + t115;
t10 = qJD(4) * qJ(5) + t13;
t3 = qJDD(5) + t108 - t122;
t7 = -qJD(4) * pkin(4) + t120;
t89 = t1 * t42 + t10 * t38 + t3 * t41 - t7 * t37 - t146;
t88 = t143 + t145;
t9 = -t41 * qJD(3) + t19 * qJD(4);
t87 = -t9 * qJD(4) - t18 * qJDD(4) - t105 * t62;
t11 = pkin(4) * t142 - t35 * qJ(5) + t46;
t86 = t11 * t35 + qJDD(5) - t91;
t85 = (t35 + t113) * qJD(4) + t93;
t84 = qJD(1) ^ 2;
t78 = -pkin(6) - qJ(3);
t66 = t83 * qJ(2);
t32 = t142 ^ 2;
t15 = t35 * pkin(4) + qJ(5) * t142;
t14 = t42 * pkin(4) + t41 * qJ(5) + t56;
t5 = t103 + 0.2e1 * t147;
t4 = t38 * pkin(4) - t37 * qJ(5) + t41 * qJD(5) + qJD(2);
t2 = -t35 * qJD(5) + t90;
t6 = [qJDD(1), t146, t105, qJDD(2) - 0.2e1 * t127 - t146, -t105 + 0.2e1 * t73 + t136, -t106 * pkin(1) - g(1) * (-t81 * pkin(1) + t66) - t112 + (t73 + t136) * qJ(2), t88 * t76, t88 * t77, t146 + t128 * (-t138 - t43), t47 * qJ(2) + t61 * qJD(2) - g(1) * (t79 * t81 + t66) - g(2) * (t83 * qJ(3) + t130) + t79 * t111 - qJD(3) * t144, t100, -t142 * t37 + t16 * t42 + t41 * t17 - t35 * t38, t131, t97, 0, qJD(2) * t142 + t56 * t17 + t46 * t38 + t40 * t42 + t87, qJD(2) * t35 - t56 * t16 + t46 * t37 - t40 * t41 + t139, t11 * t38 + t14 * t17 + t142 * t4 + t2 * t42 + t87, -t142 * t8 - t18 * t16 - t19 * t17 + t9 * t35 - t89, -t11 * t37 + t14 * t16 + t2 * t41 - t4 * t35 - t139, t1 * t19 + t10 * t8 + t2 * t14 + t11 * t4 + t3 * t18 + t7 * t9 - g(1) * t66 - t112 + (-g(1) * t94 + g(2) * t78) * t83 + (-g(1) * (-pkin(1) + t78) - g(2) * t94) * t81; 0, 0, 0, qJDD(1), -t84, -t84 * qJ(2) + t106 - t146, -t84 * t76, -t84 * t77, -t128 * qJDD(1), -t61 * qJD(1) + t111 - t146, 0, 0, 0, 0, 0, t141, -t92, t141, -t142 * t38 - t42 * t17 - t100, t92, -t11 * qJD(1) + t89; 0, 0, 0, 0, 0, 0, t119, t118, -t128 * t84, qJD(1) * t144 + t143, 0, 0, 0, 0, 0, t85, -t5, t85, -t32 - t137, t5, t10 * t142 + (-qJD(5) - t7) * t35 + t90 - t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, -t32 + t137, -t103, (t35 - t113) * qJD(4) - t93, qJDD(4), -t46 * t35 + t121 + t91, t110 * qJD(4) + t142 * t46 - t140, -t142 * t15 + t121 + 0.2e1 * t122 - t86, pkin(4) * t16 - t17 * qJ(5) + (t10 - t13) * t35 + (t7 - t120) * t142, 0.2e1 * t116 - t11 * t142 + t15 * t35 + (0.2e1 * qJD(5) - t110) * qJD(4) + t140, -t3 * pkin(4) - g(3) * t101 + t1 * qJ(5) + t120 * t10 - t11 * t15 - t7 * t13 - t146 * (pkin(4) * t63 + qJ(5) * t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t134, -t103, -qJD(4) ^ 2 - t137, -t10 * qJD(4) - t122 + t86;];
tau_reg = t6;
