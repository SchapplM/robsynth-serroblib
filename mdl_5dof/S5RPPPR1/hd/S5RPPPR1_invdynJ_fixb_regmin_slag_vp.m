% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:44
% EndTime: 2022-01-20 09:12:48
% DurationCPUTime: 0.97s
% Computational Cost: add. (816->177), mult. (1679->260), div. (0->0), fcn. (1228->14), ass. (0->121)
t77 = qJ(1) + pkin(7);
t70 = sin(t77);
t144 = g(1) * t70;
t72 = cos(t77);
t67 = g(2) * t72;
t111 = t67 - t144;
t83 = cos(pkin(7));
t142 = t83 * pkin(1);
t68 = -pkin(2) - t142;
t123 = qJDD(1) * t68;
t56 = qJDD(3) + t123;
t148 = -t111 - t56;
t79 = sin(pkin(8));
t129 = qJD(1) * t79;
t78 = sin(pkin(9));
t116 = t78 * t129;
t84 = sin(qJ(5));
t109 = t84 * t116;
t81 = cos(pkin(9));
t86 = cos(qJ(5));
t132 = t86 * t81;
t115 = qJD(5) * t132;
t96 = t86 * t78 + t84 * t81;
t11 = -qJD(5) * t109 + (qJD(1) * t115 + qJDD(1) * t96) * t79;
t82 = cos(pkin(8));
t62 = t82 * qJD(1) - qJD(5);
t147 = qJD(5) + t62;
t80 = sin(pkin(7));
t65 = t80 * pkin(1) + qJ(3);
t52 = qJD(1) * qJD(3) + qJDD(1) * t65;
t145 = pkin(6) * t79;
t143 = g(3) * t79;
t118 = t79 * t132;
t122 = qJDD(1) * t78;
t91 = t96 * qJD(5);
t10 = qJDD(1) * t118 + (-qJD(1) * t91 - t84 * t122) * t79;
t141 = t10 * t82;
t35 = t96 * t129;
t140 = t35 * t62;
t37 = qJD(1) * t118 - t109;
t139 = t37 * t62;
t138 = t65 * t82;
t137 = t70 * t82;
t136 = t72 * t82;
t135 = t82 * t11;
t88 = qJD(1) ^ 2;
t134 = t82 * t88;
t133 = t84 * t78;
t39 = (-qJD(5) * t133 + t115) * t79;
t43 = t96 * t79;
t121 = t82 * qJDD(1);
t61 = -qJDD(5) + t121;
t131 = t39 * t62 + t43 * t61;
t126 = qJD(4) * t79;
t93 = -t82 * pkin(3) - t79 * qJ(4) - pkin(2);
t51 = t93 - t142;
t19 = -qJD(1) * t126 + qJDD(1) * t51 + qJDD(3);
t32 = t79 * qJDD(2) + t82 * t52;
t6 = t78 * t19 + t81 * t32;
t33 = qJD(1) * t51 + qJD(3);
t60 = t65 * qJD(1);
t42 = t79 * qJD(2) + t82 * t60;
t9 = t78 * t33 + t81 * t42;
t18 = t81 * t138 + t78 * t51;
t74 = t79 ^ 2;
t130 = t82 ^ 2 + t74;
t128 = qJD(3) * t79;
t127 = qJD(3) * t82;
t125 = qJDD(2) - g(3);
t119 = t81 * t145;
t87 = cos(qJ(1));
t117 = t87 * pkin(1) + t72 * pkin(2) + t70 * qJ(3);
t5 = t81 * t19 - t78 * t32;
t94 = -t82 * pkin(4) - t119;
t2 = qJDD(1) * t94 + t5;
t113 = t79 * t122;
t3 = -pkin(6) * t113 + t6;
t114 = t86 * t2 - t84 * t3;
t85 = sin(qJ(1));
t112 = -t85 * pkin(1) + t72 * qJ(3);
t110 = t130 * t88;
t8 = t81 * t33 - t78 * t42;
t41 = t82 * qJD(2) - t79 * t60;
t48 = t79 * t52;
t31 = t82 * qJDD(2) - t48;
t107 = -g(1) * t72 - g(2) * t70;
t105 = g(1) * t85 - g(2) * t87;
t104 = t84 * t2 + t86 * t3;
t4 = qJD(1) * t94 + t8;
t7 = -pkin(6) * t116 + t9;
t103 = t86 * t4 - t84 * t7;
t102 = -t84 * t4 - t86 * t7;
t46 = t81 * t51;
t12 = -t119 + t46 + (-t65 * t78 - pkin(4)) * t82;
t13 = -t78 * t145 + t18;
t101 = t86 * t12 - t84 * t13;
t100 = t84 * t12 + t86 * t13;
t99 = -t31 * t79 + t32 * t82;
t38 = t79 * t91;
t95 = -t132 + t133;
t44 = t95 * t79;
t98 = t38 * t62 + t44 * t61;
t97 = t41 * t79 - t42 * t82;
t40 = qJD(4) - t41;
t26 = qJDD(4) - t31;
t92 = t95 * t62;
t90 = t26 * t79 + t52 * t74 + t107;
t76 = pkin(9) + qJ(5);
t71 = cos(t76);
t69 = sin(t76);
t50 = -t78 * t126 + t81 * t127;
t49 = -t81 * t126 - t78 * t127;
t47 = (pkin(4) * t78 + t65) * t79;
t30 = t71 * t136 + t70 * t69;
t29 = -t69 * t136 + t70 * t71;
t28 = -t71 * t137 + t72 * t69;
t27 = t69 * t137 + t72 * t71;
t21 = pkin(4) * t116 + t40;
t17 = -t78 * t138 + t46;
t16 = pkin(4) * t113 + t26;
t1 = [qJDD(1), t105, g(1) * t87 + g(2) * t85, (t105 + (t80 ^ 2 + t83 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), (-t123 + t148) * t82, t52 * t130 + t107 + t99, t56 * t68 - g(1) * (-t70 * pkin(2) + t112) - g(2) * t117 + t99 * t65 - t97 * qJD(3), (-t49 * qJD(1) - t17 * qJDD(1) - t111 * t81 - t5) * t82 + t90 * t78, (t50 * qJD(1) + t18 * qJDD(1) + t111 * t78 + t6) * t82 + t90 * t81, t6 * t18 + t9 * t50 + t5 * t17 + t8 * t49 - g(1) * t112 - g(2) * (pkin(3) * t136 + t117) + (-qJ(4) * t67 + t40 * qJD(3) + t26 * t65) * t79 - t93 * t144, -t10 * t44 - t37 * t38, -t10 * t43 + t44 * t11 + t38 * t35 - t37 * t39, t98 - t141, t131 + t135, t61 * t82, -(t86 * t49 - t84 * t50) * t62 - t101 * t61 - t114 * t82 + t35 * t128 + t47 * t11 + t16 * t43 + t21 * t39 - g(1) * t28 - g(2) * t30 + (t100 * t62 - t102 * t82) * qJD(5), (t84 * t49 + t86 * t50) * t62 + t100 * t61 + t104 * t82 + t37 * t128 + t47 * t10 - t16 * t44 - t21 * t38 - g(1) * t27 - g(2) * t29 + (t101 * t62 + t103 * t82) * qJD(5); 0, 0, 0, t125, 0, 0, t31 * t82 + t32 * t79 - g(3), 0, 0, -t26 * t82 - g(3) + (-t5 * t78 + t6 * t81) * t79, 0, 0, 0, 0, 0, t131 - t135, -t98 - t141; 0, 0, 0, 0, -t121, -t110, qJD(1) * t97 - t148, -t78 * t110 - t81 * t121, -t81 * t110 + t78 * t121, t5 * t81 + t6 * t78 + (-t40 * t79 + (t78 * t8 - t81 * t9) * t82) * qJD(1) + t111, 0, 0, 0, 0, 0, t95 * t61 + t62 * t91 + (-t62 * t82 * t96 - t79 * t35) * qJD(1), t96 * t61 - qJD(5) * t92 + (-t79 * t37 + t82 * t92) * qJD(1); 0, 0, 0, 0, 0, 0, 0, (-t81 * t134 + t122) * t79, (qJDD(1) * t81 + t78 * t134) * t79, qJDD(4) + t48 - t125 * t82 + ((t78 * t9 + t8 * t81) * qJD(1) + t107) * t79, 0, 0, 0, 0, 0, t11 - t139, t10 + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t35, -t35 ^ 2 + t37 ^ 2, t10 - t140, -t11 - t139, -t61, -g(1) * t29 + g(2) * t27 + t147 * t102 + t69 * t143 - t21 * t37 + t114, g(1) * t30 - g(2) * t28 - t147 * t103 + t71 * t143 + t21 * t35 - t104;];
tau_reg = t1;
