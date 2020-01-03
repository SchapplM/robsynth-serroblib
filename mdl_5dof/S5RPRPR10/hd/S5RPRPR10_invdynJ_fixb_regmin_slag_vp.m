% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:06
% EndTime: 2019-12-31 18:26:08
% DurationCPUTime: 0.51s
% Computational Cost: add. (896->137), mult. (1336->176), div. (0->0), fcn. (758->10), ass. (0->93)
t101 = qJD(1) - qJD(3);
t130 = t101 ^ 2;
t105 = qJ(2) * qJD(1);
t79 = -pkin(1) - pkin(2);
t129 = -qJD(3) * t105 + t79 * qJDD(1) + qJDD(2);
t71 = sin(pkin(8));
t72 = cos(pkin(8));
t74 = sin(qJ(3));
t77 = cos(qJ(3));
t36 = t71 * t74 - t72 * t77;
t115 = t101 * t36;
t104 = (qJ(2) * qJDD(1));
t103 = (qJD(1) * qJD(2));
t50 = t79 * qJD(1) + qJD(2);
t91 = qJD(3) * t50 + t103;
t128 = -t91 - t104;
t59 = t77 * t79;
t127 = -qJ(2) * t74 + t59;
t75 = sin(qJ(1));
t118 = t75 * t74;
t78 = cos(qJ(1));
t39 = -t78 * t77 - t118;
t119 = t74 * t78;
t40 = -t75 * t77 + t119;
t126 = g(1) * t40 - g(2) * t39;
t37 = t71 * t77 + t72 * t74;
t116 = t101 * t37;
t67 = qJDD(1) - qJDD(3);
t80 = qJD(5) ^ 2;
t125 = t101 * t116 - t36 * t67 + t37 * t80;
t14 = -t128 * t77 + t129 * t74;
t90 = t129 * t77;
t82 = t128 * t74 + t90;
t8 = -pkin(3) * t67 + t82;
t4 = t72 * t14 + t71 * t8;
t27 = t77 * t105 + t74 * t50;
t122 = t27 * t71;
t121 = t72 * t27;
t73 = sin(qJ(5));
t76 = cos(qJ(5));
t120 = t73 * t76;
t117 = t76 * t67;
t41 = -pkin(3) + t127;
t45 = qJ(2) * t77 + t74 * t79;
t114 = t71 * t41 + t72 * t45;
t113 = t78 * pkin(1) + t75 * qJ(2);
t112 = g(1) * t75 - g(2) * t78;
t69 = t73 ^ 2;
t111 = -t76 ^ 2 + t69;
t109 = pkin(1) * qJDD(1);
t107 = qJD(5) * t101;
t106 = qJDD(4) + g(3);
t102 = qJ(3) + pkin(8);
t100 = 2 * t103;
t98 = cos(t102);
t96 = qJDD(2) - t109;
t95 = g(1) * t78 + g(2) * t75;
t26 = -t74 * t105 + t77 * t50;
t3 = -t14 * t71 + t72 * t8;
t21 = -pkin(3) * t101 + t26;
t9 = t21 * t72 - t122;
t93 = t41 * t72 - t45 * t71;
t60 = sin(t102);
t28 = t60 * t78 - t75 * t98;
t29 = t75 * t60 + t78 * t98;
t92 = -pkin(4) * t67 + g(1) * t28 + g(2) * t29 + t3;
t5 = pkin(4) * t101 - t9;
t89 = pkin(7) * t67 + g(1) * t29 - g(2) * t28 + t101 * t5 - t4;
t24 = qJD(2) * t77 + t127 * qJD(3);
t25 = -t74 * qJD(2) - t45 * qJD(3);
t13 = t24 * t72 + t25 * t71;
t17 = pkin(4) - t93;
t18 = -pkin(7) + t114;
t88 = -qJDD(5) * t18 + (-t101 * t17 - t13 - t5) * qJD(5);
t16 = t26 * t72 - t122;
t51 = pkin(3) * t71 + pkin(7);
t52 = -pkin(3) * t72 - pkin(4);
t87 = -qJDD(5) * t51 + (-t101 * t52 + t16 + t5) * qJD(5);
t86 = -0.2e1 * qJD(5) * t115 - qJDD(5) * t37;
t12 = t24 * t71 - t25 * t72;
t85 = -t101 * t12 - t17 * t67 + t18 * t80 + t92;
t15 = t26 * t71 + t121;
t84 = -t101 * t15 - t51 * t80 + t52 * t67 + t92;
t83 = -g(1) * t39 - g(2) * t40 - t14;
t81 = qJD(1) ^ 2;
t62 = t78 * qJ(2);
t58 = pkin(3) * t77 + pkin(2);
t47 = qJDD(5) * t76 - t73 * t80;
t46 = qJDD(5) * t73 + t76 * t80;
t23 = -0.2e1 * t107 * t120 - t67 * t69;
t19 = t111 * t107 - t73 * t117;
t10 = t71 * t21 + t121;
t1 = [qJDD(1), t112, t95, -qJDD(2) + 0.2e1 * t109 + t112, t100 - t95 + (2 * t104), -t96 * pkin(1) - g(1) * (-pkin(1) * t75 + t62) - g(2) * t113 + (t100 + t104) * qJ(2), t67, -t25 * t101 - t59 * t67 + ((qJDD(1) + t67) * qJ(2) + t91) * t74 - t90 - t126, t101 * t24 + t45 * t67 - t83, t4 * t114 + t10 * t13 + t3 * t93 - t9 * t12 - g(1) * (pkin(3) * t119 + t62 + (-pkin(1) - t58) * t75) - g(2) * (pkin(3) * t118 + t58 * t78 + t113), -t23, -0.2e1 * t19, -t46, -t47, 0, t88 * t73 - t85 * t76, t85 * t73 + t88 * t76; 0, 0, 0, -qJDD(1), -t81, -qJ(2) * t81 - t112 + t96, 0, -t130 * t74 - t67 * t77, -t130 * t77 + t74 * t67, t115 * t10 + t116 * t9 - t3 * t36 + t37 * t4 - t112, 0, 0, 0, 0, 0, -t125 * t76 + t86 * t73, t125 * t73 + t86 * t76; 0, 0, 0, 0, 0, 0, -t67, -t101 * t27 + t126 + t82, -t101 * t26 + t83, -t10 * t16 + t9 * t15 + (t3 * t72 + t4 * t71 + t126) * pkin(3), t23, 0.2e1 * t19, t46, t47, 0, t87 * t73 + t84 * t76, -t84 * t73 + t87 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, 0, 0, 0, 0, 0, t47, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130 * t120, t111 * t130, -t73 * t67, -t117, qJDD(5), t106 * t76 + t89 * t73, -t106 * t73 + t89 * t76;];
tau_reg = t1;
