% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:43
% EndTime: 2019-12-31 17:15:45
% DurationCPUTime: 0.77s
% Computational Cost: add. (857->163), mult. (2060->224), div. (0->0), fcn. (1362->8), ass. (0->105)
t136 = cos(qJ(3));
t142 = pkin(5) + pkin(6);
t84 = cos(qJ(2));
t53 = t142 * t84;
t48 = qJD(1) * t53;
t81 = sin(qJ(3));
t32 = t81 * t48;
t126 = qJD(2) * pkin(2);
t82 = sin(qJ(2));
t52 = t142 * t82;
t46 = qJD(1) * t52;
t38 = -t46 + t126;
t111 = t136 * t38 - t32;
t41 = t136 * t82 + t81 * t84;
t30 = t41 * qJD(1);
t124 = t30 * qJ(4);
t146 = t124 - t111;
t83 = sin(qJ(1));
t85 = cos(qJ(1));
t104 = g(1) * t85 + g(2) * t83;
t80 = qJ(2) + qJ(3);
t72 = sin(t80);
t73 = cos(t80);
t145 = -g(3) * t73 + t104 * t72;
t77 = qJD(2) + qJD(3);
t74 = t84 * pkin(2);
t137 = pkin(1) + t74;
t129 = t136 * t53 - t81 * t52;
t143 = t30 ^ 2;
t7 = t77 * pkin(3) - t146;
t141 = t7 + t146;
t115 = t136 * t84;
t105 = qJD(1) * t115;
t123 = qJD(1) * t82;
t117 = t81 * t123;
t28 = -t105 + t117;
t51 = t137 * qJD(1);
t21 = t28 * pkin(3) + qJD(4) - t51;
t135 = t21 * t30;
t134 = t30 * t28;
t133 = t73 * t83;
t132 = t73 * t85;
t131 = t81 * t82;
t130 = -t136 * t46 - t32;
t128 = pkin(3) * t73 + t74;
t78 = t82 ^ 2;
t127 = -t84 ^ 2 + t78;
t125 = t28 * qJ(4);
t122 = qJD(3) * t81;
t121 = t82 * qJDD(1);
t120 = t84 * qJDD(1);
t119 = qJD(1) * qJD(2);
t118 = t82 * t126;
t36 = t136 * t48;
t114 = qJD(2) * t142;
t113 = t82 * t119;
t112 = t84 * t119;
t110 = t81 * t46 - t36;
t109 = -t136 * t52 - t81 * t53;
t108 = t136 * qJD(3);
t107 = qJDD(1) * t136;
t106 = -t77 * t105 - t82 * t107 - t81 * t120;
t103 = g(1) * t83 - g(2) * t85;
t102 = -t84 * t107 + t81 * t121;
t101 = t77 * t131;
t100 = -t81 * t38 - t36;
t47 = t82 * t114;
t49 = t84 * t114;
t99 = -t52 * t108 - t53 * t122 - t136 * t47 - t81 * t49;
t25 = pkin(2) * t113 - qJDD(1) * t137;
t98 = -0.2e1 * pkin(1) * t119 - pkin(5) * qJDD(2);
t86 = qJD(2) ^ 2;
t95 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t86 + t103;
t87 = qJD(1) ^ 2;
t94 = pkin(1) * t87 - pkin(5) * qJDD(1) + t104;
t19 = t77 * t41;
t11 = t19 * qJD(1) + t102;
t93 = t11 * pkin(3) + qJDD(4) + t25;
t20 = qJDD(2) * pkin(2) + t142 * (-t112 - t121);
t22 = t142 * (-t113 + t120);
t92 = t100 * qJD(3) + t136 * t20 - t81 * t22;
t91 = -t129 * qJD(3) - t136 * t49 + t81 * t47;
t90 = t38 * t108 - t48 * t122 + t136 * t22 + t81 * t20;
t89 = g(1) * t132 + g(2) * t133 + g(3) * t72 - t51 * t28 - t90;
t88 = t51 * t30 + t145 + t92;
t76 = -qJ(4) - t142;
t75 = qJDD(2) + qJDD(3);
t70 = t136 * pkin(2) + pkin(3);
t45 = pkin(1) + t128;
t40 = -t115 + t131;
t27 = t28 ^ 2;
t18 = -qJD(2) * t115 - t84 * t108 + t101;
t16 = -t40 * qJ(4) + t129;
t15 = -t41 * qJ(4) + t109;
t14 = -t27 + t143;
t13 = -t124 + t130;
t12 = t110 + t125;
t10 = qJD(1) * t101 + t106;
t9 = -t100 - t125;
t5 = -t106 + (-t117 + t28) * t77;
t4 = t18 * qJ(4) - t41 * qJD(4) + t91;
t3 = -t19 * qJ(4) - t40 * qJD(4) + t99;
t2 = -t11 * qJ(4) - t28 * qJD(4) + t90;
t1 = t75 * pkin(3) + t10 * qJ(4) - t30 * qJD(4) + t92;
t6 = [qJDD(1), t103, t104, t78 * qJDD(1) + 0.2e1 * t82 * t112, -0.2e1 * t127 * t119 + 0.2e1 * t82 * t120, qJDD(2) * t82 + t86 * t84, qJDD(2) * t84 - t86 * t82, 0, t98 * t82 + t95 * t84, -t95 * t82 + t98 * t84, -t10 * t41 - t30 * t18, t10 * t40 - t41 * t11 + t18 * t28 - t30 * t19, -t18 * t77 + t41 * t75, -t19 * t77 - t40 * t75, 0, g(1) * t133 - g(2) * t132 + t109 * t75 - t11 * t137 + t28 * t118 - t51 * t19 + t25 * t40 + t91 * t77, t10 * t137 - t103 * t72 + t30 * t118 - t129 * t75 + t51 * t18 + t25 * t41 - t99 * t77, -t1 * t41 + t15 * t10 - t16 * t11 + t7 * t18 - t9 * t19 - t2 * t40 - t3 * t28 - t4 * t30 - t104, t2 * t16 + t9 * t3 + t1 * t15 + t7 * t4 + t93 * (t40 * pkin(3) - t137) + t21 * (t19 * pkin(3) + t118) - g(1) * (-t83 * t45 - t85 * t76) - g(2) * (t85 * t45 - t83 * t76); 0, 0, 0, -t82 * t87 * t84, t127 * t87, t121, t120, qJDD(2), -g(3) * t84 + t94 * t82, g(3) * t82 + t94 * t84, t134, t14, t5, -t102, t75, -t110 * t77 + (-t77 * t122 - t28 * t123 + t136 * t75) * pkin(2) + t88, t130 * t77 + (-t77 * t108 - t30 * t123 - t81 * t75) * pkin(2) + t89, t70 * t10 + (t12 + t9) * t30 + (t13 - t7) * t28 + (-t11 * t81 + (-t136 * t28 + t30 * t81) * qJD(3)) * pkin(2), t1 * t70 - t9 * t13 - t7 * t12 - pkin(3) * t135 - g(3) * t128 - t104 * (-t82 * pkin(2) - pkin(3) * t72) + (-t21 * t123 + t2 * t81 + (t136 * t9 - t7 * t81) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t14, t5, -t102, t75, -t100 * t77 + t88, t111 * t77 + t89, pkin(3) * t10 - t141 * t28, t141 * t9 + (t1 - t135 + t145) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27 - t143, t9 * t28 + t7 * t30 - t103 + t93;];
tau_reg = t6;
