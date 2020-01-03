% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPR6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:56
% EndTime: 2019-12-31 17:04:58
% DurationCPUTime: 0.81s
% Computational Cost: add. (1000->170), mult. (2432->247), div. (0->0), fcn. (1785->10), ass. (0->107)
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t101 = sin(qJ(2));
t104 = cos(qJ(2));
t97 = sin(pkin(7));
t98 = cos(pkin(7));
t71 = -t97 * t101 + t98 * t104;
t61 = t71 * qJD(1);
t53 = t103 * t61;
t72 = t98 * t101 + t97 * t104;
t63 = t72 * qJD(1);
t27 = -t100 * t63 + t53;
t94 = qJD(2) + qJD(4);
t135 = t27 * t94;
t130 = qJD(4) * t100;
t62 = t72 * qJD(2);
t36 = -qJD(1) * t62 + t71 * qJDD(1);
t129 = qJD(1) * qJD(2);
t125 = t104 * t129;
t126 = t101 * t129;
t37 = t72 * qJDD(1) + t98 * t125 - t97 * t126;
t2 = qJD(4) * t53 + t100 * t36 + t103 * t37 - t63 * t130;
t150 = t2 - t135;
t116 = t100 * t61 + t103 * t63;
t149 = t116 * t27;
t136 = t116 * t94;
t3 = t116 * qJD(4) + t100 * t37 - t103 * t36;
t148 = -t3 + t136;
t102 = sin(qJ(1));
t105 = cos(qJ(1));
t121 = g(1) * t105 + g(2) * t102;
t147 = t116 ^ 2 - t27 ^ 2;
t141 = t61 * pkin(6);
t133 = qJ(3) + pkin(5);
t84 = t133 * t104;
t77 = qJD(1) * t84;
t134 = t98 * t77;
t131 = qJD(2) * pkin(2);
t83 = t133 * t101;
t76 = qJD(1) * t83;
t70 = -t76 + t131;
t35 = t97 * t70 + t134;
t15 = t35 + t141;
t91 = t104 * pkin(2) + pkin(1);
t78 = -t91 * qJD(1) + qJD(3);
t43 = -t61 * pkin(3) + t78;
t92 = qJ(2) + pkin(7) + qJ(4);
t88 = sin(t92);
t89 = cos(t92);
t146 = g(3) * t88 + t121 * t89 + t15 * t130 - t43 * t27;
t144 = qJD(4) - t94;
t123 = qJD(2) * t133;
t59 = -t101 * qJD(3) - t104 * t123;
t33 = qJDD(2) * pkin(2) + t59 * qJD(1) - qJDD(1) * t83;
t58 = t104 * qJD(3) - t101 * t123;
t40 = t58 * qJD(1) + qJDD(1) * t84;
t6 = t98 * t33 - t97 * t40;
t4 = qJDD(2) * pkin(3) - t37 * pkin(6) + t6;
t7 = t97 * t33 + t98 * t40;
t5 = t36 * pkin(6) + t7;
t143 = -g(3) * t89 - t100 * t5 + t103 * t4 - t43 * t116 + t121 * t88;
t142 = pkin(2) * t97;
t140 = t63 * pkin(6);
t137 = g(3) * t104;
t66 = t97 * t77;
t19 = t98 * t58 + t97 * t59;
t42 = -t98 * t76 - t66;
t45 = -t97 * t83 + t98 * t84;
t95 = t101 ^ 2;
t132 = -t104 ^ 2 + t95;
t128 = t104 * qJDD(1);
t127 = t101 * t131;
t18 = -t97 * t58 + t98 * t59;
t34 = t98 * t70 - t66;
t41 = t97 * t76 - t134;
t44 = -t98 * t83 - t97 * t84;
t120 = g(1) * t102 - g(2) * t105;
t13 = qJD(2) * pkin(3) - t140 + t34;
t119 = -t100 * t13 - t103 * t15;
t20 = -t72 * pkin(6) + t44;
t21 = t71 * pkin(6) + t45;
t118 = -t100 * t21 + t103 * t20;
t117 = t100 * t20 + t103 * t21;
t38 = t100 * t72 - t103 * t71;
t39 = t100 * t71 + t103 * t72;
t90 = t98 * pkin(2) + pkin(3);
t115 = t100 * t90 + t103 * t142;
t114 = -t100 * t142 + t103 * t90;
t113 = -0.2e1 * pkin(1) * t129 - pkin(5) * qJDD(2);
t111 = pkin(2) * t126 - t91 * qJDD(1) + qJDD(3);
t106 = qJD(2) ^ 2;
t110 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t106 + t120;
t107 = qJD(1) ^ 2;
t109 = pkin(1) * t107 - pkin(5) * qJDD(1) + t121;
t93 = qJDD(2) + qJDD(4);
t65 = t71 * qJD(2);
t49 = -t71 * pkin(3) - t91;
t47 = t62 * pkin(3) + t127;
t46 = t101 * qJD(1) * pkin(2) + t63 * pkin(3);
t17 = t42 - t140;
t16 = t41 - t141;
t14 = -t36 * pkin(3) + t111;
t11 = -t62 * pkin(6) + t19;
t10 = -t65 * pkin(6) + t18;
t9 = t39 * qJD(4) + t100 * t65 + t103 * t62;
t8 = -t38 * qJD(4) - t100 * t62 + t103 * t65;
t1 = [qJDD(1), t120, t121, t95 * qJDD(1) + 0.2e1 * t101 * t125, 0.2e1 * t101 * t128 - 0.2e1 * t132 * t129, qJDD(2) * t101 + t106 * t104, qJDD(2) * t104 - t106 * t101, 0, t113 * t101 + t110 * t104, -t110 * t101 + t113 * t104, -t18 * t63 + t19 * t61 - t34 * t65 - t35 * t62 + t45 * t36 - t44 * t37 - t6 * t72 + t7 * t71 - t121, t7 * t45 + t35 * t19 + t6 * t44 + t34 * t18 - t111 * t91 + t78 * t127 - g(1) * (-t102 * t91 + t105 * t133) - g(2) * (t102 * t133 + t105 * t91), t116 * t8 + t2 * t39, -t116 * t9 - t2 * t38 + t27 * t8 - t39 * t3, t39 * t93 + t8 * t94, -t38 * t93 - t9 * t94, 0, -t47 * t27 + t49 * t3 + t14 * t38 + t43 * t9 + (-t117 * qJD(4) + t103 * t10 - t100 * t11) * t94 + t118 * t93 + t120 * t89, t47 * t116 + t49 * t2 + t14 * t39 + t43 * t8 - (t118 * qJD(4) + t100 * t10 + t103 * t11) * t94 - t117 * t93 - t120 * t88; 0, 0, 0, -t101 * t107 * t104, t132 * t107, t101 * qJDD(1), t128, qJDD(2), t109 * t101 - t137, g(3) * t101 + t109 * t104, (t35 + t41) * t63 + (t34 - t42) * t61 + (t36 * t97 - t37 * t98) * pkin(2), -t34 * t41 - t35 * t42 + (-t137 + t6 * t98 + t7 * t97 + (-qJD(1) * t78 + t121) * t101) * pkin(2), -t149, t147, t150, t148, t93, t114 * t93 + t46 * t27 - (-t100 * t17 + t103 * t16) * t94 + (-t115 * t94 + t119) * qJD(4) + t143, -t115 * t93 - t103 * t5 - t100 * t4 - t46 * t116 + (t100 * t16 + t103 * t17) * t94 + (-t103 * t13 - t114 * t94) * qJD(4) + t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61 ^ 2 - t63 ^ 2, t34 * t63 - t35 * t61 + t111 - t120, 0, 0, 0, 0, 0, t3 + t136, t2 + t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t147, t150, t148, t93, t144 * t119 + t143, (-t15 * t94 - t4) * t100 + (-t144 * t13 - t5) * t103 + t146;];
tau_reg = t1;
