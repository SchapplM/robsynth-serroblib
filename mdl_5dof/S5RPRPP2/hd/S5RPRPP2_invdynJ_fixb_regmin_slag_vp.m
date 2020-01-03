% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:14
% EndTime: 2019-12-31 18:11:16
% DurationCPUTime: 0.95s
% Computational Cost: add. (864->214), mult. (1582->243), div. (0->0), fcn. (835->8), ass. (0->130)
t88 = cos(qJ(3));
t145 = qJ(4) * t88;
t158 = pkin(3) + pkin(4);
t86 = sin(qJ(3));
t161 = t158 * t86;
t103 = t145 - t161;
t85 = cos(pkin(7));
t153 = t85 * pkin(1);
t125 = pkin(2) + t153;
t35 = qJDD(1) * t125;
t70 = t86 * qJ(4);
t120 = pkin(2) + t70;
t130 = qJD(1) * qJD(3);
t121 = t88 * t130;
t65 = t86 * qJDD(1);
t165 = t121 + t65;
t84 = sin(pkin(7));
t53 = t84 * pkin(1) + pkin(6);
t34 = t53 * qJDD(1);
t164 = qJD(2) * qJD(3) + t34;
t36 = t53 * qJD(1);
t29 = t86 * t36;
t19 = t88 * qJD(2) - t29;
t163 = qJD(4) - t19;
t123 = t158 * qJD(3);
t134 = qJ(5) * qJD(1);
t10 = t86 * t134 + t19;
t135 = qJD(4) - t10;
t7 = -t123 + t135;
t37 = t125 * qJD(1);
t119 = t158 * qJDD(3);
t20 = t86 * qJD(2) + t88 * t36;
t11 = -t88 * t134 + t20;
t80 = qJD(3) * qJ(4);
t9 = t11 + t80;
t16 = t80 + t20;
t8 = qJD(5) + (t158 * t88 + t120 + t153) * qJD(1);
t142 = qJD(5) + t8;
t162 = t142 * t86;
t12 = -qJD(3) * pkin(3) + t163;
t77 = qJ(1) + pkin(7);
t62 = sin(t77);
t63 = cos(t77);
t110 = g(1) * t63 + g(2) * t62;
t150 = t63 * t86;
t151 = t62 * t86;
t160 = g(1) * t150 + g(2) * t151 - g(3) * t88;
t132 = qJDD(3) * t53;
t74 = t88 * pkin(3);
t104 = -t120 - t74;
t17 = (t104 - t153) * qJD(1);
t148 = t74 + t70;
t23 = -t125 - t148;
t159 = (qJD(1) * t23 + t17) * qJD(3) - t132;
t157 = g(1) * t62;
t154 = g(2) * t63;
t73 = t88 * pkin(4);
t152 = t88 * t9;
t149 = t63 * t88;
t81 = t86 ^ 2;
t82 = t88 ^ 2;
t147 = t81 - t82;
t146 = t81 + t82;
t144 = t9 * qJD(3);
t143 = qJ(5) - t53;
t141 = qJD(1) * t86;
t25 = t143 * t88;
t140 = qJD(3) * t25;
t139 = qJD(3) * t36;
t138 = qJD(3) * t86;
t137 = qJDD(3) * pkin(3);
t136 = t86 * qJD(4);
t133 = qJ(5) * qJD(3);
t67 = t88 * qJDD(1);
t131 = qJ(5) * qJDD(1);
t129 = qJD(1) * qJD(5);
t92 = qJD(1) ^ 2;
t127 = t86 * t92 * t88;
t126 = t86 * qJDD(2) + t164 * t88;
t87 = sin(qJ(1));
t124 = -t87 * pkin(1) + t63 * pkin(6);
t122 = t86 * t130;
t118 = -g(3) - t139;
t18 = -t23 + t73;
t117 = qJD(1) * t18 + t8;
t116 = (-qJDD(2) + t139) * t88 + t164 * t86;
t113 = 0.2e1 * t121;
t78 = qJDD(3) * qJ(4);
t79 = qJD(3) * qJD(4);
t112 = 0.2e1 * t78 + 0.2e1 * t79 + t126;
t89 = cos(qJ(1));
t111 = t89 * pkin(1) + pkin(3) * t149 + t62 * pkin(6) + t120 * t63;
t109 = g(1) * t87 - g(2) * t89;
t91 = qJD(3) ^ 2;
t108 = t53 * t91 + t154;
t107 = pkin(3) * t86 - t145;
t106 = -qJDD(4) - t116;
t105 = pkin(3) * t67 + t165 * qJ(4) + qJD(1) * t136 + t35;
t102 = -t116 + t160;
t101 = pkin(4) * t67 + qJDD(5) + t105;
t5 = -t106 - t137;
t4 = -t36 * t138 + t126 + t78 + t79;
t100 = t108 - 0.2e1 * t35;
t13 = t103 * qJD(3) + t136;
t3 = -t158 * t122 + t101;
t99 = qJD(1) * t13 + qJDD(1) * t18 - t154 + t3;
t98 = -0.2e1 * t37 * qJD(3) - t132;
t97 = t20 * qJD(3) + t102;
t22 = t107 * qJD(3) - t136;
t6 = pkin(3) * t122 - t105;
t96 = -qJD(1) * t22 - qJDD(1) * t23 - t108 - t6;
t95 = -t86 * t131 + qJDD(4) - t102;
t94 = t4 * t88 + t5 * t86 + (t12 * t88 - t16 * t86) * qJD(3);
t48 = qJ(5) * t122;
t47 = -t81 * t92 - t91;
t45 = t88 * t157;
t44 = g(1) * t151;
t41 = t63 * t145;
t39 = t62 * t145;
t38 = -qJDD(3) - t127;
t33 = qJDD(3) * t88 - t91 * t86;
t32 = qJDD(3) * t86 + t91 * t88;
t31 = t107 * qJD(1);
t24 = t143 * t86;
t21 = t103 * qJD(1);
t15 = -t86 * qJD(5) - t140;
t14 = -t88 * qJD(5) + t143 * t138;
t2 = t48 + (-t129 - t131) * t88 + t4;
t1 = -t165 * qJ(5) - t86 * t129 - t106 - t119;
t26 = [qJDD(1), t109, g(1) * t89 + g(2) * t87, (t109 + (t84 ^ 2 + t85 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t81 * qJDD(1) + t86 * t113, -0.2e1 * t147 * t130 + 0.2e1 * t86 * t67, t32, t33, 0, -t100 * t88 + t98 * t86 + t45, t100 * t86 + t98 * t88 - t44, t159 * t86 + t96 * t88 + t45, t146 * t34 - t110 + t94, -t159 * t88 + t96 * t86 + t44, -g(1) * t124 - g(2) * t111 - t104 * t157 + t17 * t22 + t6 * t23 + t94 * t53, t24 * qJDD(3) + t45 + (-t117 * t86 - t15) * qJD(3) + t99 * t88, -t25 * qJDD(3) + t44 + (t117 * t88 + t14) * qJD(3) + t99 * t86, (-qJD(3) * t7 + qJDD(1) * t25 - t2 + (qJD(3) * t24 - t14) * qJD(1)) * t88 + (t144 + qJDD(1) * t24 - t1 + (-t15 - t140) * qJD(1)) * t86 + t110, -t2 * t25 + t9 * t14 - t1 * t24 + t7 * t15 + t3 * t18 + t8 * t13 - g(1) * (-t63 * qJ(5) + t124) - g(2) * (pkin(4) * t149 + t111) + (-g(1) * (t104 - t73) + g(2) * qJ(5)) * t62; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, t33, -t32, t33, 0, t32, t4 * t86 - t5 * t88 - g(3) + (t12 * t86 + t16 * t88) * qJD(3), t33, t32, 0, -t1 * t88 + t2 * t86 - g(3) + (t7 * t86 + t152) * qJD(3); 0, 0, 0, 0, -t127, t147 * t92, t65, t67, qJDD(3), t37 * t141 + t97, g(3) * t86 + (t19 + t29) * qJD(3) + (qJD(1) * t37 + t110) * t88 - t126, 0.2e1 * t137 - qJDD(4) + (-t17 * t86 + t31 * t88) * qJD(1) + t97, -t107 * qJDD(1), -t19 * qJD(3) + (qJD(1) * t17 - t110) * t88 + (qJD(1) * t31 + t118) * t86 + t112, t4 * qJ(4) - t5 * pkin(3) - t17 * t31 - t12 * t20 - g(1) * (-pkin(3) * t150 + t41) - g(2) * (-pkin(3) * t151 + t39) - g(3) * t148 + t163 * t16, t11 * qJD(3) + 0.2e1 * t119 + ((-t21 + t133) * t88 + t162) * qJD(1) - t95, -t10 * qJD(3) + t48 + (-qJD(1) * t21 + t118) * t86 + (-t142 * qJD(1) - t110 - t131) * t88 + t112, -t103 * qJDD(1), t2 * qJ(4) - t1 * t158 - t7 * t11 - t8 * t21 - g(1) * t41 - g(2) * t39 - g(3) * (t73 + t148) + t135 * t9 + t110 * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t65, t47, -t16 * qJD(3) + t17 * t141 - t160 + t5, t38, t47, -t65, -t144 - t119 + (-t88 * t133 - t162) * qJD(1) + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 - 0.2e1 * t122, t65 + t113, -t146 * t92, t157 - t154 + (t152 + (t7 - t123) * t86) * qJD(1) + t101;];
tau_reg = t26;
