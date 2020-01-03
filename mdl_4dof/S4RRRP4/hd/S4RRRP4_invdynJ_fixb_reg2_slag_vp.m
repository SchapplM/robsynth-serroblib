% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:44
% EndTime: 2019-12-31 17:15:47
% DurationCPUTime: 1.33s
% Computational Cost: add. (1609->248), mult. (3928->306), div. (0->0), fcn. (2617->8), ass. (0->141)
t109 = sin(qJ(3));
t110 = sin(qJ(2));
t112 = cos(qJ(2));
t180 = cos(qJ(3));
t65 = t109 * t112 + t180 * t110;
t54 = t65 * qJD(1);
t168 = t54 * qJ(4);
t181 = pkin(6) + pkin(5);
t77 = t181 * t112;
t72 = qJD(1) * t77;
t56 = t109 * t72;
t173 = qJD(2) * pkin(2);
t76 = t181 * t110;
t70 = qJD(1) * t76;
t62 = -t70 + t173;
t29 = t180 * t62 - t56;
t19 = t29 - t168;
t105 = qJD(2) + qJD(3);
t103 = qJDD(2) + qJDD(3);
t144 = t180 * qJD(3);
t137 = pkin(2) * t144;
t179 = pkin(2) * t109;
t185 = -t103 * t179 - t105 * t137;
t100 = t103 * pkin(3);
t161 = t109 * t110;
t129 = t105 * t161;
t146 = t180 * t112;
t134 = qJD(1) * t146;
t139 = qJDD(1) * t180;
t153 = t112 * qJDD(1);
t138 = -t105 * t134 - t109 * t153 - t110 * t139;
t21 = qJD(1) * t129 + t138;
t172 = t21 * qJ(4);
t184 = t100 + t172;
t108 = qJ(2) + qJ(3);
t101 = sin(t108);
t113 = cos(qJ(1));
t164 = t101 * t113;
t111 = sin(qJ(1));
t165 = t101 * t111;
t102 = cos(t108);
t178 = g(3) * t102;
t183 = g(1) * t164 + g(2) * t165 - t178;
t182 = t54 ^ 2;
t177 = g(3) * t112;
t157 = qJD(1) * t110;
t145 = t109 * t157;
t52 = -t134 + t145;
t176 = t54 * t52;
t15 = t105 * pkin(3) + t19;
t175 = t15 - t19;
t154 = t110 * qJDD(1);
t130 = t109 * t154 - t112 * t139;
t36 = t105 * t65;
t22 = t36 * qJD(1) + t130;
t174 = -t52 * t137 - t22 * t179;
t34 = -t180 * t70 - t56;
t42 = -t109 * t76 + t180 * t77;
t171 = t22 * qJ(4);
t170 = t52 * qJ(4);
t169 = t52 * t105;
t166 = pkin(5) * qJDD(1);
t163 = t102 * t111;
t162 = t102 * t113;
t141 = t52 * pkin(3) + qJD(4);
t99 = t112 * pkin(2) + pkin(1);
t75 = t99 * qJD(1);
t39 = t141 - t75;
t160 = qJD(4) + t39;
t106 = t110 ^ 2;
t107 = t112 ^ 2;
t159 = t106 - t107;
t158 = t106 + t107;
t156 = qJD(3) * t109;
t155 = qJD(1) * qJD(2);
t152 = t180 * pkin(2);
t151 = pkin(2) * t156;
t150 = t110 * t173;
t149 = pkin(2) * t157;
t60 = t180 * t72;
t116 = qJD(1) ^ 2;
t148 = t110 * t116 * t112;
t147 = qJD(2) * t181;
t143 = t110 * t155;
t142 = t112 * t155;
t38 = qJDD(2) * pkin(2) + t181 * (-t142 - t154);
t40 = t181 * (-t143 + t153);
t140 = -t109 * t40 + t180 * t38;
t33 = t109 * t70 - t60;
t41 = -t109 * t77 - t180 * t76;
t8 = t109 * t38 + t62 * t144 - t72 * t156 + t180 * t40;
t136 = -g(1) * t165 + g(2) * t164;
t135 = g(1) * t163 - g(2) * t162;
t133 = t110 * t142;
t132 = g(1) * t113 + g(2) * t111;
t131 = g(1) * t111 - g(2) * t113;
t30 = t109 * t62 + t60;
t128 = g(1) * t162 + g(2) * t163 + g(3) * t101 - t8;
t127 = -0.2e1 * pkin(1) * t155 - pkin(5) * qJDD(2);
t71 = t110 * t147;
t73 = t112 * t147;
t13 = -t109 * t73 - t76 * t144 - t77 * t156 - t180 * t71;
t49 = pkin(2) * t143 - t99 * qJDD(1);
t124 = -t75 * t52 + t128;
t115 = qJD(2) ^ 2;
t123 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t115 + t131;
t122 = pkin(1) * t116 + t132 - t166;
t10 = t22 * pkin(3) + qJDD(4) + t49;
t121 = -t105 * t145 - t138;
t9 = -t30 * qJD(3) + t140;
t14 = -t42 * qJD(3) + t109 * t71 - t180 * t73;
t120 = t160 * t52 + t128 + t171;
t119 = t9 + t183;
t118 = t75 * t54 + t119;
t104 = -qJ(4) - t181;
t98 = t152 + pkin(3);
t96 = pkin(3) * t102;
t69 = t96 + t99;
t64 = -t146 + t161;
t51 = t52 ^ 2;
t45 = t64 * pkin(3) - t99;
t43 = t54 * pkin(3) + t149;
t35 = -qJD(2) * t146 - t112 * t144 + t129;
t28 = t36 * pkin(3) + t150;
t27 = -t64 * qJ(4) + t42;
t26 = -t65 * qJ(4) + t41;
t25 = -t51 + t182;
t24 = -t168 + t34;
t23 = t33 + t170;
t20 = t30 - t170;
t18 = -t64 * t103 - t36 * t105;
t17 = t65 * t103 - t35 * t105;
t11 = t121 + t169;
t7 = t35 * qJ(4) - t65 * qJD(4) + t14;
t6 = -t36 * qJ(4) - t64 * qJD(4) + t13;
t5 = t22 * t64 + t52 * t36;
t4 = -t21 * t65 - t54 * t35;
t3 = -t52 * qJD(4) - t171 + t8;
t2 = -t54 * qJD(4) + t184 + t9;
t1 = t21 * t64 - t65 * t22 + t35 * t52 - t54 * t36;
t12 = [0, 0, 0, 0, 0, qJDD(1), t131, t132, 0, 0, t106 * qJDD(1) + 0.2e1 * t133, 0.2e1 * t110 * t153 - 0.2e1 * t159 * t155, qJDD(2) * t110 + t115 * t112, t107 * qJDD(1) - 0.2e1 * t133, qJDD(2) * t112 - t115 * t110, 0, t127 * t110 + t123 * t112, -t123 * t110 + t127 * t112, 0.2e1 * t158 * t166 - t132, -g(1) * (-t111 * pkin(1) + t113 * pkin(5)) - g(2) * (t113 * pkin(1) + t111 * pkin(5)) + (t158 * pkin(5) ^ 2 + pkin(1) ^ 2) * qJDD(1), t4, t1, t17, t5, t18, 0, t41 * t103 + t14 * t105 + t52 * t150 - t99 * t22 - t75 * t36 + t49 * t64 + t135, -t42 * t103 - t13 * t105 + t54 * t150 + t99 * t21 + t75 * t35 + t49 * t65 + t136, -t13 * t52 - t14 * t54 + t41 * t21 - t42 * t22 + t29 * t35 - t30 * t36 - t8 * t64 - t9 * t65 - t132, t8 * t42 + t30 * t13 + t9 * t41 + t29 * t14 - t49 * t99 - t75 * t150 - g(1) * (-t111 * t99 + t113 * t181) - g(2) * (t111 * t181 + t113 * t99), t4, t1, t17, t5, t18, 0, t10 * t64 + t26 * t103 + t7 * t105 + t45 * t22 + t28 * t52 + t39 * t36 + t135, t10 * t65 - t27 * t103 - t6 * t105 - t45 * t21 + t28 * t54 - t39 * t35 + t136, t15 * t35 - t2 * t65 - t20 * t36 + t26 * t21 - t27 * t22 - t3 * t64 - t6 * t52 - t7 * t54 - t132, t3 * t27 + t20 * t6 + t2 * t26 + t15 * t7 + t10 * t45 + t39 * t28 - g(1) * (-t113 * t104 - t111 * t69) - g(2) * (-t111 * t104 + t113 * t69); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t159 * t116, t154, t148, t153, qJDD(2), t122 * t110 - t177, g(3) * t110 + t122 * t112, 0, 0, t176, t25, t11, -t176, -t130, t103, -t33 * t105 + (t180 * t103 - t105 * t156 - t52 * t157) * pkin(2) + t118, t34 * t105 - t54 * t149 + t124 + t185, t21 * t152 + (-t29 + t34) * t52 + (t30 + t33 + t151) * t54 + t174, -t29 * t33 - t30 * t34 + (t180 * t9 - t177 + t109 * t8 + (-t109 * t29 + t180 * t30) * qJD(3) + (qJD(1) * t75 + t132) * t110) * pkin(2), t176, t25, t11, -t176, -t130, t103, t98 * t103 - t23 * t105 - t43 * t52 - t160 * t54 + (-t60 + (-pkin(2) * t105 - t62) * t109) * qJD(3) + t140 + t183 + t184, t24 * t105 - t43 * t54 + t120 + t185, t98 * t21 + (-t15 + t24) * t52 + (t20 + t23 + t151) * t54 + t174, -g(3) * t96 - t15 * t23 + t2 * t98 - t20 * t24 - t39 * t43 - t132 * (-t110 * pkin(2) - pkin(3) * t101) + (-t177 + t3 * t109 + (-t109 * t15 + t180 * t20) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t25, t11, -t176, -t130, t103, t30 * t105 + t118, t29 * t105 + t124, 0, 0, t176, t25, t11, -t176, -t130, t103, t172 + t20 * t105 + 0.2e1 * t100 + (-t141 - t39) * t54 + t119, -t182 * pkin(3) + t19 * t105 + t120, t21 * pkin(3) - t175 * t52, t175 * t20 + (t101 * t132 - t39 * t54 - t178 + t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t105 + t22, t121 - t169, -t51 - t182, t15 * t54 + t20 * t52 + t10 - t131;];
tau_reg = t12;
