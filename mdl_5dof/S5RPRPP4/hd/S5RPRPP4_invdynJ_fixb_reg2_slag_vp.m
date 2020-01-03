% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:55
% EndTime: 2019-12-31 18:14:57
% DurationCPUTime: 1.21s
% Computational Cost: add. (1630->256), mult. (3112->283), div. (0->0), fcn. (1956->8), ass. (0->139)
t109 = sin(qJ(1));
t111 = cos(qJ(1));
t197 = g(1) * t109 - g(2) * t111;
t110 = cos(qJ(3));
t158 = qJD(1) * qJD(3);
t147 = t110 * t158;
t108 = sin(qJ(3));
t155 = t108 * qJDD(1);
t199 = t147 + t155;
t114 = qJD(1) ^ 2;
t198 = -t114 * qJ(2) - t197;
t100 = qJDD(1) * qJ(2);
t103 = t108 ^ 2;
t104 = t110 ^ 2;
t164 = t103 + t104;
t112 = -pkin(1) - pkin(6);
t67 = t112 * qJDD(1) + qJDD(2);
t146 = t164 * t67;
t105 = sin(pkin(7));
t106 = cos(pkin(7));
t59 = t105 * t110 + t106 * t108;
t55 = t59 * qJD(3);
t60 = -t105 * t108 + t106 * t110;
t32 = -t55 * qJD(3) + t60 * qJDD(3);
t50 = t59 * qJD(1);
t196 = -qJD(1) * t50 + t32;
t148 = t108 * t158;
t154 = t110 * qJDD(1);
t153 = t105 * t154 + t199 * t106;
t30 = t105 * t148 - t153;
t159 = qJD(3) * t110;
t160 = qJD(3) * t108;
t52 = t105 * t160 - t106 * t159;
t133 = -t59 * t30 - t50 * t52;
t162 = qJD(1) * t108;
t150 = t105 * t162;
t161 = qJD(1) * t110;
t53 = t106 * t161 - t150;
t187 = t53 ^ 2;
t47 = t50 ^ 2;
t195 = -t47 - t187;
t194 = -t47 + t187;
t132 = t105 * t155 - t106 * t154;
t31 = qJD(1) * t55 + t132;
t193 = -t30 * pkin(4) + t31 * qJ(5);
t136 = g(1) * t111 + g(2) * t109;
t101 = qJD(1) * qJD(2);
t152 = 0.2e1 * t101;
t192 = 0.2e1 * t100 + t152 - t136;
t167 = qJ(4) - t112;
t123 = -t110 * qJD(4) + t167 * t160;
t145 = t167 * t110;
t42 = -qJD(3) * t145 - t108 * qJD(4);
t23 = t105 * t42 - t106 * t123;
t24 = t105 * t123 + t106 * t42;
t64 = t167 * t108;
t34 = -t105 * t64 + t106 * t145;
t35 = -t105 * t145 - t106 * t64;
t191 = t23 * t53 - t24 * t50 + t35 * t30 - t34 * t31;
t97 = qJ(3) + pkin(7);
t89 = cos(t97);
t190 = t24 * qJD(3) + t35 * qJDD(3) + t136 * t89;
t157 = qJD(1) * qJD(4);
t61 = t110 * t67;
t68 = t112 * qJD(1) + qJD(2);
t19 = -t110 * t157 - t68 * t160 + qJDD(3) * pkin(3) + t61 + (t148 - t154) * qJ(4);
t144 = -qJ(4) * qJD(1) + t68;
t25 = t144 * t159 + (-qJ(4) * qJDD(1) - t157 + t67) * t108;
t5 = t105 * t19 + t106 * t25;
t88 = sin(t97);
t189 = g(3) * t89 + t197 * t88 - t5;
t188 = (t53 + t150) * qJD(3) - t153;
t184 = g(3) * t108;
t93 = t108 * pkin(3);
t180 = t53 * t50;
t4 = -t105 * t25 + t106 * t19;
t45 = t144 * t108;
t39 = t106 * t45;
t46 = -qJ(4) * t161 + t110 * t68;
t41 = qJD(3) * pkin(3) + t46;
t21 = t105 * t41 + t39;
t179 = -t52 * qJD(3) + t59 * qJDD(3);
t178 = (t152 + t100) * qJ(2);
t177 = t111 * pkin(1) + t109 * qJ(2);
t175 = t105 * t45;
t82 = qJ(2) + t93;
t173 = pkin(1) * qJDD(1);
t171 = qJDD(3) * pkin(4);
t26 = t105 * t46 + t39;
t169 = t26 * qJD(3);
t63 = pkin(3) * t162 + qJD(1) * qJ(2) + qJD(4);
t168 = t63 * qJD(1);
t72 = pkin(3) * t159 + qJD(2);
t27 = t106 * t46 - t175;
t166 = qJD(5) - t27;
t165 = t103 - t104;
t113 = qJD(3) ^ 2;
t163 = -t113 - t114;
t156 = qJDD(3) * t108;
t151 = t110 * t114 * t108;
t92 = t111 * qJ(2);
t149 = -t109 * pkin(1) + t92;
t143 = t164 * qJDD(1);
t142 = qJDD(2) - t173;
t140 = t108 * t147;
t138 = qJD(1) * t53 + t179;
t137 = t88 * pkin(4) - t89 * qJ(5);
t6 = -t60 * t31 - t55 * t53;
t20 = t106 * t41 - t175;
t107 = -qJ(4) - pkin(6);
t131 = t109 * t107 + t111 * t93 + t149;
t38 = t199 * pkin(3) + qJDD(4) + t100 + t101;
t130 = -t111 * t107 + t109 * t93 + t177;
t128 = -t6 - t133;
t127 = g(3) * t88 - t197 * t89 + t4;
t126 = 0.2e1 * qJ(2) * t158 + qJDD(3) * t112;
t122 = t60 * t30 + t31 * t59 + t55 * t50 + t53 * t52;
t121 = (t53 - t150) * qJD(3) + t153;
t120 = -t136 + t38;
t22 = t50 * pkin(4) - t53 * qJ(5) + t63;
t119 = -t22 * t53 - qJDD(5) + t127;
t14 = -qJD(3) * pkin(4) + qJD(5) - t20;
t15 = qJD(3) * qJ(5) + t21;
t98 = qJDD(3) * qJ(5);
t2 = qJD(3) * qJD(5) + t5 + t98;
t3 = qJDD(5) - t171 - t4;
t118 = t14 * t55 - t15 * t52 + t2 * t59 - t3 * t60 - t197;
t117 = -t20 * t55 - t21 * t52 + t4 * t60 + t5 * t59 - t197;
t116 = -t23 * qJD(3) - t34 * qJDD(3) - t136 * t88;
t115 = -t112 * t113 + t192;
t90 = qJDD(3) * t110;
t81 = -t106 * pkin(3) - pkin(4);
t78 = t105 * pkin(3) + qJ(5);
t29 = t59 * pkin(4) - t60 * qJ(5) + t82;
t28 = pkin(3) * t161 + t53 * pkin(4) + t50 * qJ(5);
t10 = 0.2e1 * t50 * qJD(3) + t132;
t9 = -t52 * pkin(4) + t55 * qJ(5) - t60 * qJD(5) + t72;
t1 = -t53 * qJD(5) + t193 + t38;
t7 = [0, 0, 0, 0, 0, qJDD(1), t197, t136, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - 0.2e1 * t173 - t197, t192, -t142 * pkin(1) - g(1) * t149 - g(2) * t177 + t178, t104 * qJDD(1) - 0.2e1 * t140, -0.2e1 * t108 * t154 + 0.2e1 * t165 * t158, -t113 * t108 + t90, t103 * qJDD(1) + 0.2e1 * t140, -t113 * t110 - t156, 0, t115 * t108 + t126 * t110, -t126 * t108 + t115 * t110, -t112 * t143 - t146 + t197, -g(1) * (t112 * t109 + t92) - g(2) * (t111 * pkin(6) + t177) + t112 * t146 + t178, t6, t122, t32, t133, -t179, 0, -t82 * t30 + t38 * t59 + t72 * t50 - t63 * t52 + t116, -t82 * t31 + t38 * t60 + t72 * t53 - t63 * t55 - t190, -t117 + t191, -g(1) * t131 - g(2) * t130 - t20 * t23 + t21 * t24 - t4 * t34 + t5 * t35 + t38 * t82 + t63 * t72, t6, t32, -t122, 0, t179, t133, t1 * t59 - t22 * t52 - t29 * t30 + t9 * t50 + t116, -t118 + t191, -t1 * t60 + t22 * t55 + t29 * t31 - t9 * t53 + t190, t2 * t35 + t15 * t24 + t1 * t29 + t22 * t9 + t3 * t34 + t14 * t23 - g(1) * (t111 * t137 + t131) - g(2) * (t109 * t137 + t130); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t114, t198 + t142, 0, 0, 0, 0, 0, 0, t163 * t108 + t90, t163 * t110 - t156, -t143, t146 + t198, 0, 0, 0, 0, 0, 0, t196, -t138, t128, t117 - t168, 0, 0, 0, 0, 0, 0, t196, t128, t138, -t22 * qJD(1) + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, -t165 * t114, t154, -t151, -t155, qJDD(3), t110 * t198 + t184 + t61, g(3) * t110 + (-t198 - t67) * t108, 0, 0, t180, t194, -t132, -t180, t188, qJDD(3), t169 - t63 * t53 + (qJDD(3) * t106 - t50 * t161) * pkin(3) + t127, t27 * qJD(3) + t63 * t50 + (-qJDD(3) * t105 - t53 * t161) * pkin(3) + t189, (t21 - t26) * t53 + (-t20 + t27) * t50 + (t105 * t30 + t106 * t31) * pkin(3), t20 * t26 - t21 * t27 + (t184 + t105 * t5 + t106 * t4 + (-t197 - t168) * t110) * pkin(3), t180, -t132, -t194, qJDD(3), -t188, -t180, t169 - t28 * t50 + (pkin(4) - t81) * qJDD(3) + t119, t78 * t30 - t81 * t31 + (t15 - t26) * t53 + (t14 - t166) * t50, t78 * qJDD(3) - t22 * t50 + t28 * t53 + t98 + (0.2e1 * qJD(5) - t27) * qJD(3) - t189, t2 * t78 + t3 * t81 - t22 * t28 - t14 * t26 - g(3) * (-t137 - t93) + t166 * t15 - t197 * (pkin(3) * t110 + pkin(4) * t89 + qJ(5) * t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, -t10, t195, t20 * t53 + t21 * t50 + t120, 0, 0, 0, 0, 0, 0, t121, t195, t10, t15 * t50 + (-qJD(5) - t14) * t53 + t120 + t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) + t180, -t132, -t187 - t113, -t15 * qJD(3) - t119 - t171;];
tau_reg = t7;
