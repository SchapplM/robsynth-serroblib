% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:29:13
% EndTime: 2019-12-05 17:29:17
% DurationCPUTime: 1.73s
% Computational Cost: add. (1859->227), mult. (3932->329), div. (0->0), fcn. (2778->14), ass. (0->146)
t102 = sin(pkin(8));
t108 = sin(qJ(5));
t104 = cos(pkin(9));
t110 = cos(qJ(5));
t160 = t110 * t104;
t145 = t102 * t160;
t101 = sin(pkin(9));
t150 = qJDD(1) * t101;
t71 = t110 * t101 + t108 * t104;
t188 = t71 * qJD(5);
t19 = -qJDD(1) * t145 + (qJD(1) * t188 + t108 * t150) * t102;
t105 = cos(pkin(8));
t84 = t105 * qJD(1) - qJD(5);
t159 = qJD(1) * t102;
t144 = t101 * t159;
t135 = t108 * t144;
t143 = qJD(5) * t160;
t20 = -qJD(5) * t135 + t102 * (qJD(1) * t143 + t71 * qJDD(1));
t163 = pkin(1) * qJDD(1);
t103 = sin(pkin(7));
t86 = t103 * pkin(1) + qJ(3);
t68 = qJD(1) * qJD(3) + qJDD(1) * t86;
t183 = pkin(6) * t102;
t146 = t104 * t183;
t123 = -t105 * pkin(4) - t146;
t156 = qJD(4) * t102;
t106 = cos(pkin(7));
t88 = -t106 * pkin(1) - pkin(2);
t67 = -t105 * pkin(3) - t102 * qJ(4) + t88;
t29 = -qJD(1) * t156 + t67 * qJDD(1) + qJDD(3);
t42 = t102 * qJDD(2) + t105 * t68;
t14 = -t101 * t42 + t104 * t29;
t10 = t123 * qJDD(1) + t14;
t148 = t102 * qJDD(1);
t140 = t101 * t148;
t15 = t101 * t29 + t104 * t42;
t11 = -pkin(6) * t140 + t15;
t43 = t67 * qJD(1) + qJD(3);
t78 = t86 * qJD(1);
t56 = t102 * qJD(2) + t105 * t78;
t17 = -t101 * t56 + t104 * t43;
t13 = t123 * qJD(1) + t17;
t18 = t101 * t43 + t104 * t56;
t16 = -pkin(6) * t144 + t18;
t125 = t108 * t16 - t110 * t13;
t1 = -t125 * qJD(5) + t108 * t10 + t110 * t11;
t49 = qJD(1) * t145 - t135;
t186 = t49 ^ 2;
t100 = qJ(1) + pkin(7);
t93 = cos(t100);
t185 = g(3) * t93;
t184 = pkin(4) * t101;
t182 = g(1) * t102;
t111 = cos(qJ(1));
t181 = t111 * pkin(1);
t118 = qJD(1) * t71;
t46 = t102 * t118;
t180 = t46 * t84;
t179 = t49 * t46;
t178 = t49 * t84;
t52 = t102 * t188;
t161 = t108 * t101;
t70 = t160 - t161;
t58 = t70 * t102;
t177 = -t58 * t20 + t52 * t46;
t53 = (-qJD(5) * t161 + t143) * t102;
t57 = t71 * t102;
t147 = t105 * qJDD(1);
t83 = -qJDD(5) + t147;
t176 = t53 * t84 + t57 * t83;
t175 = -t105 * t118 + t188;
t174 = t84 * t70;
t167 = t105 * t86;
t28 = t101 * t67 + t104 * t167;
t169 = t102 * t93;
t91 = sin(t100);
t170 = t102 * t91;
t173 = g(2) * t169 + g(3) * t170;
t95 = t101 ^ 2;
t97 = t104 ^ 2;
t172 = -t95 - t97;
t96 = t102 ^ 2;
t98 = t105 ^ 2;
t171 = t96 + t98;
t168 = t105 * t20;
t166 = t105 * t91;
t165 = t105 * t93;
t164 = t19 * t105;
t112 = qJD(1) ^ 2;
t162 = t105 * t112;
t158 = qJD(3) * t102;
t157 = qJD(3) * t105;
t154 = qJDD(1) * t88;
t152 = t96 * qJDD(1);
t149 = qJDD(1) * t104;
t109 = sin(qJ(1));
t142 = -t109 * pkin(1) + t93 * qJ(3);
t141 = t171 * t112;
t138 = t104 * t147;
t55 = t105 * qJD(2) - t102 * t78;
t74 = qJDD(3) + t154;
t137 = t74 + t154;
t41 = t105 * qJDD(2) - t102 * t68;
t136 = 0.2e1 * t102 * t147;
t133 = g(2) * t93 + g(3) * t91;
t132 = g(2) * t91 - t185;
t131 = g(2) * t111 + g(3) * t109;
t130 = -t57 * t19 + t49 * t53;
t129 = -t52 * t84 + t58 * t83;
t54 = qJD(4) - t55;
t127 = -t41 * t102 + t42 * t105;
t126 = t55 * t102 - t56 * t105;
t6 = t108 * t13 + t110 * t16;
t60 = t104 * t67;
t21 = -t146 + t60 + (-t101 * t86 - pkin(4)) * t105;
t23 = -t101 * t183 + t28;
t7 = -t108 * t23 + t110 * t21;
t8 = t108 * t21 + t110 * t23;
t36 = qJDD(4) - t41;
t124 = -t91 * pkin(2) + t142;
t122 = t102 * (-pkin(6) - qJ(4)) - t105 * (t104 * pkin(4) + pkin(3)) - pkin(2);
t27 = -t101 * t167 + t60;
t63 = -t101 * t157 - t104 * t156;
t121 = -t63 * qJD(1) - t27 * qJDD(1) - t14;
t64 = -t101 * t156 + t104 * t157;
t120 = t64 * qJD(1) + t28 * qJDD(1) + t15;
t119 = g(1) * t105 + g(2) * t170 + t36;
t116 = -t93 * pkin(2) - t91 * qJ(3) - t181;
t2 = -t6 * qJD(5) + t110 * t10 - t108 * t11;
t115 = t36 * t102 + t68 * t96 + t132;
t99 = pkin(9) + qJ(5);
t92 = cos(t99);
t90 = sin(t99);
t89 = t98 * qJDD(1);
t82 = pkin(4) * t140;
t61 = (t86 + t184) * t102;
t45 = t46 ^ 2;
t40 = -t92 * t165 - t91 * t90;
t39 = t90 * t165 - t91 * t92;
t38 = t92 * t166 - t93 * t90;
t37 = t90 * t166 + t93 * t92;
t31 = pkin(4) * t144 + t54;
t26 = t36 + t82;
t4 = -t8 * qJD(5) - t108 * t64 + t110 * t63;
t3 = t7 * qJD(5) + t108 * t63 + t110 * t64;
t5 = [0, 0, 0, 0, 0, qJDD(1), t131, -g(2) * t109 + g(3) * t111, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t106 * t163 + t133, -0.2e1 * t103 * t163 - t132, 0, (t131 + (t103 ^ 2 + t106 ^ 2) * t163) * pkin(1), t152, t136, 0, t89, 0, 0, (t133 - t137) * t105, t137 * t102 - t173, t171 * t68 + t127 + t132, -g(2) * t116 - g(3) * t124 - t126 * qJD(3) + t127 * t86 + t74 * t88, t97 * t152, -0.2e1 * t96 * t101 * t149, -0.2e1 * t102 * t138, t95 * t152, t101 * t136, t89, (t133 * t104 + t121) * t105 + t115 * t101, (-t133 * t101 + t120) * t105 + t115 * t104, (-t120 * t101 + t121 * t104) * t102 + t173, t15 * t28 + t18 * t64 + t14 * t27 + t17 * t63 - g(2) * (-pkin(3) * t165 + t116) - g(3) * (-pkin(3) * t166 + t124) + (t133 * qJ(4) + t54 * qJD(3) + t36 * t86) * t102, -t19 * t58 - t49 * t52, -t130 + t177, -t129 + t164, t20 * t57 + t46 * t53, t168 + t176, t83 * t105, -g(2) * t40 + g(3) * t38 - t2 * t105 + t46 * t158 + t61 * t20 + t26 * t57 + t31 * t53 - t4 * t84 - t7 * t83, -g(2) * t39 - g(3) * t37 + t1 * t105 + t49 * t158 - t61 * t19 + t26 * t58 + t3 * t84 - t31 * t52 + t8 * t83, -t1 * t57 - t125 * t52 + t7 * t19 - t2 * t58 - t8 * t20 - t3 * t46 - t4 * t49 - t6 * t53 + t173, t1 * t8 + t6 * t3 + t2 * t7 - t125 * t4 + t26 * t61 + t31 * t158 + g(2) * t181 - g(3) * t142 + (-g(2) * t122 - g(3) * t184) * t93 + (-g(2) * (-qJ(3) - t184) - g(3) * t122) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t102 + t41 * t105 - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36 * t105 - g(1) + (-t101 * t14 + t104 * t15) * t102, 0, 0, 0, 0, 0, 0, -t168 + t176, t129 + t164, t130 + t177, t1 * t58 - t26 * t105 + t125 * t53 - t2 * t57 - t6 * t52 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t148, -t141, t126 * qJD(1) - t133 + t74, 0, 0, 0, 0, 0, 0, -t101 * t141 - t138, t101 * t147 - t104 * t141, t172 * t148, t15 * t101 + t14 * t104 + (-t102 * t54 + (t101 * t17 - t104 * t18) * t105) * qJD(1) - t133, 0, 0, 0, 0, 0, 0, -t46 * t159 + t175 * t84 - t70 * t83, -t49 * t159 - t174 * t84 + t71 * t83, t174 * t46 + t175 * t49 + t70 * t19 - t71 * t20, t1 * t71 + t125 * t175 - t31 * t159 - t174 * t6 + t2 * t70 - t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t104 * t162 + t150) * t102, (t101 * t162 + t149) * t102, t172 * t96 * t112, (-t185 + (t101 * t18 + t104 * t17) * qJD(1)) * t102 + t119, 0, 0, 0, 0, 0, 0, t20 - t178, -t19 + t180, -t45 - t186, -g(3) * t169 - t125 * t49 + t46 * t6 + t119 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t45 + t186, -t19 - t180, -t179, -t178 - t20, -t83, -g(2) * t37 + g(3) * t39 + t90 * t182 - t31 * t49 - t6 * t84 + t2, -g(2) * t38 - g(3) * t40 + t125 * t84 + t92 * t182 + t31 * t46 - t1, 0, 0;];
tau_reg = t5;
