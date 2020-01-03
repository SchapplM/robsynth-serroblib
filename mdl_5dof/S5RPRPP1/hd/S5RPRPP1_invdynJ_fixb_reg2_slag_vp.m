% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:19
% EndTime: 2019-12-31 18:09:21
% DurationCPUTime: 1.38s
% Computational Cost: add. (1748->272), mult. (3687->319), div. (0->0), fcn. (2431->12), ass. (0->150)
t106 = sin(pkin(8));
t108 = cos(pkin(8));
t113 = cos(qJ(3));
t154 = t113 * qJDD(1);
t111 = sin(qJ(3));
t155 = t111 * qJDD(1);
t133 = t106 * t155 - t108 * t154;
t66 = t106 * t113 + t108 * t111;
t60 = t66 * qJD(3);
t35 = qJD(1) * t60 + t133;
t167 = t108 * t113;
t149 = qJD(1) * t167;
t163 = qJD(1) * t111;
t58 = t106 * t163 - t149;
t161 = qJD(3) * t113;
t162 = qJD(3) * t111;
t63 = -t106 * t162 + t108 * t161;
t176 = -t35 * t66 - t58 * t63;
t190 = t66 * qJD(1);
t158 = qJD(1) * qJD(3);
t148 = t111 * t158;
t125 = qJDD(1) * t66 - t106 * t148;
t147 = t113 * t158;
t36 = t108 * t147 + t125;
t65 = t106 * t111 - t167;
t177 = t190 * t60 + t36 * t65;
t196 = t176 + t177;
t195 = t177 - t176;
t194 = 0.2e1 * qJD(3);
t107 = sin(pkin(7));
t87 = pkin(1) * t107 + pkin(6);
t75 = t87 * qJDD(1);
t143 = -qJD(2) * qJD(3) - t75;
t102 = qJ(1) + pkin(7);
t93 = sin(t102);
t95 = cos(t102);
t140 = g(1) * t95 + g(2) * t93;
t109 = cos(pkin(7));
t89 = -pkin(1) * t109 - pkin(2);
t99 = t113 * pkin(3);
t193 = t89 - t99;
t187 = t190 ^ 2;
t54 = t58 ^ 2;
t192 = -t54 - t187;
t191 = -t54 + t187;
t186 = g(1) * t93;
t150 = g(2) * t95 - t186;
t169 = pkin(1) * qJDD(1);
t171 = qJ(4) + t87;
t146 = t171 * t111;
t64 = t171 * t113;
t34 = -t106 * t146 + t108 * t64;
t101 = qJ(3) + pkin(8);
t92 = sin(t101);
t189 = -t34 * qJDD(3) + t150 * t92;
t77 = t87 * qJD(1);
t144 = qJ(4) * qJD(1) + t77;
t157 = qJD(1) * qJD(4);
t97 = t113 * qJDD(2);
t15 = qJDD(3) * pkin(3) + t97 - t144 * t161 + (-qJ(4) * qJDD(1) + t143 - t157) * t111;
t153 = -t111 * qJDD(2) + t113 * t143;
t31 = -t162 * t77 - t153;
t18 = t113 * t157 + (-t148 + t154) * qJ(4) + t31;
t4 = t106 * t15 + t108 * t18;
t94 = cos(t101);
t188 = g(3) * t92 + t140 * t94 - t4;
t183 = pkin(3) * t111;
t182 = g(3) * t113;
t112 = sin(qJ(1));
t181 = t112 * pkin(1);
t180 = t190 * t58;
t178 = t94 * t95;
t3 = -t106 * t18 + t108 * t15;
t160 = t111 * qJD(2);
t49 = t113 * t144 + t160;
t43 = t108 * t49;
t98 = t113 * qJD(2);
t48 = -t111 * t144 + t98;
t45 = qJD(3) * pkin(3) + t48;
t17 = t106 * t45 + t43;
t175 = t106 * t49;
t174 = t111 * t77;
t173 = t113 * t77;
t172 = t92 * qJ(5);
t114 = cos(qJ(1));
t100 = t114 * pkin(1);
t91 = t99 + pkin(2);
t170 = t91 * t95 + t100;
t78 = qJD(1) * t89;
t168 = qJDD(3) * pkin(4);
t20 = t106 * t48 + t43;
t166 = t20 * qJD(3);
t21 = t108 * t48 - t175;
t165 = qJD(5) - t21;
t104 = t111 ^ 2;
t105 = t113 ^ 2;
t164 = t104 - t105;
t76 = qJDD(1) * t89;
t152 = pkin(3) * t162;
t116 = qJD(1) ^ 2;
t151 = t111 * t116 * t113;
t145 = qJD(3) * t171;
t141 = t111 * t147;
t138 = -pkin(4) * t94 - t172;
t137 = g(1) * t112 - g(2) * t114;
t136 = t35 * t65 + t58 * t60;
t110 = -qJ(4) - pkin(6);
t134 = -t95 * t110 - t181;
t16 = t108 * t45 - t175;
t52 = t160 + t173;
t39 = qJD(3) * t60 + qJDD(3) * t65;
t131 = -g(3) * t94 + t140 * t92 + t3;
t33 = t106 * t64 + t108 * t146;
t130 = -g(2) * t178 - t33 * qJDD(3) + t186 * t94;
t127 = -qJD(1) * t78 + t140;
t126 = -t111 * qJD(4) - t113 * t145;
t124 = -qJDD(3) * t87 + t194 * t78;
t57 = qJD(1) * t193 + qJD(4);
t26 = t58 * pkin(4) - qJ(5) * t190 + t57;
t123 = -t190 * t26 - qJDD(5) + t131;
t47 = pkin(3) * t148 + qJDD(1) * t193 + qJDD(4);
t115 = qJD(3) ^ 2;
t122 = -t115 * t87 - t150 - 0.2e1 * t76;
t50 = t113 * qJD(4) - t111 * t145;
t24 = t106 * t50 - t108 * t126;
t25 = t106 * t126 + t108 * t50;
t121 = t190 * t24 - t25 * t58 + t33 * t36 - t34 * t35 - t140;
t32 = -t52 * qJD(3) - t111 * t75 + t97;
t51 = t98 - t174;
t120 = -t32 * t111 + t31 * t113 + (-t111 * t52 - t113 * t51) * qJD(3);
t119 = t35 * pkin(4) - t36 * qJ(5) + t47;
t118 = t190 * t194 + t133;
t103 = qJDD(3) * qJ(5);
t88 = -pkin(3) * t108 - pkin(4);
t84 = pkin(3) * t106 + qJ(5);
t74 = qJDD(3) * t113 - t111 * t115;
t73 = qJDD(3) * t111 + t113 * t115;
t38 = qJD(3) * t63 + qJDD(3) * t66;
t30 = pkin(4) * t65 - qJ(5) * t66 + t193;
t29 = pkin(3) * t163 + pkin(4) * t190 + qJ(5) * t58;
t23 = (t58 + t149) * qJD(3) + t125;
t22 = (-t58 + t149) * qJD(3) + t125;
t19 = pkin(4) * t60 - qJ(5) * t63 - qJD(5) * t66 + t152;
t12 = qJD(3) * qJ(5) + t17;
t11 = -qJD(3) * pkin(4) + qJD(5) - t16;
t6 = t190 * t63 + t36 * t66;
t5 = -qJD(5) * t190 + t119;
t2 = qJDD(5) - t168 - t3;
t1 = qJD(3) * qJD(5) + t103 + t4;
t7 = [0, 0, 0, 0, 0, qJDD(1), t137, g(1) * t114 + g(2) * t112, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t109 * t169 - t150, -0.2e1 * t107 * t169 + t140, 0, (t137 + (t107 ^ 2 + t109 ^ 2) * t169) * pkin(1), qJDD(1) * t104 + 0.2e1 * t141, 0.2e1 * t111 * t154 - 0.2e1 * t158 * t164, t73, qJDD(1) * t105 - 0.2e1 * t141, t74, 0, t111 * t124 + t113 * t122, -t111 * t122 + t113 * t124, (t104 + t105) * t75 + t120 - t140, t76 * t89 - g(1) * (-pkin(2) * t93 + pkin(6) * t95 - t181) - g(2) * (pkin(2) * t95 + pkin(6) * t93 + t100) + t120 * t87, t6, -t195, t38, t136, -t39, 0, t193 * t35 + t47 * t65 + t57 * t60 + (t183 * t58 - t24) * qJD(3) + t130, t193 * t36 + t47 * t66 + t57 * t63 + (t183 * t190 - t25) * qJD(3) + t189, -t16 * t63 - t17 * t60 - t3 * t66 - t4 * t65 + t121, t4 * t34 + t17 * t25 - t3 * t33 - t16 * t24 + t47 * t193 + t57 * t152 - g(1) * (-t91 * t93 + t134) - g(2) * (-t93 * t110 + t170), t6, t38, t195, 0, t39, t136, -qJD(3) * t24 + t19 * t58 + t26 * t60 + t30 * t35 + t5 * t65 + t130, -t1 * t65 + t11 * t63 - t12 * t60 + t2 * t66 + t121, t25 * qJD(3) - t19 * t190 - t26 * t63 - t30 * t36 - t5 * t66 - t189, t1 * t34 + t12 * t25 + t5 * t30 + t26 * t19 + t2 * t33 + t11 * t24 - g(1) * t134 - g(2) * (pkin(4) * t178 + t172 * t95 + t170) + (-g(1) * (t138 - t91) + g(2) * t110) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t74, -t73, 0, t31 * t111 + t32 * t113 - g(3) + (-t111 * t51 + t113 * t52) * qJD(3), 0, 0, 0, 0, 0, 0, -t39, -t38, t196, -t16 * t60 + t17 * t63 - t3 * t65 + t4 * t66 - g(3), 0, 0, 0, 0, 0, 0, -t39, t196, t38, t1 * t66 + t11 * t60 + t12 * t63 + t2 * t65 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, t164 * t116, t155, t151, t154, qJDD(3), -t182 + t97 + (t52 - t173) * qJD(3) + (t127 + t143) * t111, g(3) * t111 + (t51 + t174) * qJD(3) + t127 * t113 + t153, 0, 0, t180, t191, t23, -t180, -t133, qJDD(3), t166 - t57 * t190 + (qJDD(3) * t108 - t163 * t58) * pkin(3) + t131, t21 * qJD(3) + t57 * t58 + (-qJDD(3) * t106 - t163 * t190) * pkin(3) + t188, (t17 - t20) * t190 + (-t16 + t21) * t58 + (-t106 * t35 - t108 * t36) * pkin(3), t16 * t20 - t17 * t21 + (-t182 + t106 * t4 + t108 * t3 + (-qJD(1) * t57 + t140) * t111) * pkin(3), t180, t23, -t191, qJDD(3), t133, -t180, t166 - t29 * t58 + (pkin(4) - t88) * qJDD(3) + t123, -t84 * t35 + t88 * t36 + (t12 - t20) * t190 + (t11 - t165) * t58, t84 * qJDD(3) - t26 * t58 + t29 * t190 + t103 + (0.2e1 * qJD(5) - t21) * qJD(3) - t188, t1 * t84 + t2 * t88 - t26 * t29 - t11 * t20 - g(3) * (-t138 + t99) + t165 * t12 + t140 * (pkin(4) * t92 - qJ(5) * t94 + t183); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t22, t192, t16 * t190 + t17 * t58 + t150 + t47, 0, 0, 0, 0, 0, 0, t118, t192, -t22, t12 * t58 + (-qJD(5) - t11) * t190 + t119 + t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) + t180, t23, -t187 - t115, -qJD(3) * t12 - t123 - t168;];
tau_reg = t7;
