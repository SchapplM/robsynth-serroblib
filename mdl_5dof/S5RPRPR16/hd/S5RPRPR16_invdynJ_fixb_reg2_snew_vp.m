% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR16
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR16_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR16_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:32
% EndTime: 2019-12-31 18:39:36
% DurationCPUTime: 1.07s
% Computational Cost: add. (2301->211), mult. (4680->233), div. (0->0), fcn. (2440->6), ass. (0->143)
t162 = -pkin(6) - pkin(1);
t101 = sin(qJ(3));
t141 = qJD(1) * qJD(3);
t89 = t101 * t141;
t104 = cos(qJ(3));
t90 = t104 * qJDD(1);
t68 = t90 - t89;
t174 = t68 - t89;
t107 = qJD(1) ^ 2;
t145 = t104 * t107;
t82 = t101 * t145;
t75 = qJDD(3) - t82;
t148 = t104 * t75;
t106 = qJD(3) ^ 2;
t98 = t101 ^ 2;
t91 = t98 * t107;
t77 = -t91 - t106;
t43 = t101 * t77 + t148;
t173 = t162 * t43;
t74 = qJDD(3) + t82;
t154 = t101 * t74;
t99 = t104 ^ 2;
t92 = t99 * t107;
t79 = -t92 - t106;
t45 = -t104 * t79 + t154;
t172 = t162 * t45;
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t143 = qJD(1) * t101;
t60 = t100 * qJD(3) - t103 * t143;
t62 = t103 * qJD(3) + t100 * t143;
t41 = t62 * t60;
t59 = qJDD(5) + t68;
t168 = -t41 + t59;
t171 = t100 * t168;
t170 = t103 * t168;
t102 = sin(qJ(1));
t105 = cos(qJ(1));
t129 = t105 * g(1) + t102 * g(2);
t97 = qJDD(1) * qJ(2);
t122 = t129 - t97;
t133 = t104 * t141;
t113 = -pkin(3) * t133 + t174 * qJ(4) + t122;
t136 = t162 * t107;
t138 = t101 * qJDD(1);
t67 = t133 + t138;
t169 = t67 * pkin(3) - t113 + t136;
t142 = t104 * qJD(1);
t140 = qJD(2) * qJD(1);
t94 = 0.2e1 * t140;
t167 = -0.2e1 * qJD(4) * t142 + t94;
t146 = t104 * qJ(4);
t126 = t101 * pkin(3) - t146;
t63 = t126 * qJD(1);
t166 = qJDD(3) * qJ(4) - t63 * t143;
t165 = (pkin(4) * t141 + pkin(7) * t145 - g(3)) * t101 + t68 * pkin(4);
t164 = t63 * t142 + qJDD(4);
t57 = t60 ^ 2;
t58 = t62 ^ 2;
t85 = qJD(5) + t142;
t80 = t85 ^ 2;
t163 = pkin(3) + pkin(7);
t160 = t101 * g(3);
t132 = t100 * qJDD(3) - t103 * t67;
t115 = (-qJD(5) + t85) * t62 - t132;
t123 = t103 * qJDD(3) + t100 * t67;
t34 = -t60 * qJD(5) + t123;
t53 = t85 * t60;
t26 = t34 + t53;
t9 = t100 * t115 - t103 * t26;
t159 = t104 * t9;
t158 = t106 * pkin(3);
t134 = t102 * g(1) - t105 * g(2);
t130 = qJDD(2) - t134;
t119 = -t107 * qJ(2) + t130;
t52 = t162 * qJDD(1) + t119;
t48 = t101 * t52;
t38 = t104 * g(3) - t48;
t116 = -t38 + t166;
t111 = -t116 + t158;
t121 = pkin(4) * t142 - qJD(3) * pkin(7);
t139 = qJD(4) * qJD(3);
t14 = t67 * pkin(4) + pkin(7) * t91 - qJD(3) * t121 + t111 - 0.2e1 * t139;
t157 = t100 * t14;
t31 = t41 + t59;
t156 = t100 * t31;
t155 = t100 * t85;
t12 = -t121 * t142 + t163 * t67 + (-t98 * pkin(4) + t162) * t107 - t113 + t167;
t153 = t103 * t12;
t152 = t103 * t14;
t151 = t103 * t31;
t150 = t103 * t85;
t149 = t104 * t52;
t147 = qJDD(1) * pkin(1);
t144 = qJD(5) + t85;
t137 = t104 * t41;
t117 = -t106 * qJ(4) - t149 + t164;
t4 = t100 * t12 - t103 * (-t163 * qJDD(3) + t117 + t165);
t71 = (t98 + t99) * qJDD(1);
t73 = -t92 - t91;
t131 = qJ(2) * t73 - t162 * t71;
t112 = qJDD(3) * pkin(3) - t117;
t109 = -qJDD(3) * pkin(7) - t112;
t5 = t153 + (t109 + t165) * t100;
t2 = t100 * t5 - t103 * t4;
t128 = t100 * t4 + t103 * t5;
t37 = t149 + t160;
t127 = t101 * t14 + t104 * t2;
t20 = -t101 * t38 + t104 * t37;
t125 = -t104 * (t91 - t106) + t154;
t124 = t101 * (-t92 + t106) - t148;
t120 = qJ(2) + t126;
t118 = t101 * t163 + qJ(2) - t146;
t108 = t167 + t169;
t93 = 0.2e1 * t139;
t72 = t92 - t91;
t69 = t90 - 0.2e1 * t89;
t66 = 0.2e1 * t133 + t138;
t54 = -t119 + t147;
t51 = -t58 + t80;
t50 = t57 - t80;
t49 = t122 - t136 - 0.2e1 * t140;
t47 = t174 * t104;
t42 = (t67 + t133) * t101;
t40 = t58 - t57;
t39 = -t58 - t80;
t36 = -t101 * t69 - t104 * t66;
t35 = -t80 - t57;
t33 = -t62 * qJD(5) - t132;
t29 = -t57 - t58;
t28 = t112 + t160;
t27 = -t111 + t93;
t25 = t34 - t53;
t24 = -t144 * t60 + t123;
t21 = t144 * t62 + t132;
t18 = -t100 * t39 - t151;
t17 = t103 * t39 - t156;
t15 = t100 * t35 + t170;
t11 = t101 * t27 + t104 * t28;
t8 = t101 * t24 - t104 * t17;
t7 = t101 * t21 - t104 * t15;
t6 = t101 * t29 - t159;
t1 = [0, 0, 0, 0, 0, qJDD(1), t134, t129, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t130 - 0.2e1 * t147, -t129 + t94 + 0.2e1 * t97, pkin(1) * t54 + qJ(2) * (-t107 * pkin(1) - t122 + t94), t47, t36, -t124, t42, -t125, 0, qJ(2) * t66 - t101 * t49 + t173, qJ(2) * t69 - t104 * t49 - t172, t131 - t20, -qJ(2) * t49 + t162 * t20, 0, t124, t125, t47, t36, t42, t104 * (-qJ(4) * t73 - t112) - t101 * (-pkin(3) * t73 - t158 + t166 + t48 + t93) + t131, -t101 * t108 - t120 * t66 - t173, t104 * (0.2e1 * (qJD(4) * t104 - qJD(2)) * qJD(1) - t169) - t120 * t69 + t172, t120 * t108 + t162 * t11, t137 - t101 * (-t100 * t34 - t62 * t150), t104 * t40 - t101 * (t100 * t21 - t103 * t25), t104 * t26 - t101 * (-t103 * t51 - t171), -t137 - t101 * (-t103 * t33 - t60 * t155), t104 * t115 - t101 * (-t100 * t50 - t151), t104 * t59 - t101 * (t100 * t60 + t103 * t62) * t85, t104 * (pkin(4) * t15 - t4) - t101 * (pkin(4) * t21 - t152) + t162 * t7 + t118 * (t103 * t35 - t171), t104 * (-t153 - t100 * (pkin(7) * t82 + t109 - t160) - qJ(4) * t18 + (-t100 * (t68 + t89) + t17) * pkin(4)) - t101 * (pkin(4) * t24 - t163 * t18 + t157) + qJ(2) * t18 + t162 * t8, pkin(4) * t159 - t101 * (pkin(4) * t29 - t128) + t162 * t6 + t118 * (t100 * t26 + t103 * t115), t118 * t128 + (pkin(4) - t162) * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t107, -t54, 0, 0, 0, 0, 0, 0, t43, -t45, -t71, t20, 0, 0, 0, 0, 0, 0, -t71, -t43, t45, t11, 0, 0, 0, 0, 0, 0, t7, t8, t6, -t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t72, t90, -t82, -t138, qJDD(3), t37, t38, 0, 0, qJDD(3), -t90, t138, t82, t72, -t82, (-pkin(3) * t104 - qJ(4) * t101) * qJDD(1), (-t106 - t77) * qJ(4) + (-qJDD(3) - t75) * pkin(3) - t37 + t164, qJ(4) * t74 + t93 + (-t106 - t79) * pkin(3) + t116, pkin(3) * t28 + qJ(4) * t27, t103 * t34 - t62 * t155, -t100 * t25 - t103 * t21, -t100 * t51 + t170, -t100 * t33 + t60 * t150, t103 * t50 - t156, (t100 * t62 - t103 * t60) * t85, qJ(4) * t21 - t163 * t15 - t157, qJ(4) * t24 - t163 * t17 - t152, qJ(4) * t29 - t163 * t9 - t2, -qJ(4) * t14 - t163 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t75, t79, -t28, 0, 0, 0, 0, 0, 0, t15, t17, t9, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40, t26, -t41, t115, t59, -t4, -t5, 0, 0;];
tauJ_reg = t1;
