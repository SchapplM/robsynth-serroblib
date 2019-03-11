% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% tau_reg [6x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:05
% EndTime: 2019-03-09 01:32:09
% DurationCPUTime: 1.42s
% Computational Cost: add. (1750->244), mult. (3324->313), div. (0->0), fcn. (2438->14), ass. (0->143)
t106 = sin(pkin(10));
t108 = cos(pkin(10));
t111 = sin(qJ(5));
t114 = cos(qJ(5));
t64 = t114 * t106 + t111 * t108;
t56 = t64 * qJD(1);
t194 = qJD(6) + t56;
t110 = sin(qJ(6));
t113 = cos(qJ(6));
t160 = t113 * qJD(5);
t165 = qJD(1) * t108;
t151 = t114 * t165;
t166 = qJD(1) * t106;
t152 = t111 * t166;
t58 = t151 - t152;
t37 = t110 * t58 - t160;
t195 = t194 * t37;
t107 = sin(pkin(9));
t77 = t107 * pkin(1) + qJ(3);
t168 = t77 * qJDD(1);
t101 = qJ(1) + pkin(9);
t90 = sin(t101);
t92 = cos(t101);
t153 = -g(1) * t90 + g(2) * t92;
t100 = pkin(10) + qJ(5);
t89 = sin(t100);
t91 = cos(t100);
t121 = g(3) * t89 + t153 * t91;
t109 = cos(pkin(9));
t81 = -t109 * pkin(1) - pkin(2);
t74 = -qJ(4) + t81;
t62 = t74 * qJD(1) + qJD(3);
t40 = -t106 * qJD(2) + t108 * t62;
t33 = -pkin(7) * t165 + t40;
t41 = t108 * qJD(2) + t106 * t62;
t34 = -pkin(7) * t166 + t41;
t13 = t111 * t33 + t114 * t34;
t157 = t108 * qJDD(1);
t159 = qJD(1) * qJD(4);
t49 = t74 * qJDD(1) + qJDD(3) - t159;
t35 = -t106 * qJDD(2) + t108 * t49;
t31 = -pkin(7) * t157 + t35;
t158 = t106 * qJDD(1);
t36 = t108 * qJDD(2) + t106 * t49;
t32 = -pkin(7) * t158 + t36;
t131 = t111 * t32 - t114 * t31;
t2 = -qJDD(5) * pkin(5) + t13 * qJD(5) + t131;
t193 = t121 - (t58 * pkin(5) + t194 * pkin(8)) * t194 - t2;
t144 = t113 * t194;
t122 = -qJD(5) * t152 + t64 * qJDD(1);
t28 = qJD(5) * t151 + t122;
t25 = qJDD(6) + t28;
t173 = t110 * t25;
t192 = -t144 * t194 - t173;
t190 = t56 * qJD(5);
t189 = t81 * qJDD(1);
t163 = qJD(5) * t114;
t164 = qJD(5) * t111;
t59 = -t106 * t163 - t108 * t164;
t63 = t111 * t106 - t114 * t108;
t135 = -t194 * t59 + t25 * t63;
t162 = qJD(6) * t110;
t154 = t63 * t162;
t187 = -t113 * t135 + t154 * t194;
t129 = t111 * t34 - t114 * t33;
t10 = -qJD(5) * pkin(5) + t129;
t182 = -pkin(7) + t74;
t54 = t182 * t106;
t55 = t182 * t108;
t19 = t111 * t54 - t114 * t55;
t14 = -t64 * qJD(4) - t19 * qJD(5);
t130 = t111 * t31 + t114 * t32;
t67 = t77 * qJD(1);
t65 = qJD(4) + t67;
t53 = pkin(4) * t166 + t65;
t16 = t56 * pkin(5) - t58 * pkin(8) + t53;
t148 = qJDD(5) * pkin(8) - t129 * qJD(5) + qJD(6) * t16 + t130;
t66 = t106 * pkin(4) + t77;
t18 = t64 * pkin(5) + t63 * pkin(8) + t66;
t20 = t111 * t55 + t114 * t54;
t186 = t10 * t59 - (qJD(6) * t18 + t14) * t194 - t148 * t64 - t2 * t63 - t20 * t25;
t183 = g(3) * t91;
t39 = t110 * qJD(5) + t113 * t58;
t60 = -t106 * t164 + t108 * t163;
t134 = -t111 * t158 + t114 * t157;
t27 = t134 - t190;
t8 = qJD(6) * t160 + t110 * qJDD(5) + t113 * t27 - t58 * t162;
t181 = t39 * t60 + t8 * t64;
t180 = t10 * t63;
t178 = t18 * t25;
t177 = t39 * t58;
t176 = t58 * t37;
t175 = t8 * t110;
t174 = t106 ^ 2 + t108 ^ 2;
t21 = t113 * t25;
t172 = t90 * t110;
t171 = t90 * t113;
t170 = t92 * t110;
t169 = t92 * t113;
t161 = qJD(6) * t113;
t103 = qJD(3) * qJD(1);
t115 = cos(qJ(1));
t156 = t115 * pkin(1) + t92 * pkin(2) + t90 * qJ(3);
t155 = -t103 - t168;
t112 = sin(qJ(1));
t150 = -t112 * pkin(1) + t92 * qJ(3);
t11 = qJD(5) * pkin(8) + t13;
t61 = qJDD(4) - t155;
t48 = pkin(4) * t158 + t61;
t6 = t28 * pkin(5) - t27 * pkin(8) + t48;
t147 = qJD(6) * t11 - t6;
t146 = -t113 * qJDD(5) + t110 * t27;
t143 = t174 * qJDD(1);
t142 = qJD(6) * t64 + qJD(1);
t140 = -g(1) * t92 - g(2) * t90;
t9 = t39 * qJD(6) + t146;
t138 = -t60 * t37 - t64 * t9;
t137 = t59 * t39 - t63 * t8;
t136 = g(1) * t112 - g(2) * t115;
t133 = t36 * t106 + t35 * t108;
t132 = t41 * t106 + t40 * t108;
t29 = t59 * qJD(5) - t63 * qJDD(5);
t30 = -t60 * qJD(5) - t64 * qJDD(5);
t128 = t21 + (-t110 * t56 - t162) * t194;
t127 = -t148 + t183;
t126 = qJDD(3) + t189;
t124 = t140 + t168;
t123 = -t133 - t153;
t120 = -pkin(8) * t25 + (t10 - t129) * t194;
t119 = t124 + t61 + t103;
t118 = t161 * t194 * t63 + t135 * t110;
t116 = qJD(1) ^ 2;
t105 = qJDD(2) - g(3);
t47 = t89 * t169 - t172;
t46 = t89 * t170 + t171;
t45 = t89 * t171 + t170;
t44 = -t89 * t172 + t169;
t23 = t60 * pkin(5) - t59 * pkin(8) + qJD(3);
t15 = -t63 * qJD(4) + t20 * qJD(5);
t5 = t113 * t6;
t4 = t113 * t11 + t110 * t16;
t3 = -t110 * t11 + t113 * t16;
t1 = [qJDD(1), t136, g(1) * t115 + g(2) * t112 (t136 + (t107 ^ 2 + t109 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(3) + t153 + 0.2e1 * t189, 0.2e1 * t103 + t124 + t168, -t155 * t77 + t67 * qJD(3) + t126 * t81 - g(1) * (-t90 * pkin(2) + t150) - g(2) * t156, t119 * t106, t119 * t108, -t74 * t143 + t174 * t159 + t123, t61 * t77 + t65 * qJD(3) - g(1) * ((-pkin(2) - qJ(4)) * t90 + t150) - g(2) * (t92 * qJ(4) + t156) + t133 * t74 - t132 * qJD(4), -t27 * t63 + t58 * t59, -t27 * t64 + t63 * t28 - t59 * t56 - t58 * t60, t29, t30, 0, qJD(3) * t56 - t15 * qJD(5) - t19 * qJDD(5) + t140 * t89 + t66 * t28 + t48 * t64 + t53 * t60, qJD(3) * t58 - t14 * qJD(5) - t20 * qJDD(5) + t140 * t91 + t66 * t27 - t48 * t63 + t53 * t59, t113 * t137 + t154 * t39 (-t110 * t39 - t113 * t37) * t59 + (t175 + t113 * t9 + (-t110 * t37 + t113 * t39) * qJD(6)) * t63, t181 + t187, t118 + t138, t194 * t60 + t25 * t64, -g(1) * t47 - g(2) * t45 + t15 * t37 + t19 * t9 + t3 * t60 + t5 * t64 + (t178 + t23 * t194 + (-t11 * t64 - t194 * t20 - t180) * qJD(6)) * t113 + t186 * t110, g(1) * t46 - g(2) * t44 + t15 * t39 + t19 * t8 - t4 * t60 + (-(-qJD(6) * t20 + t23) * t194 - t178 + t147 * t64 + qJD(6) * t180) * t110 + t186 * t113; 0, 0, 0, t105, 0, 0, t105, 0, 0, 0, -t35 * t106 + t36 * t108 - g(3), 0, 0, 0, 0, 0, t30, -t29, 0, 0, 0, 0, 0, t118 - t138, t181 - t187; 0, 0, 0, 0, qJDD(1), -t116, -t67 * qJD(1) + t126 + t153, -t116 * t106, -t116 * t108, -t143, -t65 * qJD(1) - t123, 0, 0, 0, 0, 0, -qJD(1) * t56 + t29, -qJD(1) * t58 + t30, 0, 0, 0, 0, 0, -t64 * t173 - t59 * t37 + t63 * t9 + (-t110 * t60 - t113 * t142) * t194, -t64 * t21 + (t110 * t142 - t113 * t60) * t194 - t137; 0, 0, 0, 0, 0, 0, 0, t158, t157, -t174 * t116, t132 * qJD(1) + t140 + t61, 0, 0, 0, 0, 0 (t58 + t151) * qJD(5) + t122, t134 - 0.2e1 * t190, 0, 0, 0, 0, 0, t128 - t176, -t177 + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t56, -t56 ^ 2 + t58 ^ 2, t134 (t58 - t151) * qJD(5) - t122, qJDD(5), -t53 * t58 + t121 - t131, -t153 * t89 + t53 * t56 - t130 + t183, t144 * t39 + t175 (t8 - t195) * t113 + (-t194 * t39 - t9) * t110, -t177 - t192, t128 + t176, -t194 * t58, -pkin(5) * t9 + t120 * t110 + t193 * t113 - t13 * t37 - t3 * t58, -pkin(5) * t8 - t193 * t110 + t120 * t113 - t13 * t39 + t4 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t37, -t37 ^ 2 + t39 ^ 2, t8 + t195, -t146 + (-qJD(6) + t194) * t39, t25, -g(1) * t44 - g(2) * t46 - t10 * t39 - t11 * t161 + t110 * t127 + t194 * t4 + t5, g(1) * t45 - g(2) * t47 + t10 * t37 + t110 * t147 + t113 * t127 + t194 * t3;];
tau_reg  = t1;
