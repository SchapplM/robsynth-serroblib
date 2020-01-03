% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP6
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:16
% EndTime: 2019-12-31 18:43:21
% DurationCPUTime: 1.51s
% Computational Cost: add. (3812->230), mult. (7440->294), div. (0->0), fcn. (4651->8), ass. (0->156)
t145 = sin(qJ(4));
t148 = cos(qJ(4));
t146 = sin(qJ(3));
t178 = qJD(1) * t146;
t109 = -t148 * qJD(3) + t145 * t178;
t149 = cos(qJ(3));
t174 = qJD(1) * qJD(3);
t166 = t149 * t174;
t173 = t146 * qJDD(1);
t116 = t166 + t173;
t84 = -qJD(4) * t109 + qJDD(3) * t145 + t116 * t148;
t127 = qJD(1) * t149 - qJD(4);
t97 = t109 * t127;
t68 = t84 - t97;
t213 = qJ(5) * t68;
t131 = t146 * t174;
t172 = t149 * qJDD(1);
t117 = -t131 + t172;
t108 = -qJDD(4) + t117;
t111 = qJD(3) * t145 + t148 * t178;
t89 = t111 * t109;
t206 = -t108 - t89;
t212 = pkin(4) * t206;
t106 = t109 ^ 2;
t179 = -g(3) + qJDD(2);
t130 = t149 * t179;
t151 = qJD(1) ^ 2;
t196 = pkin(3) * t149;
t160 = -pkin(7) * t146 - t196;
t147 = sin(qJ(1));
t150 = cos(qJ(1));
t165 = t147 * g(1) - g(2) * t150;
t112 = qJDD(1) * pkin(1) + t165;
t159 = g(1) * t150 + g(2) * t147;
t113 = -pkin(1) * t151 - t159;
t141 = sin(pkin(8));
t142 = cos(pkin(8));
t187 = t141 * t112 + t142 * t113;
t79 = -pkin(2) * t151 + qJDD(1) * pkin(6) + t187;
t162 = t151 * t160 + t79;
t204 = qJD(3) ^ 2;
t54 = -qJDD(3) * pkin(3) - t204 * pkin(7) + t162 * t146 - t130;
t161 = -t148 * qJDD(3) + t116 * t145;
t83 = -qJD(4) * t111 - t161;
t92 = -pkin(4) * t127 - qJ(5) * t111;
t20 = -t83 * pkin(4) - t106 * qJ(5) + t111 * t92 + qJDD(5) + t54;
t211 = t145 * t206;
t210 = t148 * t206;
t157 = -t117 + t131;
t158 = t116 + t166;
t164 = t142 * t112 - t141 * t113;
t78 = -qJDD(1) * pkin(2) - t151 * pkin(6) - t164;
t48 = t157 * pkin(3) - t158 * pkin(7) + t78;
t163 = t146 * t179;
t55 = -t204 * pkin(3) + qJDD(3) * pkin(7) + t162 * t149 + t163;
t23 = t145 * t48 + t148 * t55;
t154 = t83 * qJ(5) - 0.2e1 * qJD(5) * t109 + t127 * t92 + t23;
t107 = t111 ^ 2;
t125 = t127 ^ 2;
t86 = -t107 - t125;
t209 = -t154 + (t106 + t86) * pkin(4);
t207 = t84 + t97;
t64 = (qJD(4) + t127) * t111 + t161;
t85 = -t125 - t106;
t46 = t145 * t85 + t210;
t203 = pkin(3) * t46;
t75 = t108 - t89;
t190 = t145 * t75;
t51 = t148 * t86 + t190;
t202 = pkin(3) * t51;
t176 = qJD(5) * t111;
t102 = -0.2e1 * t176;
t22 = t145 * t55 - t148 * t48;
t153 = t212 - t22 - t213;
t11 = t102 + t153;
t201 = pkin(4) * t11;
t200 = pkin(4) * t68;
t36 = -t145 * t64 - t148 * t68;
t199 = pkin(7) * t36;
t198 = pkin(7) * t46;
t197 = pkin(7) * t51;
t37 = t145 * t68 - t148 * t64;
t74 = -t106 - t107;
t195 = -pkin(3) * t74 + pkin(7) * t37;
t47 = t148 * t85 - t211;
t63 = (qJD(4) - t127) * t111 + t161;
t194 = -pkin(3) * t63 + pkin(7) * t47;
t188 = t148 * t75;
t52 = -t145 * t86 + t188;
t193 = -pkin(3) * t207 + pkin(7) * t52;
t191 = t145 * t54;
t189 = t148 * t54;
t185 = qJ(5) * t145;
t184 = qJ(5) * t148;
t183 = t127 * t145;
t182 = t127 * t148;
t126 = t149 * t151 * t146;
t121 = qJDD(3) + t126;
t181 = t146 * t121;
t122 = qJDD(3) - t126;
t180 = t149 * t122;
t19 = t146 * t74 + t149 * t37;
t171 = pkin(1) * (t141 * t19 - t142 * t36) + pkin(6) * t19 - pkin(2) * t36;
t27 = t146 * t63 + t149 * t47;
t170 = pkin(1) * (t141 * t27 - t142 * t46) + pkin(6) * t27 - pkin(2) * t46;
t30 = t146 * t207 + t149 * t52;
t169 = pkin(1) * (t141 * t30 - t142 * t51) + pkin(6) * t30 - pkin(2) * t51;
t168 = t149 * t89;
t9 = t145 * t22 + t148 * t23;
t70 = t146 * t79 - t130;
t71 = t149 * t79 + t163;
t38 = t146 * t70 + t149 * t71;
t156 = t145 * t23 - t148 * t22;
t152 = t153 + t212;
t138 = t149 ^ 2;
t137 = t146 ^ 2;
t135 = t138 * t151;
t133 = t137 * t151;
t124 = -t135 - t204;
t123 = -t133 - t204;
t120 = t133 + t135;
t119 = (t137 + t138) * qJDD(1);
t118 = -0.2e1 * t131 + t172;
t115 = 0.2e1 * t166 + t173;
t103 = 0.2e1 * t176;
t94 = -t107 + t125;
t93 = t106 - t125;
t91 = -t123 * t146 - t180;
t90 = t124 * t149 - t181;
t87 = t107 - t106;
t72 = (t109 * t145 + t111 * t148) * t127;
t60 = -t111 * t182 + t145 * t84;
t59 = -t109 * t183 + t148 * t83;
t58 = t149 * t108 + t146 * (t109 * t148 - t111 * t145) * t127;
t57 = t145 * t93 - t188;
t56 = t148 * t94 + t211;
t41 = t146 * (t111 * t183 + t148 * t84) - t168;
t40 = t146 * (-t109 * t182 - t145 * t83) + t168;
t39 = -pkin(4) * t207 + qJ(5) * t75;
t35 = -t145 * t63 + t148 * t207;
t32 = t146 * (t148 * t93 + t190) + t149 * t64;
t31 = t146 * (-t145 * t94 + t210) - t149 * t68;
t29 = t146 * t52 - t149 * t207;
t26 = t146 * t47 - t149 * t63;
t24 = t146 * (-t145 * t207 - t148 * t63) - t149 * t87;
t18 = t146 * t37 - t149 * t74;
t16 = -qJ(5) * t86 + t20;
t13 = -pkin(4) * t106 + t154;
t12 = -pkin(4) * t63 + qJ(5) * t85 - t20;
t7 = t103 - t153 + t213;
t6 = -qJ(5) * t64 + (-t106 - t74) * pkin(4) + t154;
t4 = -pkin(4) * t20 + qJ(5) * t13;
t3 = -t11 * t145 + t13 * t148;
t2 = t11 * t148 + t13 * t145;
t1 = t146 * t20 + t149 * t3;
t5 = [0, 0, 0, 0, 0, qJDD(1), t165, t159, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t142 - t141 * t151) + t164, pkin(1) * (-qJDD(1) * t141 - t142 * t151) - t187, 0, pkin(1) * (t141 * t187 + t142 * t164), t158 * t146, t115 * t149 + t118 * t146, t181 + t149 * (-t133 + t204), -t157 * t149, t146 * (t135 - t204) + t180, 0, -t149 * t78 + pkin(2) * t118 + pkin(6) * t90 + pkin(1) * (t118 * t142 + t141 * t90), t146 * t78 - pkin(2) * t115 + pkin(6) * t91 + pkin(1) * (-t115 * t142 + t141 * t91), pkin(2) * t120 + pkin(6) * t119 + pkin(1) * (t119 * t141 + t120 * t142) + t38, -pkin(2) * t78 + pkin(6) * t38 + pkin(1) * (t141 * t38 - t142 * t78), t41, t24, t31, t40, t32, t58, t146 * (t191 - t198) + t149 * (t22 - t203) + t170, t146 * (t189 - t197) + t149 * (t23 - t202) + t169, t146 * (-t156 - t199) - t36 * t196 + t171, (pkin(1) * t141 + pkin(6)) * (t146 * t54 + t149 * t9) + (-pkin(1) * t142 - pkin(2) + t160) * t156, t41, t24, t31, t40, t32, t58, t146 * (-t12 * t145 - t184 * t206 - t198) + t149 * (t103 - t152 - t203) + t170, t146 * (-t145 * t39 + t148 * t16 - t197) + t149 * (-t202 - t209) + t169, t146 * (-t145 * t6 + t148 * t7 - t199) + t149 * (-pkin(3) * t36 + t200) + t171, t146 * (-pkin(7) * t2 - t11 * t184 - t145 * t4) + t149 * (-pkin(3) * t2 - t201) - pkin(2) * t2 + pkin(6) * t1 + pkin(1) * (t1 * t141 - t142 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, 0, 0, 0, 0, 0, 0, t121 * t149 + t124 * t146, -t122 * t146 + t123 * t149, 0, t146 * t71 - t149 * t70, 0, 0, 0, 0, 0, 0, t26, t29, t18, t146 * t9 - t149 * t54, 0, 0, 0, 0, 0, 0, t26, t29, t18, t146 * t3 - t149 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, t133 - t135, t173, t126, t172, qJDD(3), -t70, -t71, 0, 0, t60, t35, t56, t59, t57, t72, -t189 + t194, t191 + t193, t9 + t195, -pkin(3) * t54 + pkin(7) * t9, t60, t35, t56, t59, t57, t72, t12 * t148 - t185 * t206 + t194, t145 * t16 + t148 * t39 + t193, t145 * t7 + t148 * t6 + t195, -pkin(3) * t20 + pkin(7) * t3 - t11 * t185 + t148 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t87, t68, -t89, -t64, -t108, -t22, -t23, 0, 0, t89, t87, t68, -t89, -t64, -t108, t102 + t152, t209, -t200, t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t207, t74, t20;];
tauJ_reg = t5;
