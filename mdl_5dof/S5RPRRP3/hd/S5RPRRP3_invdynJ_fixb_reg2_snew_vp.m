% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP3
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:47:43
% EndTime: 2020-01-03 11:47:49
% DurationCPUTime: 1.63s
% Computational Cost: add. (4103->227), mult. (8434->299), div. (0->0), fcn. (5393->8), ass. (0->141)
t133 = qJDD(3) + qJDD(4);
t142 = sin(qJ(4));
t144 = cos(qJ(4));
t145 = cos(qJ(3));
t143 = sin(qJ(3));
t170 = qJD(1) * t143;
t106 = -t144 * t145 * qJD(1) + t142 * t170;
t108 = (t145 * t142 + t143 * t144) * qJD(1);
t84 = t108 * t106;
t193 = t84 - t133;
t194 = t193 * pkin(4);
t192 = t142 * t193;
t191 = t144 * t193;
t167 = qJD(1) * qJD(3);
t158 = t145 * t167;
t166 = t143 * qJDD(1);
t113 = t158 + t166;
t129 = t145 * qJDD(1);
t159 = t143 * t167;
t114 = t129 - t159;
t73 = -t106 * qJD(4) + t144 * t113 + t142 * t114;
t134 = qJD(3) + qJD(4);
t95 = t134 * t106;
t190 = t73 - t95;
t147 = qJD(1) ^ 2;
t124 = t145 * t147 * t143;
t118 = qJDD(3) + t124;
t137 = -g(1) + qJDD(2);
t183 = sin(qJ(1));
t184 = cos(qJ(1));
t151 = g(2) * t183 - g(3) * t184;
t111 = -t147 * pkin(1) - t151;
t139 = sin(pkin(8));
t140 = cos(pkin(8));
t152 = -g(2) * t184 - g(3) * t183;
t149 = qJDD(1) * pkin(1) + t152;
t171 = t140 * t111 + t139 * t149;
t76 = -t147 * pkin(2) + qJDD(1) * pkin(6) + t171;
t70 = -t145 * t137 + t143 * t76;
t46 = (-t113 + t158) * pkin(7) + t118 * pkin(3) - t70;
t120 = qJD(3) * pkin(3) - pkin(7) * t170;
t136 = t145 ^ 2;
t131 = t136 * t147;
t71 = t143 * t137 + t145 * t76;
t49 = -pkin(3) * t131 + t114 * pkin(7) - qJD(3) * t120 + t71;
t24 = t142 * t49 - t144 * t46;
t188 = qJ(5) * t95 + 0.2e1 * qJD(5) * t108 + t194 + t24;
t101 = t106 ^ 2;
t25 = t142 * t46 + t144 * t49;
t156 = t142 * t113 - t144 * t114;
t72 = -t108 * qJD(4) - t156;
t90 = t134 * pkin(4) - t108 * qJ(5);
t18 = -t101 * pkin(4) + t72 * qJ(5) - 0.2e1 * qJD(5) * t106 - t134 * t90 + t25;
t102 = t108 ^ 2;
t132 = t134 ^ 2;
t150 = (-qJD(4) + t134) * t108 - t156;
t63 = t73 + t95;
t34 = t142 * t150 - t144 * t63;
t187 = pkin(7) * t34;
t78 = -t132 - t101;
t50 = t142 * t78 - t191;
t186 = pkin(7) * t50;
t80 = t84 + t133;
t180 = t142 * t80;
t89 = -t102 - t132;
t64 = t144 * t89 - t180;
t185 = pkin(7) * t64;
t9 = t142 * t25 - t144 * t24;
t182 = t143 * t9;
t157 = -t139 * t111 + t140 * t149;
t75 = -qJDD(1) * pkin(2) - t147 * pkin(6) - t157;
t52 = -t114 * pkin(3) - pkin(7) * t131 + t120 * t170 + t75;
t181 = t142 * t52;
t179 = t144 * t52;
t178 = t144 * t80;
t177 = qJ(5) * t142;
t176 = qJ(5) * t144;
t175 = t134 * t142;
t174 = t134 * t144;
t173 = t143 * t118;
t119 = qJDD(3) - t124;
t172 = t145 * t119;
t35 = t142 * t63 + t144 * t150;
t16 = -t143 * t34 + t145 * t35;
t74 = -t101 - t102;
t165 = pkin(6) * t16 + pkin(1) * (t139 * t16 - t140 * t74) - pkin(2) * t74;
t51 = t144 * t78 + t192;
t28 = -t143 * t50 + t145 * t51;
t58 = (qJD(4) + t134) * t108 + t156;
t164 = pkin(1) * (t139 * t28 - t140 * t58) + pkin(6) * t28 - pkin(2) * t58;
t65 = -t142 * t89 - t178;
t37 = -t143 * t64 + t145 * t65;
t163 = pkin(1) * (t139 * t37 - t140 * t190) + pkin(6) * t37 - pkin(2) * t190;
t162 = -pkin(3) * t74 + pkin(7) * t35;
t161 = -pkin(3) * t58 + pkin(7) * t51;
t160 = -pkin(3) * t190 + pkin(7) * t65;
t10 = t142 * t24 + t144 * t25;
t41 = t143 * t70 + t145 * t71;
t153 = pkin(4) * t89 - t18;
t13 = -t73 * qJ(5) - t188;
t148 = t13 - t194;
t22 = -t72 * pkin(4) - t101 * qJ(5) + t108 * t90 + qJDD(5) + t52;
t146 = qJD(3) ^ 2;
t135 = t143 ^ 2;
t130 = t135 * t147;
t122 = -t131 - t146;
t121 = -t130 - t146;
t117 = t130 + t131;
t116 = (t135 + t136) * qJDD(1);
t115 = t129 - 0.2e1 * t159;
t112 = 0.2e1 * t158 + t166;
t92 = -t102 + t132;
t91 = t101 - t132;
t88 = -t143 * t121 - t172;
t87 = t145 * t122 - t173;
t82 = t102 - t101;
t57 = pkin(3) * t64;
t53 = pkin(4) * t63;
t48 = pkin(3) * t50;
t42 = (t143 * (-t106 * t144 + t108 * t142) + t145 * (-t106 * t142 - t108 * t144)) * t134;
t40 = -pkin(4) * t190 - qJ(5) * t80;
t39 = t143 * (t144 * t91 - t180) + t145 * (t142 * t91 + t178);
t38 = t143 * (-t142 * t92 - t191) + t145 * (t144 * t92 - t192);
t36 = t143 * t65 + t145 * t64;
t32 = pkin(3) * t34;
t30 = t143 * (-t108 * t175 + t144 * t73) + t145 * (t108 * t174 + t142 * t73);
t29 = t143 * (t106 * t174 - t142 * t72) + t145 * (t106 * t175 + t144 * t72);
t27 = t143 * t51 + t145 * t50;
t21 = -qJ(5) * t89 + t22;
t17 = -pkin(4) * t58 + qJ(5) * t78 - t22;
t15 = t143 * t35 + t145 * t34;
t14 = t143 * (-t142 * t190 - t144 * t58) + t145 * (-t142 * t58 + t144 * t190);
t12 = pkin(4) * t13;
t7 = (t63 + t73) * qJ(5) + t188;
t6 = -pkin(4) * t74 + qJ(5) * t150 + t18;
t5 = -pkin(4) * t22 + qJ(5) * t18;
t4 = -t142 * t13 + t144 * t18;
t3 = t144 * t13 + t142 * t18;
t2 = t145 * t10 - t182;
t1 = -t143 * t3 + t145 * t4;
t8 = [0, 0, 0, 0, 0, qJDD(1), t152, t151, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t140 * qJDD(1) - t139 * t147) + t157, pkin(1) * (-t139 * qJDD(1) - t140 * t147) - t171, 0, pkin(1) * (t139 * t171 + t140 * t157), (t113 + t158) * t143, t145 * t112 + t143 * t115, t173 + t145 * (-t130 + t146), (t114 - t159) * t145, t143 * (t131 - t146) + t172, 0, -t145 * t75 + pkin(2) * t115 + pkin(6) * t87 + pkin(1) * (t140 * t115 + t139 * t87), t143 * t75 - pkin(2) * t112 + pkin(6) * t88 + pkin(1) * (-t140 * t112 + t139 * t88), pkin(2) * t117 + pkin(6) * t116 + pkin(1) * (t139 * t116 + t140 * t117) + t41, -pkin(2) * t75 + pkin(6) * t41 + pkin(1) * (t139 * t41 - t140 * t75), t30, t14, t38, t29, t39, t42, t143 * (t181 - t186) + t145 * (t161 - t179) + t164, t143 * (t179 - t185) + t145 * (t160 + t181) + t163, t143 * (-t9 - t187) + t145 * (t10 + t162) + t165, -pkin(7) * t182 + t145 * (-pkin(3) * t52 + pkin(7) * t10) - pkin(2) * t52 + pkin(6) * t2 + pkin(1) * (t139 * t2 - t140 * t52), t30, t14, t38, t29, t39, t42, t143 * (-t142 * t17 + t176 * t193 - t186) + t145 * (t144 * t17 + t177 * t193 + t161) + t164, t143 * (-t142 * t40 + t144 * t21 - t185) + t145 * (t142 * t21 + t144 * t40 + t160) + t163, t143 * (-t142 * t6 + t144 * t7 - t187) + t145 * (t142 * t7 + t144 * t6 + t162) + t165, t143 * (-pkin(7) * t3 - t13 * t176 - t142 * t5) + t145 * (-pkin(3) * t22 + pkin(7) * t4 - t13 * t177 + t144 * t5) - pkin(2) * t22 + pkin(6) * t1 + pkin(1) * (t139 * t1 - t140 * t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, 0, 0, 0, 0, 0, 0, t145 * t118 + t143 * t122, -t143 * t119 + t145 * t121, 0, t143 * t71 - t145 * t70, 0, 0, 0, 0, 0, 0, t27, t36, t15, t143 * t10 + t145 * t9, 0, 0, 0, 0, 0, 0, t27, t36, t15, t143 * t4 + t145 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, t130 - t131, t166, t124, t129, qJDD(3), -t70, -t71, 0, 0, t84, t82, t63, -t84, t150, t133, -t24 + t48, t57 - t25, t32, pkin(3) * t9, t84, t82, t63, -t84, t150, t133, t148 + t48, t57 + t153, -t53 + t32, pkin(3) * t3 + t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t82, t63, -t84, t150, t133, -t24, -t25, 0, 0, t84, t82, t63, -t84, t150, t133, t148, t153, -t53, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t190, t74, t22;];
tauJ_reg = t8;
