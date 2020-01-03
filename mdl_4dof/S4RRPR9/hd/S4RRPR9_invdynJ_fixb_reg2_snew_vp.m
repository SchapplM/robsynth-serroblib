% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRPR9
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
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRPR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:04
% EndTime: 2019-12-31 17:10:09
% DurationCPUTime: 1.68s
% Computational Cost: add. (5182->258), mult. (11402->361), div. (0->0), fcn. (7593->8), ass. (0->159)
t149 = qJD(1) ^ 2;
t143 = sin(qJ(4));
t146 = cos(qJ(4));
t144 = sin(qJ(2));
t133 = t144 * qJDD(1);
t147 = cos(qJ(2));
t166 = qJD(1) * qJD(2);
t162 = t147 * t166;
t121 = t133 + t162;
t141 = sin(pkin(7));
t142 = cos(pkin(7));
t158 = -t142 * qJDD(2) + t141 * t121;
t169 = qJD(1) * t144;
t115 = -t142 * qJD(2) + t141 * t169;
t117 = t141 * qJD(2) + t142 * t169;
t89 = t146 * t115 + t143 * t117;
t99 = t141 * qJDD(2) + t142 * t121;
t58 = -t89 * qJD(4) - t143 * t158 + t146 * t99;
t168 = t147 * qJD(1);
t128 = -qJD(4) + t168;
t79 = t89 * t128;
t194 = t58 + t79;
t104 = t115 * t168;
t83 = t104 - t99;
t130 = t144 * t166;
t165 = t147 * qJDD(1);
t122 = -t130 + t165;
t174 = t117 * t115;
t151 = -t122 - t174;
t193 = t141 * t151;
t192 = t142 * t151;
t118 = -qJDD(4) + t122;
t91 = -t143 * t115 + t146 * t117;
t184 = t91 * t89;
t150 = -t118 - t184;
t191 = t143 * t150;
t190 = t146 * t150;
t105 = t117 * t168;
t81 = -t158 - t105;
t160 = t143 * t99 + t146 * t158;
t44 = (qJD(4) + t128) * t91 + t160;
t87 = t89 ^ 2;
t88 = t91 ^ 2;
t113 = t115 ^ 2;
t114 = t117 ^ 2;
t126 = t128 ^ 2;
t189 = qJD(2) ^ 2;
t145 = sin(qJ(1));
t148 = cos(qJ(1));
t161 = t145 * g(1) - t148 * g(2);
t109 = qJDD(1) * pkin(1) + t149 * pkin(5) + t161;
t154 = t121 + t162;
t70 = -t154 * qJ(3) + (-t122 + t130) * pkin(2) - t109;
t156 = t148 * g(1) + t145 * g(2);
t175 = qJDD(1) * pkin(5);
t110 = -t149 * pkin(1) - t156 + t175;
t155 = -t147 * pkin(2) - t144 * qJ(3);
t157 = t149 * t155 + t110;
t186 = t144 * g(3);
t74 = -t189 * pkin(2) + qJDD(2) * qJ(3) + t157 * t147 - t186;
t49 = 0.2e1 * qJD(3) * t117 + t141 * t74 - t142 * t70;
t28 = t151 * pkin(3) + t83 * pkin(6) - t49;
t100 = -pkin(3) * t168 - t117 * pkin(6);
t50 = -0.2e1 * qJD(3) * t115 + t141 * t70 + t142 * t74;
t29 = -t113 * pkin(3) - pkin(6) * t158 + t100 * t168 + t50;
t12 = t143 * t29 - t146 * t28;
t13 = t143 * t28 + t146 * t29;
t6 = -t146 * t12 + t143 * t13;
t188 = t141 * t6;
t187 = t142 * t6;
t185 = t147 * g(3);
t73 = -qJDD(2) * pkin(2) - t189 * qJ(3) + t157 * t144 + qJDD(3) + t185;
t183 = t141 * t73;
t84 = t122 - t174;
t182 = t141 * t84;
t181 = t142 * t73;
t180 = t142 * t84;
t51 = pkin(3) * t158 - t113 * pkin(6) + t117 * t100 + t73;
t179 = t143 * t51;
t59 = t118 - t184;
t178 = t143 * t59;
t177 = t146 * t51;
t176 = t146 * t59;
t173 = t128 * t143;
t172 = t128 * t146;
t127 = t147 * t149 * t144;
t171 = t144 * (qJDD(2) + t127);
t170 = t147 * (qJDD(2) - t127);
t164 = t147 * t184;
t163 = t147 * t174;
t7 = t143 * t12 + t146 * t13;
t26 = t141 * t49 + t142 * t50;
t95 = t144 * t110 + t185;
t96 = t147 * t110 - t186;
t159 = t144 * t95 + t147 * t96;
t153 = t141 * t50 - t142 * t49;
t152 = -pkin(1) + t155;
t139 = t147 ^ 2;
t138 = t144 ^ 2;
t135 = t139 * t149;
t134 = t138 * t149;
t123 = -0.2e1 * t130 + t165;
t120 = t133 + 0.2e1 * t162;
t111 = t147 * t122;
t103 = -t114 - t135;
t102 = -t114 + t135;
t101 = t113 - t135;
t92 = -t135 - t113;
t82 = t104 + t99;
t80 = -t105 + t158;
t77 = -t113 - t114;
t76 = -t88 + t126;
t75 = t87 - t126;
t72 = -t88 - t126;
t67 = -t141 * t103 + t180;
t66 = t142 * t103 + t182;
t65 = t88 - t87;
t64 = -t126 - t87;
t63 = t142 * t92 - t193;
t62 = t141 * t92 + t192;
t57 = -t91 * qJD(4) - t160;
t56 = -t141 * t83 + t142 * t81;
t54 = (-t143 * t91 + t146 * t89) * t128;
t53 = (t143 * t89 + t146 * t91) * t128;
t52 = -t87 - t88;
t48 = t58 - t79;
t43 = (qJD(4) - t128) * t91 + t160;
t41 = t146 * t75 + t178;
t40 = -t143 * t76 + t190;
t39 = t143 * t75 - t176;
t38 = t146 * t76 + t191;
t37 = t146 * t58 + t91 * t173;
t36 = t143 * t58 - t91 * t172;
t35 = -t143 * t57 - t89 * t172;
t34 = t146 * t57 - t89 * t173;
t33 = -t143 * t72 + t176;
t32 = t146 * t72 + t178;
t31 = t146 * t64 - t191;
t30 = t143 * t64 + t190;
t24 = t143 * t48 - t146 * t44;
t23 = -t143 * t194 - t146 * t43;
t22 = -t143 * t44 - t146 * t48;
t21 = -t143 * t43 + t146 * t194;
t20 = -pkin(6) * t32 + t177;
t19 = -t141 * t32 + t142 * t33;
t18 = t141 * t33 + t142 * t32;
t17 = -pkin(6) * t30 + t179;
t16 = -t141 * t30 + t142 * t31;
t15 = t141 * t31 + t142 * t30;
t14 = -pkin(3) * t194 + pkin(6) * t33 + t179;
t10 = -pkin(3) * t43 + pkin(6) * t31 - t177;
t9 = -t141 * t22 + t142 * t24;
t8 = t141 * t24 + t142 * t22;
t5 = -pkin(3) * t51 + pkin(6) * t7;
t4 = -pkin(6) * t22 - t6;
t3 = -pkin(3) * t52 + pkin(6) * t24 + t7;
t2 = t142 * t7 - t188;
t1 = t141 * t7 + t187;
t11 = [0, 0, 0, 0, 0, qJDD(1), t161, t156, 0, 0, t154 * t144, t147 * t120 + t144 * t123, t171 + t147 * (-t134 + t189), -t144 * t162 + t111, t144 * (t135 - t189) + t170, 0, t147 * t109 + pkin(1) * t123 + pkin(5) * (t147 * (-t135 - t189) - t171), -t144 * t109 - pkin(1) * t120 + pkin(5) * (-t170 - t144 * (-t134 - t189)), pkin(1) * (t134 + t135) + (t138 + t139) * t175 + t159, pkin(1) * t109 + pkin(5) * t159, t144 * (t141 * t105 + t142 * t99) - t163, t144 * (-t141 * t82 - t142 * t80) + t147 * (-t114 + t113), t144 * (-t141 * t102 + t192) + t147 * t83, t144 * (-t142 * t104 + t141 * t158) + t163, t144 * (t142 * t101 + t182) - t147 * t81, t111 + t144 * (t115 * t142 - t117 * t141) * t168, t144 * (-qJ(3) * t62 + t183) + t147 * (-pkin(2) * t62 + t49) - pkin(1) * t62 + pkin(5) * (t144 * t80 + t147 * t63), t144 * (-qJ(3) * t66 + t181) + t147 * (-pkin(2) * t66 + t50) - pkin(1) * t66 + pkin(5) * (t144 * t82 + t147 * t67), -t144 * t153 + pkin(5) * (t144 * t77 + t147 * t56) + t152 * (t141 * t81 + t142 * t83), pkin(5) * (t144 * t73 + t147 * t26) + t152 * t153, t144 * (-t141 * t36 + t142 * t37) - t164, t144 * (-t141 * t21 + t142 * t23) - t147 * t65, t144 * (-t141 * t38 + t142 * t40) - t147 * t48, t144 * (-t141 * t34 + t142 * t35) + t164, t144 * (-t141 * t39 + t142 * t41) + t147 * t44, t144 * (-t141 * t53 + t142 * t54) + t147 * t118, t144 * (-qJ(3) * t15 - t141 * t10 + t142 * t17) + t147 * (-pkin(2) * t15 - pkin(3) * t30 + t12) - pkin(1) * t15 + pkin(5) * (t144 * t43 + t147 * t16), t144 * (-qJ(3) * t18 - t141 * t14 + t142 * t20) + t147 * (-pkin(2) * t18 - pkin(3) * t32 + t13) - pkin(1) * t18 + pkin(5) * (t144 * t194 + t147 * t19), t144 * (-qJ(3) * t8 - t141 * t3 + t142 * t4) + t147 * (-pkin(2) * t8 - pkin(3) * t22) - pkin(1) * t8 + pkin(5) * (t144 * t52 + t147 * t9), t144 * (-pkin(6) * t187 - qJ(3) * t1 - t141 * t5) + t147 * (-pkin(2) * t1 - pkin(3) * t6) - pkin(1) * t1 + pkin(5) * (t144 * t51 + t147 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, t134 - t135, t133, t127, t165, qJDD(2), -t95, -t96, 0, 0, -t142 * t105 + t141 * t99, -t141 * t80 + t142 * t82, t142 * t102 + t193, -t141 * t104 - t142 * t158, t141 * t101 - t180, (t115 * t141 + t117 * t142) * t168, -pkin(2) * t80 + qJ(3) * t63 - t181, -pkin(2) * t82 + qJ(3) * t67 + t183, -pkin(2) * t77 + qJ(3) * t56 + t26, -pkin(2) * t73 + qJ(3) * t26, t141 * t37 + t142 * t36, t141 * t23 + t142 * t21, t141 * t40 + t142 * t38, t141 * t35 + t142 * t34, t141 * t41 + t142 * t39, t141 * t54 + t142 * t53, -pkin(2) * t43 + qJ(3) * t16 + t142 * t10 + t141 * t17, -pkin(2) * t194 + qJ(3) * t19 + t142 * t14 + t141 * t20, -pkin(2) * t52 + qJ(3) * t9 + t141 * t4 + t142 * t3, -pkin(2) * t51 - pkin(6) * t188 + qJ(3) * t2 + t142 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t82, t77, t73, 0, 0, 0, 0, 0, 0, t43, t194, t52, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, t65, t48, -t184, -t44, -t118, -t12, -t13, 0, 0;];
tauJ_reg = t11;
