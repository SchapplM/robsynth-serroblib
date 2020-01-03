% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRRR3
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:37
% EndTime: 2019-12-31 17:24:42
% DurationCPUTime: 1.76s
% Computational Cost: add. (6475->274), mult. (14293->377), div. (0->0), fcn. (9888->8), ass. (0->171)
t141 = sin(qJ(4));
t137 = qJDD(2) + qJDD(3);
t129 = qJDD(4) + t137;
t142 = sin(qJ(3));
t146 = cos(qJ(3));
t147 = cos(qJ(2));
t143 = sin(qJ(2));
t169 = qJD(1) * t143;
t113 = -t146 * t147 * qJD(1) + t142 * t169;
t115 = (t147 * t142 + t143 * t146) * qJD(1);
t145 = cos(qJ(4));
t94 = t145 * t113 + t141 * t115;
t96 = -t141 * t113 + t145 * t115;
t70 = t96 * t94;
t193 = -t70 + t129;
t197 = t141 * t193;
t102 = t115 * t113;
t191 = -t102 + t137;
t196 = t142 * t191;
t195 = t145 * t193;
t194 = t146 * t191;
t138 = qJD(2) + qJD(3);
t109 = t138 * t113;
t131 = t143 * qJDD(1);
t166 = qJD(1) * qJD(2);
t164 = t147 * t166;
t121 = t131 + t164;
t132 = t147 * qJDD(1);
t165 = t143 * t166;
t122 = t132 - t165;
t156 = t146 * t121 + t142 * t122;
t85 = -t113 * qJD(3) + t156;
t192 = -t85 - t109;
t92 = t94 ^ 2;
t93 = t96 ^ 2;
t111 = t113 ^ 2;
t112 = t115 ^ 2;
t130 = qJD(4) + t138;
t128 = t130 ^ 2;
t136 = t138 ^ 2;
t190 = cos(qJ(1));
t149 = qJD(1) ^ 2;
t171 = t143 * t149;
t144 = sin(qJ(1));
t155 = t190 * g(1) + t144 * g(2);
t178 = qJDD(1) * pkin(5);
t117 = -t149 * pkin(1) - t155 + t178;
t173 = t143 * t117;
t80 = qJDD(2) * pkin(2) - t121 * pkin(6) - t173 + (pkin(2) * t171 + pkin(6) * t166 - g(3)) * t147;
t106 = -t143 * g(3) + t147 * t117;
t140 = t147 ^ 2;
t134 = t140 * t149;
t153 = qJD(2) * pkin(2) - pkin(6) * t169;
t81 = -pkin(2) * t134 + t122 * pkin(6) - qJD(2) * t153 + t106;
t59 = t142 * t81 - t146 * t80;
t30 = t191 * pkin(3) + t192 * pkin(7) - t59;
t158 = t138 * pkin(3) - t115 * pkin(7);
t60 = t142 * t80 + t146 * t81;
t160 = t142 * t121 - t146 * t122;
t84 = -t115 * qJD(3) - t160;
t31 = -t111 * pkin(3) + t84 * pkin(7) - t138 * t158 + t60;
t14 = t141 * t31 - t145 * t30;
t15 = t141 * t30 + t145 * t31;
t7 = -t145 * t14 + t141 * t15;
t189 = t142 * t7;
t188 = t146 * t7;
t163 = t144 * g(1) - t190 * g(2);
t152 = qJDD(1) * pkin(1) + t163;
t87 = t122 * pkin(2) - t153 * t169 + (pkin(6) * t140 + pkin(5)) * t149 + t152;
t45 = t84 * pkin(3) + t111 * pkin(7) - t115 * t158 + t87;
t187 = t141 * t45;
t65 = t70 + t129;
t186 = t141 * t65;
t185 = t142 * t87;
t99 = t102 + t137;
t184 = t142 * t99;
t27 = t142 * t60 - t146 * t59;
t183 = t143 * t27;
t182 = t145 * t45;
t181 = t145 * t65;
t180 = t146 * t87;
t179 = t146 * t99;
t177 = t130 * t141;
t176 = t130 * t145;
t175 = t138 * t142;
t174 = t138 * t146;
t127 = t147 * t171;
t172 = t143 * (qJDD(2) + t127);
t170 = t147 * (qJDD(2) - t127);
t168 = qJD(3) + t138;
t167 = qJD(4) + t130;
t8 = t141 * t14 + t145 * t15;
t162 = t141 * t85 - t145 * t84;
t28 = t142 * t59 + t146 * t60;
t105 = t147 * g(3) + t173;
t161 = t143 * t105 + t147 * t106;
t63 = -t128 - t92;
t43 = t141 * t63 + t195;
t159 = pkin(3) * t43 - t14;
t157 = t141 * t84 + t145 * t85;
t83 = -t93 - t128;
t51 = t145 * t83 - t186;
t154 = pkin(3) * t51 - t15;
t151 = (-qJD(4) + t130) * t96 - t162;
t47 = -t94 * qJD(4) + t157;
t150 = (-qJD(3) + t138) * t115 - t160;
t148 = qJD(2) ^ 2;
t139 = t143 ^ 2;
t133 = t139 * t149;
t123 = t132 - 0.2e1 * t165;
t120 = t131 + 0.2e1 * t164;
t116 = t149 * pkin(5) + t152;
t108 = -t112 + t136;
t107 = t111 - t136;
t104 = -t112 - t136;
t101 = t112 - t111;
t97 = -t136 - t111;
t90 = t130 * t94;
t89 = -t93 + t128;
t88 = t92 - t128;
t86 = -t111 - t112;
t78 = -t142 * t104 - t179;
t77 = t146 * t104 - t184;
t75 = -t109 + t85;
t74 = -t168 * t113 + t156;
t71 = t168 * t115 + t160;
t69 = t93 - t92;
t68 = t146 * t97 - t196;
t67 = t142 * t97 + t194;
t62 = (t141 * t96 - t145 * t94) * t130;
t61 = (-t141 * t94 - t145 * t96) * t130;
t58 = -t92 - t93;
t56 = t145 * t88 - t186;
t55 = -t141 * t89 + t195;
t54 = t141 * t88 + t181;
t53 = t145 * t89 + t197;
t52 = -t141 * t83 - t181;
t49 = -t142 * t192 + t146 * t150;
t48 = t142 * t150 + t146 * t192;
t46 = -t96 * qJD(4) - t162;
t44 = t145 * t63 - t197;
t41 = -t167 * t94 + t157;
t40 = t47 + t90;
t39 = t47 - t90;
t36 = t167 * t96 + t162;
t35 = t145 * t47 - t96 * t177;
t34 = t141 * t47 + t96 * t176;
t33 = -t141 * t46 + t94 * t176;
t32 = t145 * t46 + t94 * t177;
t26 = -t142 * t51 + t146 * t52;
t25 = t142 * t52 + t146 * t51;
t24 = -pkin(7) * t51 - t182;
t23 = -pkin(7) * t43 - t187;
t22 = -t142 * t43 + t146 * t44;
t21 = t142 * t44 + t146 * t43;
t20 = t141 * t40 + t145 * t151;
t19 = -t141 * t39 - t145 * t36;
t18 = t141 * t151 - t145 * t40;
t17 = -t141 * t36 + t145 * t39;
t16 = pkin(3) * t18;
t12 = -pkin(3) * t41 + pkin(7) * t52 - t187;
t11 = -pkin(3) * t36 + pkin(7) * t44 + t182;
t10 = -t142 * t18 + t146 * t20;
t9 = t142 * t20 + t146 * t18;
t6 = pkin(3) * t7;
t5 = pkin(3) * t45 + pkin(7) * t8;
t4 = -pkin(7) * t18 - t7;
t3 = -pkin(3) * t58 + pkin(7) * t20 + t8;
t2 = t146 * t8 - t189;
t1 = t142 * t8 + t188;
t13 = [0, 0, 0, 0, 0, qJDD(1), t163, t155, 0, 0, (t121 + t164) * t143, t147 * t120 + t143 * t123, t172 + t147 * (-t133 + t148), (t122 - t165) * t147, t143 * (t134 - t148) + t170, 0, t147 * t116 + pkin(1) * t123 + pkin(5) * (t147 * (-t134 - t148) - t172), -t143 * t116 - pkin(1) * t120 + pkin(5) * (-t170 - t143 * (-t133 - t148)), pkin(1) * (t133 + t134) + (t139 + t140) * t178 + t161, pkin(1) * t116 + pkin(5) * t161, t143 * (-t115 * t175 + t146 * t85) + t147 * (t115 * t174 + t142 * t85), t143 * (-t142 * t75 - t146 * t71) + t147 * (-t142 * t71 + t146 * t75), t143 * (-t142 * t108 + t194) + t147 * (t146 * t108 + t196), t143 * (t113 * t174 - t142 * t84) + t147 * (t113 * t175 + t146 * t84), t143 * (t146 * t107 - t184) + t147 * (t142 * t107 + t179), (t143 * (-t113 * t146 + t115 * t142) + t147 * (-t113 * t142 - t115 * t146)) * t138, t143 * (-pkin(6) * t67 - t185) + t147 * (-pkin(2) * t71 + pkin(6) * t68 + t180) - pkin(1) * t71 + pkin(5) * (-t143 * t67 + t147 * t68), t143 * (-pkin(6) * t77 - t180) + t147 * (-pkin(2) * t74 + pkin(6) * t78 - t185) - pkin(1) * t74 + pkin(5) * (-t143 * t77 + t147 * t78), t143 * (-pkin(6) * t48 - t27) + t147 * (-pkin(2) * t86 + pkin(6) * t49 + t28) - pkin(1) * t86 + pkin(5) * (-t143 * t48 + t147 * t49), -pkin(6) * t183 + t147 * (pkin(2) * t87 + pkin(6) * t28) + pkin(1) * t87 + pkin(5) * (t147 * t28 - t183), t143 * (-t142 * t34 + t146 * t35) + t147 * (t142 * t35 + t146 * t34), t143 * (-t142 * t17 + t146 * t19) + t147 * (t142 * t19 + t146 * t17), t143 * (-t142 * t53 + t146 * t55) + t147 * (t142 * t55 + t146 * t53), t143 * (-t142 * t32 + t146 * t33) + t147 * (t142 * t33 + t146 * t32), t143 * (-t142 * t54 + t146 * t56) + t147 * (t142 * t56 + t146 * t54), t143 * (-t142 * t61 + t146 * t62) + t147 * (t142 * t62 + t146 * t61), t143 * (-pkin(6) * t21 - t142 * t11 + t146 * t23) + t147 * (-pkin(2) * t36 + pkin(6) * t22 + t146 * t11 + t142 * t23) - pkin(1) * t36 + pkin(5) * (-t143 * t21 + t147 * t22), t143 * (-pkin(6) * t25 - t142 * t12 + t146 * t24) + t147 * (-pkin(2) * t41 + pkin(6) * t26 + t146 * t12 + t142 * t24) - pkin(1) * t41 + pkin(5) * (-t143 * t25 + t147 * t26), t143 * (-pkin(6) * t9 - t142 * t3 + t146 * t4) + t147 * (-pkin(2) * t58 + pkin(6) * t10 + t142 * t4 + t146 * t3) - pkin(1) * t58 + pkin(5) * (t147 * t10 - t143 * t9), t143 * (-pkin(6) * t1 - pkin(7) * t188 - t142 * t5) + t147 * (pkin(2) * t45 + pkin(6) * t2 - pkin(7) * t189 + t146 * t5) + pkin(1) * t45 + pkin(5) * (-t143 * t1 + t147 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, t133 - t134, t131, t127, t132, qJDD(2), -t105, -t106, 0, 0, t102, t101, -t192, -t102, t150, t137, pkin(2) * t67 - t59, pkin(2) * t77 - t60, pkin(2) * t48, pkin(2) * t27, t70, t69, t40, -t70, t151, t129, pkin(2) * t21 + t159, pkin(2) * t25 + t154, pkin(2) * t9 + t16, pkin(2) * t1 + t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t101, -t192, -t102, t150, t137, -t59, -t60, 0, 0, t70, t69, t40, -t70, t151, t129, t159, t154, t16, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t69, t40, -t70, t151, t129, -t14, -t15, 0, 0;];
tauJ_reg = t13;
