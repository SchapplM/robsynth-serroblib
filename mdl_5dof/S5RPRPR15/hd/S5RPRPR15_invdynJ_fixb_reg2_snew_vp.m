% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR15
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR15_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:25
% EndTime: 2019-12-31 18:37:31
% DurationCPUTime: 2.06s
% Computational Cost: add. (6837->271), mult. (14452->370), div. (0->0), fcn. (9199->8), ass. (0->178)
t147 = sin(pkin(8));
t148 = cos(pkin(8));
t153 = cos(qJ(3));
t184 = qJD(1) * t153;
t125 = -t148 * qJD(3) + t147 * t184;
t127 = t147 * qJD(3) + t148 * t184;
t106 = t127 * t125;
t181 = qJD(1) * qJD(3);
t172 = t153 * t181;
t150 = sin(qJ(3));
t179 = t150 * qJDD(1);
t130 = t172 + t179;
t202 = -t106 + t130;
t209 = t147 * t202;
t208 = t148 * t202;
t149 = sin(qJ(5));
t128 = qJDD(5) + t130;
t152 = cos(qJ(5));
t100 = t152 * t125 + t149 * t127;
t102 = -t149 * t125 + t152 * t127;
t74 = t102 * t100;
t204 = t128 - t74;
t207 = t149 * t204;
t206 = t152 * t204;
t140 = t153 * qJDD(1);
t173 = t150 * t181;
t131 = t140 - t173;
t113 = t147 * qJDD(3) + t148 * t131;
t169 = -t148 * qJDD(3) + t147 * t131;
t66 = -t100 * qJD(5) + t152 * t113 - t149 * t169;
t183 = t150 * qJD(1);
t137 = qJD(5) + t183;
t89 = t137 * t100;
t205 = t66 - t89;
t156 = qJD(1) ^ 2;
t201 = pkin(6) + pkin(1);
t203 = t201 * t156;
t170 = t149 * t113 + t152 * t169;
t52 = (qJD(5) - t137) * t102 + t170;
t98 = t100 ^ 2;
t99 = t102 ^ 2;
t123 = t125 ^ 2;
t124 = t127 ^ 2;
t136 = t137 ^ 2;
t180 = qJD(2) * qJD(1);
t142 = 0.2e1 * t180;
t144 = qJDD(1) * qJ(2);
t151 = sin(qJ(1));
t154 = cos(qJ(1));
t168 = t154 * g(1) + t151 * g(2);
t162 = -t144 + t168;
t158 = t142 - t162;
t164 = -t131 + t173;
t165 = t130 + t172;
t80 = t165 * pkin(3) + t164 * qJ(4) + t158 - t203;
t171 = t151 * g(1) - t154 * g(2);
t167 = qJDD(2) - t171;
t157 = -t156 * qJ(2) + t167;
t114 = -t201 * qJDD(1) + t157;
t105 = t153 * g(3) - t150 * t114;
t155 = qJD(3) ^ 2;
t166 = pkin(3) * t150 - qJ(4) * t153;
t159 = t156 * t166;
t83 = -t155 * pkin(3) + qJDD(3) * qJ(4) - t150 * t159 - t105;
t40 = 0.2e1 * qJD(4) * t127 + t147 * t83 - t148 * t80;
t175 = t125 * t183;
t93 = -t113 - t175;
t33 = pkin(4) * t202 + t93 * pkin(7) - t40;
t161 = pkin(4) * t183 - t127 * pkin(7);
t41 = -0.2e1 * qJD(4) * t125 + t147 * t80 + t148 * t83;
t34 = -t123 * pkin(4) - pkin(7) * t169 - t161 * t183 + t41;
t14 = t149 * t34 - t152 * t33;
t15 = t149 * t33 + t152 * t34;
t7 = -t152 * t14 + t149 * t15;
t200 = t147 * t7;
t199 = t148 * t7;
t104 = t150 * g(3) + t153 * t114;
t82 = qJDD(3) * pkin(3) + t155 * qJ(4) - t153 * t159 - qJDD(4) + t104;
t198 = t147 * t82;
t95 = t106 + t130;
t197 = t147 * t95;
t196 = t148 * t82;
t195 = t148 * t95;
t58 = -pkin(4) * t169 + t123 * pkin(7) - t127 * t161 + t82;
t194 = t149 * t58;
t68 = t128 + t74;
t193 = t149 * t68;
t192 = t152 * t58;
t191 = t152 * t68;
t190 = qJDD(1) * pkin(1);
t189 = t137 * t149;
t188 = t137 * t152;
t146 = t153 ^ 2;
t187 = t146 * t156;
t176 = t150 * t156 * t153;
t186 = t150 * (qJDD(3) + t176);
t185 = t153 * (qJDD(3) - t176);
t178 = t150 * t74;
t177 = t150 * t106;
t174 = t127 * t183;
t8 = t149 * t14 + t152 * t15;
t26 = t147 * t40 + t148 * t41;
t163 = t147 * t41 - t148 * t40;
t79 = t153 * t104 - t150 * t105;
t160 = qJ(2) + t166;
t91 = -t169 + t174;
t145 = t150 ^ 2;
t141 = t145 * t156;
t133 = (t145 + t146) * qJDD(1);
t132 = t140 - 0.2e1 * t173;
t129 = 0.2e1 * t172 + t179;
t118 = -t157 + t190;
t117 = -t124 - t141;
t116 = -t124 + t141;
t115 = t123 - t141;
t111 = t162 - 0.2e1 * t180 + t203;
t109 = -t186 + t153 * (-t155 - t187);
t108 = t150 * (-t141 - t155) + t185;
t103 = -t141 - t123;
t92 = t113 - t175;
t90 = t169 + t174;
t87 = -t123 - t124;
t86 = -t99 + t136;
t85 = t98 - t136;
t84 = -t99 - t136;
t76 = -t147 * t117 - t195;
t75 = t148 * t117 - t197;
t73 = t99 - t98;
t72 = -t136 - t98;
t71 = t148 * t103 - t209;
t70 = t147 * t103 + t208;
t65 = -t102 * qJD(5) - t170;
t64 = -t147 * t93 + t148 * t91;
t62 = (-t100 * t152 + t102 * t149) * t137;
t61 = (-t100 * t149 - t102 * t152) * t137;
t60 = -t98 - t99;
t59 = t150 * t76 - t153 * t92;
t57 = t150 * t71 - t153 * t90;
t56 = t66 + t89;
t51 = (qJD(5) + t137) * t102 + t170;
t50 = t152 * t85 - t193;
t49 = -t149 * t86 + t206;
t48 = t149 * t85 + t191;
t47 = t152 * t86 + t207;
t46 = -t102 * t189 + t152 * t66;
t45 = t102 * t188 + t149 * t66;
t44 = t100 * t188 - t149 * t65;
t43 = t100 * t189 + t152 * t65;
t42 = t150 * t64 - t153 * t87;
t39 = -t149 * t84 - t191;
t38 = t152 * t84 - t193;
t36 = t152 * t72 - t207;
t35 = t149 * t72 + t206;
t31 = t149 * t56 - t152 * t52;
t30 = -t149 * t205 - t152 * t51;
t29 = -t149 * t52 - t152 * t56;
t28 = -t149 * t51 + t152 * t205;
t27 = -pkin(7) * t38 - t192;
t24 = -t147 * t38 + t148 * t39;
t23 = t147 * t39 + t148 * t38;
t22 = -pkin(7) * t35 - t194;
t21 = -t147 * t35 + t148 * t36;
t20 = t147 * t36 + t148 * t35;
t19 = t150 * t26 + t153 * t82;
t18 = -pkin(4) * t205 + pkin(7) * t39 - t194;
t17 = t150 * t24 - t153 * t205;
t16 = -pkin(4) * t51 + pkin(7) * t36 + t192;
t12 = t150 * t21 - t153 * t51;
t11 = -t147 * t29 + t148 * t31;
t10 = t147 * t31 + t148 * t29;
t9 = t150 * t11 - t153 * t60;
t6 = pkin(4) * t58 + pkin(7) * t8;
t5 = -pkin(7) * t29 - t7;
t4 = -pkin(4) * t60 + pkin(7) * t31 + t8;
t3 = t148 * t8 - t200;
t2 = t147 * t8 + t199;
t1 = t150 * t3 + t153 * t58;
t13 = [0, 0, 0, 0, 0, qJDD(1), t171, t168, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t167 - 0.2e1 * t190, t142 + 0.2e1 * t144 - t168, pkin(1) * t118 + qJ(2) * (-t156 * pkin(1) + t158), -t164 * t153, -t153 * t129 - t150 * t132, t185 - t150 * (t155 - t187), t165 * t150, t153 * (t141 - t155) - t186, 0, qJ(2) * t129 - t201 * t108 - t150 * t111, qJ(2) * t132 - t201 * t109 - t153 * t111, qJ(2) * (-t141 - t187) + t201 * t133 - t79, -qJ(2) * t111 - t201 * t79, t153 * (t148 * t113 - t147 * t174) + t177, t153 * (-t147 * t92 - t148 * t90) - t150 * (-t124 + t123), t153 * (-t147 * t116 + t208) - t150 * t93, t153 * (t147 * t169 + t148 * t175) - t177, t153 * (t148 * t115 - t197) + t150 * t91, (t130 + (-t125 * t148 + t127 * t147) * t184) * t150, t153 * (-qJ(4) * t70 - t198) - t150 * (-pkin(3) * t70 + t40) + qJ(2) * t70 - t201 * t57, t153 * (-qJ(4) * t75 - t196) - t150 * (-pkin(3) * t75 + t41) + qJ(2) * t75 - t201 * t59, -t153 * t163 + t160 * (t147 * t91 + t148 * t93) - t201 * t42, t160 * t163 - t201 * t19, t153 * (-t147 * t45 + t148 * t46) + t178, t153 * (-t147 * t28 + t148 * t30) + t150 * t73, t153 * (-t147 * t47 + t148 * t49) + t150 * t56, t153 * (-t147 * t43 + t148 * t44) - t178, t153 * (-t147 * t48 + t148 * t50) - t150 * t52, t153 * (-t147 * t61 + t148 * t62) + t150 * t128, t153 * (-qJ(4) * t20 - t147 * t16 + t148 * t22) - t150 * (-pkin(3) * t20 - pkin(4) * t35 + t14) + qJ(2) * t20 - t201 * t12, t153 * (-qJ(4) * t23 - t147 * t18 + t148 * t27) - t150 * (-pkin(3) * t23 - pkin(4) * t38 + t15) + qJ(2) * t23 - t201 * t17, t153 * (-qJ(4) * t10 - t147 * t4 + t148 * t5) - t150 * (-pkin(3) * t10 - pkin(4) * t29) + qJ(2) * t10 - t201 * t9, t153 * (-pkin(7) * t199 - qJ(4) * t2 - t147 * t6) - t150 * (-pkin(3) * t2 - pkin(4) * t7) + qJ(2) * t2 - t201 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t156, -t118, 0, 0, 0, 0, 0, 0, t108, t109, -t133, t79, 0, 0, 0, 0, 0, 0, t57, t59, t42, t19, 0, 0, 0, 0, 0, 0, t12, t17, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, -t141 + t187, t140, -t176, -t179, qJDD(3), t104, t105, 0, 0, t147 * t113 + t148 * t174, -t147 * t90 + t148 * t92, t148 * t116 + t209, t147 * t175 - t148 * t169, t147 * t115 + t195, (-t125 * t147 - t127 * t148) * t183, -pkin(3) * t90 + qJ(4) * t71 + t196, -pkin(3) * t92 + qJ(4) * t76 - t198, -pkin(3) * t87 + qJ(4) * t64 + t26, pkin(3) * t82 + qJ(4) * t26, t147 * t46 + t148 * t45, t147 * t30 + t148 * t28, t147 * t49 + t148 * t47, t147 * t44 + t148 * t43, t147 * t50 + t148 * t48, t147 * t62 + t148 * t61, -pkin(3) * t51 + qJ(4) * t21 + t147 * t22 + t148 * t16, -pkin(3) * t205 + qJ(4) * t24 + t147 * t27 + t148 * t18, -pkin(3) * t60 + qJ(4) * t11 + t147 * t5 + t148 * t4, pkin(3) * t58 - pkin(7) * t200 + qJ(4) * t3 + t148 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t92, t87, -t82, 0, 0, 0, 0, 0, 0, t51, t205, t60, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t73, t56, -t74, -t52, t128, -t14, -t15, 0, 0;];
tauJ_reg = t13;
