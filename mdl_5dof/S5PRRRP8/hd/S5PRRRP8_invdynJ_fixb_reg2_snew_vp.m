% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRP8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRP8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:37
% EndTime: 2019-12-05 17:00:45
% DurationCPUTime: 2.24s
% Computational Cost: add. (3919->253), mult. (7623->323), div. (0->0), fcn. (5323->10), ass. (0->173)
t142 = cos(qJ(4));
t139 = sin(qJ(4));
t140 = sin(qJ(3));
t181 = qJD(2) * qJD(3);
t126 = t140 * t181;
t143 = cos(qJ(3));
t179 = t143 * qJDD(2);
t109 = -t126 + t179;
t102 = -qJDD(4) + t109;
t183 = qJD(2) * t140;
t103 = -t142 * qJD(3) + t139 * t183;
t105 = t139 * qJD(3) + t142 * t183;
t189 = t105 * t103;
t153 = t102 - t189;
t203 = t139 * t153;
t101 = t105 ^ 2;
t123 = t143 * qJD(2) - qJD(4);
t210 = t123 ^ 2;
t219 = -t101 - t210;
t34 = -t142 * t219 - t203;
t258 = pkin(2) * t34;
t257 = pkin(3) * t34;
t256 = pkin(8) * t34;
t197 = t142 * t153;
t41 = -t139 * t219 + t197;
t255 = pkin(8) * t41;
t254 = t140 * t41;
t253 = t143 * t41;
t144 = cos(qJ(2));
t252 = t144 * t34;
t62 = t102 + t189;
t196 = t142 * t62;
t211 = t103 ^ 2;
t214 = -t210 - t211;
t223 = t139 * t214 - t196;
t175 = t143 * t181;
t180 = t140 * qJDD(2);
t108 = t175 + t180;
t170 = t142 * qJDD(3) - t108 * t139;
t157 = t105 * qJD(4) - t170;
t91 = t105 * t123;
t220 = t157 - t91;
t202 = t139 * t62;
t222 = t142 * t214 + t202;
t237 = t140 * t220 + t143 * t222;
t251 = -pkin(2) * t223 + pkin(7) * t237;
t49 = t157 + t91;
t85 = t211 - t210;
t250 = t140 * (t142 * t85 + t203) + t143 * t49;
t159 = -t139 * qJDD(3) - t142 * t108;
t152 = -t103 * qJD(4) - t159;
t190 = t103 * t123;
t215 = t152 + t190;
t205 = t139 * t215;
t218 = t101 - t211;
t249 = t140 * (t142 * t220 + t205) + t143 * t218;
t135 = sin(pkin(5));
t136 = cos(pkin(5));
t141 = sin(qJ(2));
t248 = t136 * (t140 * t222 - t143 * t220) + (t141 * t237 - t144 * t223) * t135;
t246 = pkin(3) * t223;
t245 = pkin(8) * t222;
t244 = pkin(8) * t223;
t86 = -t101 + t210;
t243 = t142 * t86 - t202;
t242 = qJ(5) * t215;
t239 = t139 * t85 - t197;
t216 = t152 - t190;
t235 = t140 * (-t139 * t86 - t196) - t143 * t216;
t217 = t101 + t211;
t234 = pkin(3) * t217;
t145 = qJD(2) ^ 2;
t165 = -t143 * pkin(3) - t140 * pkin(8);
t192 = sin(pkin(9));
t193 = cos(pkin(9));
t113 = -t193 * g(1) - t192 * g(2);
t154 = t192 * g(1) - t193 * g(2);
t151 = t136 * t154;
t184 = -g(3) + qJDD(1);
t224 = t135 * t184 + t151;
t67 = t144 * t113 + t224 * t141;
t57 = -t145 * pkin(2) + qJDD(2) * pkin(7) + t67;
t171 = t145 * t165 + t57;
t209 = qJD(3) ^ 2;
t147 = -t135 * t154 + t136 * t184;
t83 = t143 * t147;
t29 = -qJDD(3) * pkin(3) - t209 * pkin(8) + t171 * t140 - t83;
t233 = pkin(4) * t157 - t242 + t29;
t231 = t140 * t217;
t227 = t143 * t217;
t221 = -t139 * t220 + t142 * t215;
t146 = t140 * t147;
t30 = -t209 * pkin(3) + qJDD(3) * pkin(8) + t171 * t143 + t146;
t161 = -t109 + t126;
t162 = t108 + t175;
t166 = t141 * t113 - t224 * t144;
t56 = -qJDD(2) * pkin(2) - t145 * pkin(7) + t166;
t33 = t161 * pkin(3) - t162 * pkin(8) + t56;
t16 = t139 * t33 + t142 * t30;
t74 = pkin(4) * t103 - qJ(5) * t105;
t173 = t102 * qJ(5) + t103 * t74 - t16;
t212 = -(t210 + t219) * pkin(4) - qJ(5) * t153 - t173;
t208 = pkin(4) * t142;
t207 = t139 * t29;
t204 = t139 * t216;
t200 = t142 * t29;
t198 = t142 * t216;
t191 = qJ(5) * t142;
t188 = t123 * t139;
t187 = t123 * t142;
t122 = t140 * t145 * t143;
t114 = qJDD(3) + t122;
t186 = t140 * t114;
t115 = qJDD(3) - t122;
t185 = t143 * t115;
t182 = qJD(5) * t123;
t177 = t103 * t187;
t176 = t143 * t189;
t174 = -qJ(5) * t139 - pkin(3);
t15 = t139 * t30 - t142 * t33;
t6 = t139 * t15 + t142 * t16;
t43 = t140 * t57 - t83;
t44 = t143 * t57 + t146;
t19 = t140 * t43 + t143 * t44;
t116 = -0.2e1 * t182;
t169 = t116 - t173;
t84 = t105 * t188;
t168 = t140 * (t142 * t152 + t84) - t176;
t167 = -t103 * t188 - t142 * t157;
t11 = -pkin(4) * t210 + t169;
t12 = t102 * pkin(4) - qJ(5) * t210 + t105 * t74 + qJDD(5) + t15;
t164 = -pkin(4) * t12 + qJ(5) * t11;
t163 = -pkin(4) * t216 - qJ(5) * t49;
t5 = t139 * t16 - t142 * t15;
t160 = -pkin(2) + t165;
t156 = (t103 * t139 + t105 * t142) * t123;
t155 = t140 * (-t84 + t177) + t143 * t102;
t150 = t140 * (t139 * t157 - t177) + t176;
t149 = 0.2e1 * qJD(5) * t105 - t233;
t148 = -pkin(4) * t62 + qJ(5) * t214 - t12;
t132 = t143 ^ 2;
t131 = t140 ^ 2;
t130 = t132 * t145;
t128 = t131 * t145;
t120 = -t130 - t209;
t119 = -t128 - t209;
t112 = t128 + t130;
t111 = (t131 + t132) * qJDD(2);
t110 = -0.2e1 * t126 + t179;
t107 = 0.2e1 * t175 + t180;
t82 = -t140 * t119 - t185;
t81 = t143 * t120 - t186;
t55 = (qJD(4) - t123) * t103 + t159;
t50 = (-qJD(4) - t123) * t105 + t170;
t46 = -t105 * t187 + t139 * t152;
t27 = t142 * t50 + t204;
t26 = -t142 * t49 + t204;
t25 = t139 * t50 - t198;
t24 = -t139 * t49 - t198;
t22 = -t140 * t55 + t253;
t20 = -t140 * t215 - t253;
t18 = t143 * t27 - t231;
t17 = t143 * t26 - t231;
t13 = (-pkin(4) * t123 - 0.2e1 * qJD(5)) * t105 + t233;
t10 = (-t220 + t91) * pkin(4) + t149;
t9 = pkin(4) * t91 + t149 + t242;
t8 = qJ(5) * t217 + t12;
t7 = (-t210 + t217) * pkin(4) + t169;
t4 = t140 * t29 + t143 * t6;
t3 = t142 * t11 + t139 * t12;
t2 = t139 * t11 - t142 * t12;
t1 = t140 * t13 + t143 * t3;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t184, 0, 0, 0, 0, 0, 0, (qJDD(2) * t144 - t141 * t145) * t135, (-qJDD(2) * t141 - t144 * t145) * t135, 0, t136 ^ 2 * t184 + (t141 * t67 - t144 * t166 - t151) * t135, 0, 0, 0, 0, 0, 0, t136 * (t143 * t114 + t140 * t120) + (t144 * t110 + t141 * t81) * t135, t136 * (-t140 * t115 + t143 * t119) + (-t144 * t107 + t141 * t82) * t135, (t111 * t141 + t112 * t144) * t135, t136 * (t140 * t44 - t143 * t43) + (t141 * t19 - t144 * t56) * t135, 0, 0, 0, 0, 0, 0, t248, t136 * (t143 * t55 + t254) + (t141 * t22 + t252) * t135, t136 * (t140 * t27 + t227) + (t141 * t18 - t144 * t25) * t135, t136 * (t140 * t6 - t143 * t29) + (t141 * t4 - t144 * t5) * t135, 0, 0, 0, 0, 0, 0, t248, t136 * (t140 * t26 + t227) + (t141 * t17 - t144 * t24) * t135, t136 * (t143 * t215 - t254) + (t141 * t20 - t252) * t135, t136 * (-t143 * t13 + t140 * t3) + (t141 * t1 - t144 * t2) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t166, -t67, 0, 0, t162 * t140, t143 * t107 + t140 * t110, t186 + t143 * (-t128 + t209), -t161 * t143, t140 * (t130 - t209) + t185, 0, pkin(2) * t110 + pkin(7) * t81 - t143 * t56, -pkin(2) * t107 + pkin(7) * t82 + t140 * t56, pkin(2) * t112 + pkin(7) * t111 + t19, -pkin(2) * t56 + pkin(7) * t19, t168, -t249, t235, t150, t250, t155, t140 * (t207 - t244) + t143 * (t15 - t246) + t251, t140 * (t200 + t256) + t143 * (t16 + t257) + t258 + pkin(7) * t22, pkin(7) * t18 - t140 * t5 + t160 * t25, pkin(7) * t4 + t160 * t5, t168, t235, t249, t155, -t250, t150, t140 * (-t139 * t10 - t191 * t220 - t244) + t143 * (-t148 - t246) + t251, t140 * (-pkin(8) * t24 - t139 * t7 + t142 * t8) + t143 * (-pkin(3) * t24 - t163) - pkin(2) * t24 + pkin(7) * t17, t140 * (-pkin(4) * t205 + t142 * t9 - t256) + t143 * (0.2e1 * t182 - t212 - t257) - t258 + pkin(7) * t20, t140 * (-pkin(8) * t2 + (pkin(4) * t139 - t191) * t13) + t143 * (-pkin(3) * t2 - t164) - pkin(2) * t2 + pkin(7) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, t128 - t130, t180, t122, t179, qJDD(3), -t43, -t44, 0, 0, t46, t221, t243, t167, t239, t156, -pkin(3) * t220 - t200 + t245, pkin(3) * t55 + t207 + t255, pkin(8) * t27 + t234 + t6, -pkin(3) * t29 + pkin(8) * t6, t46, t243, -t221, t156, -t239, t167, t142 * t10 + t174 * t220 + t245, pkin(8) * t26 + t139 * t8 + t142 * t7 + t234, -t255 + t139 * t9 + (pkin(3) + t208) * t215, pkin(8) * t3 + (t174 - t208) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t218, t216, -t189, -t49, -t102, -t15, -t16, 0, 0, t189, t216, -t218, -t102, t49, -t189, t148, t163, t116 + t212, t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t216, t219, t12;];
tauJ_reg = t14;
