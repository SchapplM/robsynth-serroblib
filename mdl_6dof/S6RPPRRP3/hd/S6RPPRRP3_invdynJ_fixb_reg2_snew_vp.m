% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:51:55
% EndTime: 2019-05-05 14:52:02
% DurationCPUTime: 2.09s
% Computational Cost: add. (4792->236), mult. (8973->284), div. (0->0), fcn. (5357->8), ass. (0->171)
t141 = sin(qJ(5));
t144 = cos(qJ(5));
t145 = cos(qJ(4));
t187 = qJD(1) * qJD(4);
t181 = t145 * t187;
t142 = sin(qJ(4));
t186 = t142 * qJDD(1);
t112 = -t181 - t186;
t105 = qJDD(5) - t112;
t189 = qJD(1) * t145;
t106 = -t144 * qJD(4) + t141 * t189;
t108 = t141 * qJD(4) + t144 * t189;
t198 = t108 * t106;
t222 = t105 + t198;
t203 = t144 * t222;
t104 = t108 ^ 2;
t126 = t142 * qJD(1) + qJD(5);
t217 = t126 ^ 2;
t226 = -t104 - t217;
t33 = t141 * t226 + t203;
t259 = pkin(8) * t33;
t258 = t142 * t33;
t257 = t145 * t33;
t208 = t141 * t222;
t128 = t145 * qJDD(1);
t182 = t142 * t187;
t113 = t128 - t182;
t176 = t144 * qJDD(4) - t141 * t113;
t158 = t108 * qJD(5) - t176;
t94 = t126 * t108;
t47 = t158 - t94;
t218 = t106 ^ 2;
t87 = t218 - t217;
t256 = t142 * t47 + t145 * (-t144 * t87 + t208);
t169 = t142 * pkin(4) - t145 * pkin(8);
t138 = sin(pkin(9));
t180 = pkin(1) * t138 + qJ(3);
t155 = t169 + t180;
t255 = t155 * (-t144 * t226 + t208);
t254 = t141 * t87 + t203;
t161 = -t141 * qJDD(4) - t144 * t113;
t154 = -t106 * qJD(5) - t161;
t199 = t106 * t126;
t221 = t154 - t199;
t210 = t141 * t221;
t225 = t104 - t218;
t227 = t158 + t94;
t252 = -t142 * t225 + t145 * (t144 * t227 + t210);
t139 = cos(pkin(9));
t215 = -pkin(2) - pkin(7);
t174 = pkin(1) * t139 - t215;
t219 = -t217 - t218;
t223 = t105 - t198;
t207 = t141 * t223;
t230 = t144 * t219 - t207;
t245 = t142 * t230 - t145 * t227;
t59 = t144 * t223;
t251 = -t174 * t245 + t155 * (t141 * t219 + t59);
t250 = pkin(8) * t230;
t244 = t142 * t227 + t145 * t230;
t220 = t154 + t199;
t88 = -t104 + t217;
t243 = t142 * t220 + t145 * (-t141 * t88 + t59);
t224 = t104 + t218;
t242 = pkin(4) * t224;
t240 = t144 * t88 + t207;
t239 = qJ(6) * t221;
t236 = t142 * t224;
t232 = t145 * t224;
t229 = -t141 * t227 + t144 * t221;
t216 = qJD(4) ^ 2;
t214 = pkin(5) * t144;
t213 = t217 * pkin(5);
t163 = -t113 + t182;
t164 = -t112 + t181;
t147 = qJD(1) ^ 2;
t130 = 0.2e1 * qJD(3) * qJD(1);
t131 = qJDD(1) * qJ(3);
t143 = sin(qJ(1));
t146 = cos(qJ(1));
t166 = t146 * g(1) + t143 * g(2);
t109 = -t147 * pkin(1) - t166;
t165 = t143 * g(1) - t146 * g(2);
t156 = qJDD(1) * pkin(1) + t165;
t212 = t139 * t109 + t138 * t156;
t175 = t130 + t131 + t212;
t57 = t215 * t147 + t175;
t30 = pkin(4) * t164 + pkin(8) * t163 + t57;
t110 = t169 * qJD(1);
t137 = qJDD(1) * pkin(2);
t177 = -t138 * t109 + t139 * t156;
t61 = -t147 * qJ(3) + qJDD(3) - t137 - t177;
t152 = -qJDD(1) * pkin(7) + t61;
t135 = -g(3) + qJDD(2);
t191 = t145 * t135;
t37 = -t216 * pkin(4) + qJDD(4) * pkin(8) + t191 + (-qJD(1) * t110 + t152) * t142;
t16 = t141 * t30 + t144 * t37;
t209 = t141 * t220;
t54 = t142 * t135 - t145 * t152;
t36 = -qJDD(4) * pkin(4) - t216 * pkin(8) + t110 * t189 + t54;
t150 = t158 * pkin(5) - t239 + t36;
t13 = (pkin(5) * t126 - (2 * qJD(6))) * t108 + t150;
t202 = t145 * t13;
t201 = t145 * t36;
t200 = qJ(6) * t144;
t197 = t126 * t141;
t196 = t126 * t144;
t133 = t142 ^ 2;
t195 = t133 * t147;
t134 = t145 ^ 2;
t194 = t134 * t147;
t184 = t142 * t147 * t145;
t119 = qJDD(4) + t184;
t193 = t142 * t119;
t120 = qJDD(4) - t184;
t192 = t145 * t120;
t190 = t133 + t134;
t188 = qJD(6) * t126;
t185 = t142 * t198;
t183 = t106 * t196;
t179 = -qJ(6) * t141 - pkin(4);
t15 = t141 * t37 - t144 * t30;
t6 = t141 * t15 + t144 * t16;
t75 = t106 * pkin(5) - t108 * qJ(6);
t178 = t105 * qJ(6) - t106 * t75 + t16;
t173 = -pkin(1) * (t138 * qJDD(1) + t139 * t147) - t212;
t121 = 0.2e1 * t188;
t172 = t121 + t178;
t86 = t108 * t197;
t171 = t145 * (t144 * t154 - t86) + t185;
t170 = t106 * t197 - t144 * t158;
t11 = t172 - t213;
t12 = -t105 * pkin(5) - qJ(6) * t217 + t108 * t75 + qJDD(6) + t15;
t168 = pkin(5) * t12 - qJ(6) * t11;
t167 = pkin(5) * t220 + qJ(6) * t47;
t162 = t141 * t16 - t144 * t15;
t55 = t142 * t152 + t191;
t27 = t142 * t55 - t145 * t54;
t160 = -pkin(1) * (-t139 * qJDD(1) + t138 * t147) + t177;
t159 = qJ(6) * t222 + t178;
t157 = (-t106 * t141 - t108 * t144) * t126;
t153 = t145 * (t86 - t183) + t142 * t105;
t151 = t145 * (t141 * t158 + t183) - t185;
t149 = 0.2e1 * qJD(6) * t108 - t150;
t148 = pkin(5) * t223 + qJ(6) * t219 - t12;
t124 = -t194 - t216;
t123 = -t195 - t216;
t117 = t190 * qJDD(1);
t114 = t128 - 0.2e1 * t182;
t111 = 0.2e1 * t181 + t186;
t85 = t145 * t124 - t193;
t84 = t142 * t123 + t192;
t60 = -t147 * pkin(2) + t175;
t53 = (qJD(5) + t126) * t106 + t161;
t48 = (-qJD(5) + t126) * t108 + t176;
t44 = t144 * t220;
t43 = t108 * t196 + t141 * t154;
t26 = t144 * t48 + t209;
t25 = -t144 * t47 + t209;
t21 = t145 * t53 - t258;
t19 = t145 * t221 + t258;
t18 = t142 * t26 + t232;
t17 = t142 * t25 + t232;
t10 = (-t227 - t94) * pkin(5) + t149;
t9 = -pkin(5) * t94 + t149 + t239;
t8 = qJ(6) * t224 + t12;
t7 = (-t217 + t224) * pkin(5) + t172;
t4 = t142 * t6 - t201;
t3 = t144 * t11 + t141 * t12;
t1 = t142 * t3 - t202;
t2 = [0, 0, 0, 0, 0, qJDD(1), t165, t166, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t160, t173, 0, pkin(1) * (t138 * t212 + t139 * t177), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - 0.2e1 * t137 - t160, t130 + 0.2e1 * t131 - t173, pkin(1) * (t138 * t60 - t139 * t61) - pkin(2) * t61 + qJ(3) * t60, -t163 * t145, -t145 * t111 - t142 * t114, t192 - t142 * (-t194 + t216), t164 * t142, t145 * (t195 - t216) - t193, 0, t180 * t111 + t142 * t57 - t174 * t84, t180 * t114 + t145 * t57 - t174 * t85, -t180 * t190 * t147 + t174 * t117 - t27, -t174 * t27 + t180 * t57, t171, -t252, t243, t151, -t256, t153, t141 * t201 - t142 * t15 + t251, -t142 * t16 + t144 * t201 - t174 * t21 - t255, -t145 * t162 + t155 * (t141 * t48 - t44) - t174 * t18, t155 * t162 - t174 * t4, t171, t243, t252, t153, t256, t151, t145 * (-t141 * t10 - t200 * t227) + t142 * t148 + t251, t145 * (-t141 * t7 + t144 * t8) - t142 * t167 + t155 * (-t141 * t47 - t44) - t174 * t17, t145 * (-pkin(5) * t210 + t144 * t9) - t142 * (pkin(5) * t226 - t159 - 0.2e1 * t188 + t213) + t255 - t174 * t19, -t142 * t168 + (pkin(5) * t141 - t200) * t202 + t155 * (t141 * t11 - t144 * t12) - t174 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, 0, 0, 0, 0, 0, 0, -t142 * t120 + t145 * t123, -t145 * t119 - t142 * t124, 0, t142 * t54 + t145 * t55, 0, 0, 0, 0, 0, 0, t244, -t142 * t53 - t257, t145 * t26 - t236, t142 * t36 + t145 * t6, 0, 0, 0, 0, 0, 0, t244, t145 * t25 - t236, -t142 * t221 + t257, t142 * t13 + t145 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t147, t61, 0, 0, 0, 0, 0, 0, t84, t85, -t117, t27, 0, 0, 0, 0, 0, 0, t245, t21, t18, t4, 0, 0, 0, 0, 0, 0, t245, t17, t19, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, (-t133 + t134) * t147, t128, -t184, -t186, qJDD(4), -t54, -t55, 0, 0, t43, t229, t240, t170, t254, t157, -pkin(4) * t227 - t144 * t36 + t250, pkin(4) * t53 + t141 * t36 - t259, pkin(8) * t26 + t242 + t6, -pkin(4) * t36 + pkin(8) * t6, t43, t240, -t229, t157, -t254, t170, t144 * t10 + t179 * t227 + t250, pkin(8) * t25 + t141 * t8 + t144 * t7 + t242, t259 + t141 * t9 + (pkin(4) + t214) * t221, pkin(8) * t3 + (t179 - t214) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t225, t220, -t198, -t47, t105, -t15, -t16, 0, 0, t198, t220, -t225, t105, t47, -t198, t148, -t167, t121 + (-t217 - t226) * pkin(5) + t159, -t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t223, t220, t226, t12;];
tauJ_reg  = t2;
