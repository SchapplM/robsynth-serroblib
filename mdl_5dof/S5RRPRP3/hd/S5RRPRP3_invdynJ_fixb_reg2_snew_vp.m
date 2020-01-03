% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:14
% EndTime: 2019-12-31 19:51:20
% DurationCPUTime: 2.08s
% Computational Cost: add. (6239->206), mult. (8880->268), div. (0->0), fcn. (5749->8), ass. (0->139)
t167 = sin(pkin(8));
t168 = cos(pkin(8));
t171 = cos(qJ(4));
t169 = sin(qJ(4));
t207 = qJD(1) + qJD(2);
t194 = t171 * t207;
t195 = t169 * t207;
t137 = t167 * t195 - t168 * t194;
t139 = t167 * t194 + t168 * t195;
t220 = t139 * t137;
t232 = qJDD(4) + t220;
t216 = t169 * t232;
t136 = t139 ^ 2;
t173 = qJD(4) ^ 2;
t236 = -t136 - t173;
t61 = -t171 * t236 + t216;
t211 = t171 * t232;
t63 = t169 * t236 + t211;
t49 = t167 * t61 - t168 * t63;
t266 = qJ(3) * t49;
t170 = sin(qJ(2));
t265 = t170 * t49;
t264 = pkin(7) * t61;
t263 = pkin(7) * t63;
t231 = t137 ^ 2;
t122 = t231 - t173;
t45 = t167 * (-t171 * t122 + t216) - t168 * (t169 * t122 + t211);
t131 = qJD(4) * t137;
t163 = qJDD(1) + qJDD(2);
t218 = t167 * t171;
t133 = (t168 * t169 + t218) * t163;
t112 = t133 - t131;
t241 = t131 - t112;
t262 = t241 * qJ(5);
t217 = t168 * t163;
t219 = t167 * t163;
t132 = t169 * t219 - t171 * t217;
t208 = t139 * qJD(4);
t109 = t132 + 0.2e1 * t208;
t233 = qJDD(4) - t220;
t215 = t169 * t233;
t234 = -t231 - t173;
t237 = t171 * t234 - t215;
t91 = t171 * t233;
t240 = t169 * t234 + t91;
t247 = -t167 * t240 + t168 * t237;
t259 = -pkin(2) * t109 + qJ(3) * t247;
t238 = -t171 * t132 + t169 * t133;
t239 = -t169 * t132 - t171 * t133;
t246 = -t167 * t239 + t168 * t238;
t84 = -t231 - t136;
t258 = -pkin(2) * t84 + qJ(3) * t246;
t172 = cos(qJ(2));
t257 = pkin(1) * (-t172 * t109 + t170 * t247);
t256 = pkin(1) * (t170 * t246 - t172 * t84);
t254 = pkin(7) * t237;
t253 = pkin(7) * t239;
t252 = pkin(7) * t240;
t248 = -pkin(3) * t84 + pkin(7) * t238;
t123 = -t136 + t173;
t245 = t167 * (-t169 * t123 + t91) + t168 * (t171 * t123 + t215);
t228 = sin(qJ(1));
t229 = cos(qJ(1));
t185 = t229 * g(1) + t228 * g(2);
t149 = -qJD(1) ^ 2 * pkin(1) - t185;
t184 = t228 * g(1) - t229 * g(2);
t183 = qJDD(1) * pkin(1) + t184;
t118 = t172 * t149 + t170 * t183;
t206 = t207 ^ 2;
t187 = -t206 * pkin(2) + t163 * qJ(3) + 0.2e1 * qJD(3) * t207 + t118;
t174 = t167 ^ 2;
t176 = t168 ^ 2;
t192 = t176 * t206;
t148 = t174 * t206 + t192;
t191 = -t167 * g(3) + t187 * t168;
t226 = t168 * g(3);
t55 = t168 * t191 + t167 * (t187 * t167 + t226);
t235 = t136 - t231;
t230 = 2 * qJD(5);
t227 = pkin(4) * t171;
t178 = -t226 + (t206 * t168 * pkin(3) - t163 * pkin(7) - t187) * t167;
t69 = -pkin(3) * t192 + pkin(7) * t217 + t191;
t41 = t169 * t178 + t171 * t69;
t117 = -t170 * t149 + t172 * t183;
t101 = -t163 * pkin(2) - t206 * qJ(3) + qJDD(3) - t117;
t225 = -pkin(2) * t101 + qJ(3) * t55;
t40 = t169 * t69 - t171 * t178;
t18 = t169 * t41 - t171 * t40;
t224 = t167 * t18;
t75 = -pkin(3) * t217 - t148 * pkin(7) + t101;
t223 = t169 * t75;
t222 = t169 * t241;
t221 = t171 * t75;
t214 = t169 * t109;
t210 = t171 * t109;
t209 = t172 * t163;
t143 = t148 * t168;
t205 = pkin(2) * t217 - qJ(3) * t143 - t168 * t101;
t100 = t137 * pkin(4) - t139 * qJ(5);
t189 = qJDD(4) * qJ(5) + qJD(4) * t230 - t137 * t100 + t41;
t25 = -t173 * pkin(4) + t189;
t26 = -qJDD(4) * pkin(4) - t173 * qJ(5) + t139 * t100 + qJDD(5) + t40;
t11 = t169 * t25 - t171 * t26;
t12 = t169 * t26 + t171 * t25;
t201 = -qJ(5) * t169 - pkin(3);
t110 = -t132 - t208;
t182 = -t110 * pkin(4) + t262 + t75;
t28 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t139 + t182;
t3 = -t167 * t11 + t168 * t12;
t204 = t167 * (-pkin(7) * t11 + (pkin(4) * t169 - qJ(5) * t171) * t28) + t168 * (pkin(7) * t12 + (t201 - t227) * t28) - pkin(2) * t28 + qJ(3) * t3;
t20 = (-t173 - t84) * pkin(4) + t189;
t21 = -qJ(5) * t84 + t26;
t203 = t167 * (-t169 * t20 + t171 * t21 - t253) + t168 * (t169 * t21 + t171 * t20 + t248) + t258;
t181 = t139 * t230 - t182;
t22 = -pkin(4) * t208 + t181 - t262;
t202 = t167 * (pkin(4) * t222 + t171 * t22 - t264) + t168 * (t263 + t169 * t22 - (pkin(3) + t227) * t241) - pkin(2) * t241 - t266;
t19 = t169 * t40 + t171 * t41;
t200 = t167 * (-t18 - t253) + t168 * (t19 + t248) + t258;
t23 = (-t109 - t208) * pkin(4) + t181;
t199 = t167 * (-qJ(5) * t210 - t169 * t23 - t252) + t168 * (t201 * t109 + t171 * t23 + t254) + t259;
t198 = t167 * (t223 - t252) + t168 * (-pkin(3) * t109 - t221 + t254) + t259;
t111 = t133 - 0.2e1 * t131;
t197 = t167 * (t221 + t264) + t168 * (-pkin(3) * t111 + t223 - t263) - pkin(2) * t111 + t266;
t158 = t174 * t163;
t159 = t176 * t163;
t146 = t159 + t158;
t196 = pkin(2) * t148 + qJ(3) * t146 + t55;
t142 = t148 * t167;
t188 = -pkin(2) * t219 + qJ(3) * t142 + t167 * t101;
t7 = t168 * t19 - t224;
t186 = -pkin(7) * t224 + qJ(3) * t7 + t168 * (-pkin(3) * t75 + pkin(7) * t19) - pkin(2) * t75;
t180 = t167 * (-t169 * t110 + t171 * t131) + t168 * (t171 * t110 + t169 * t131);
t120 = t169 * t208;
t179 = t167 * t120 + (-t137 * t218 + t168 * (-t137 * t169 - t139 * t171)) * qJD(4);
t150 = 0.2e1 * t167 * t217;
t52 = t167 * (t171 * t112 - t120) + t168 * (t169 * t112 + t171 * t208);
t36 = t167 * (-t169 * t111 - t210) + t168 * (t171 * t111 - t214);
t24 = t167 * (t210 - t222) + t168 * (t171 * t241 + t214);
t1 = [0, 0, 0, 0, 0, qJDD(1), t184, t185, 0, 0, 0, 0, 0, 0, 0, t163, pkin(1) * (-t170 * t206 + t209) + t117, pkin(1) * (-t170 * t163 - t172 * t206) - t118, 0, pkin(1) * (t172 * t117 + t170 * t118), t158, t150, 0, t159, 0, 0, pkin(1) * (-t170 * t143 + t168 * t209) + t205, pkin(1) * (t170 * t142 - t167 * t209) + t188, pkin(1) * (t170 * t146 + t172 * t148) + t196, pkin(1) * (-t172 * t101 + t170 * t55) + t225, t52, t36, t245, t180, -t45, t179, t257 + t198, pkin(1) * (-t172 * t111 + t265) + t197, t256 + t200, pkin(1) * (t170 * t7 - t172 * t75) + t186, t52, t245, t24, t179, t45, t180, t257 + t199, t256 + t203, pkin(1) * (-t172 * t241 - t265) + t202, pkin(1) * (t170 * t3 - t172 * t28) + t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, t117, -t118, 0, 0, t158, t150, 0, t159, 0, 0, t205, t188, t196, t225, t52, t36, t245, t180, -t45, t179, t198, t197, t200, t186, t52, t245, t24, t179, t45, t180, t199, t203, t202, t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t217, t219, -t148, t101, 0, 0, 0, 0, 0, 0, t109, t111, t84, t75, 0, 0, 0, 0, 0, 0, t109, t84, t241, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220, t235, t133, -t220, -t132, qJDD(4), -t40, -t41, 0, 0, t220, t131 + t112, -t235, qJDD(4), t132, -t220, pkin(4) * t233 + qJ(5) * t234 - t26, -pkin(4) * t133 - qJ(5) * t132, qJ(5) * t232 + (-t236 - t173) * pkin(4) + t189, -pkin(4) * t26 + qJ(5) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t233, t133, t236, t26;];
tauJ_reg = t1;
