% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:38:10
% EndTime: 2019-03-09 02:38:21
% DurationCPUTime: 6.88s
% Computational Cost: add. (6127->489), mult. (15300->687), div. (0->0), fcn. (10788->10), ass. (0->230)
t180 = sin(pkin(10));
t187 = cos(qJ(3));
t237 = cos(pkin(10));
t205 = t237 * t187;
t185 = sin(qJ(3));
t224 = qJD(1) * t185;
t148 = -qJD(1) * t205 + t180 * t224;
t271 = -t148 / 0.2e1;
t162 = t180 * t187 + t185 * t237;
t150 = t162 * qJD(1);
t179 = sin(pkin(11));
t182 = cos(pkin(11));
t126 = qJD(3) * t179 + t150 * t182;
t211 = -cos(pkin(9)) * pkin(1) - pkin(2);
t165 = -pkin(3) * t187 + t211;
t225 = qJD(1) * t165;
t147 = qJD(4) + t225;
t254 = Ifges(6,2) * t179;
t256 = Ifges(6,4) * t182;
t196 = -t254 + t256;
t257 = Ifges(6,4) * t179;
t197 = Ifges(6,1) * t182 - t257;
t198 = mrSges(6,1) * t179 + mrSges(6,2) * t182;
t203 = t182 * qJD(3) - t150 * t179;
t235 = Ifges(5,5) * qJD(3);
t267 = t182 / 0.2e1;
t268 = -t179 / 0.2e1;
t173 = sin(pkin(9)) * pkin(1) + pkin(7);
t166 = t173 * qJD(1);
t201 = qJ(4) * qJD(1) + t166;
t220 = t185 * qJD(2);
t129 = t187 * t201 + t220;
t119 = t180 * t129;
t178 = t187 * qJD(2);
t128 = -t185 * t201 + t178;
t122 = qJD(3) * pkin(3) + t128;
t72 = t122 * t237 - t119;
t66 = -qJD(3) * pkin(4) + qJD(5) - t72;
t301 = (t126 * Ifges(6,1) + Ifges(6,4) * t203 + t148 * Ifges(6,5)) * t267 + (t126 * Ifges(6,4) + Ifges(6,2) * t203 + t148 * Ifges(6,6)) * t268 + t147 * mrSges(5,2) + t66 * t198 + t126 * t197 / 0.2e1 + t203 * t196 / 0.2e1 + t235 / 0.2e1;
t184 = sin(qJ(6));
t186 = cos(qJ(6));
t300 = -t126 * t184 + t186 * t203;
t76 = t126 * t186 + t184 * t203;
t246 = t150 * Ifges(5,4);
t251 = t126 * Ifges(6,5);
t252 = t203 * Ifges(6,6);
t206 = t237 * t129;
t73 = t180 * t122 + t206;
t67 = qJD(3) * qJ(5) + t73;
t90 = t148 * pkin(4) - t150 * qJ(5) + t147;
t29 = -t179 * t67 + t182 * t90;
t30 = t179 * t90 + t182 * t67;
t16 = pkin(5) * t148 - pkin(8) * t126 + t29;
t20 = pkin(8) * t203 + t30;
t5 = t16 * t186 - t184 * t20;
t6 = t16 * t184 + t186 * t20;
t299 = t30 * mrSges(6,2) + t6 * mrSges(7,2) + Ifges(5,2) * t271 + Ifges(5,6) * qJD(3) + t246 / 0.2e1 - t147 * mrSges(5,1) - t29 * mrSges(6,1) - t5 * mrSges(7,1) - t251 / 0.2e1 - t252 / 0.2e1;
t189 = -t180 * t185 + t205;
t151 = t189 * qJD(3);
t141 = qJD(1) * t151;
t191 = t179 * t184 - t182 * t186;
t38 = qJD(6) * t300 - t141 * t191;
t280 = t38 / 0.2e1;
t163 = t179 * t186 + t182 * t184;
t39 = -qJD(6) * t76 - t141 * t163;
t279 = t39 / 0.2e1;
t149 = t162 * qJD(3);
t140 = qJD(1) * t149;
t274 = t140 / 0.2e1;
t265 = pkin(3) * t180;
t171 = qJ(5) + t265;
t261 = pkin(8) + t171;
t156 = t261 * t179;
t157 = t261 * t182;
t113 = -t156 * t184 + t157 * t186;
t229 = t148 * t182;
t216 = pkin(3) * t224;
t102 = pkin(4) * t150 + qJ(5) * t148 + t216;
t80 = t128 * t237 - t119;
t40 = t182 * t102 - t179 * t80;
t21 = pkin(5) * t150 + pkin(8) * t229 + t40;
t230 = t148 * t179;
t41 = t179 * t102 + t182 * t80;
t28 = pkin(8) * t230 + t41;
t297 = -qJD(5) * t163 - qJD(6) * t113 + t184 * t28 - t186 * t21;
t112 = -t156 * t186 - t157 * t184;
t296 = -qJD(5) * t191 + qJD(6) * t112 - t184 * t21 - t186 * t28;
t13 = -t39 * mrSges(7,1) + t38 * mrSges(7,2);
t231 = t141 * t182;
t232 = t141 * t179;
t91 = mrSges(6,1) * t232 + mrSges(6,2) * t231;
t295 = t13 + t91;
t248 = t150 * mrSges(5,3);
t294 = qJD(3) * mrSges(5,1) + mrSges(6,1) * t203 - mrSges(6,2) * t126 - t248;
t152 = t191 * qJD(6);
t97 = t191 * t148;
t293 = t152 + t97;
t153 = t163 * qJD(6);
t96 = t163 * t148;
t292 = t153 + t96;
t218 = qJD(1) * qJD(3);
t208 = t185 * t218;
t98 = -mrSges(6,2) * t140 - mrSges(6,3) * t232;
t99 = mrSges(6,1) * t140 - mrSges(6,3) * t231;
t289 = -t179 * t99 + t182 * t98;
t92 = -mrSges(6,2) * t148 + mrSges(6,3) * t203;
t93 = mrSges(6,1) * t148 - mrSges(6,3) * t126;
t288 = -t179 * t93 + t182 * t92;
t222 = qJD(3) * t185;
t136 = qJD(3) * t178 - t166 * t222;
t221 = qJD(4) * t187;
t117 = (-qJ(4) * t222 + t221) * qJD(1) + t136;
t219 = t185 * qJD(4);
t188 = -qJD(1) * t219 - qJD(3) * t129;
t55 = t237 * t117 + t180 * t188;
t50 = qJD(3) * qJD(5) + t55;
t200 = pkin(3) * t208;
t64 = pkin(4) * t140 - qJ(5) * t141 - qJD(5) * t150 + t200;
t18 = -t179 * t50 + t182 * t64;
t19 = t179 * t64 + t182 * t50;
t193 = -t179 * t18 + t182 * t19;
t14 = pkin(5) * t140 - pkin(8) * t231 + t18;
t15 = -pkin(8) * t232 + t19;
t1 = qJD(6) * t5 + t14 * t184 + t15 * t186;
t2 = -qJD(6) * t6 + t14 * t186 - t15 * t184;
t287 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t38 + Ifges(7,6) * t39;
t284 = Ifges(7,4) * t280 + Ifges(7,2) * t279 + Ifges(7,6) * t274;
t283 = Ifges(7,1) * t280 + Ifges(7,4) * t279 + Ifges(7,5) * t274;
t146 = qJD(6) + t148;
t266 = Ifges(7,4) * t76;
t23 = Ifges(7,2) * t300 + Ifges(7,6) * t146 + t266;
t282 = t23 / 0.2e1;
t71 = Ifges(7,4) * t300;
t24 = Ifges(7,1) * t76 + Ifges(7,5) * t146 + t71;
t281 = t24 / 0.2e1;
t278 = -t300 / 0.2e1;
t277 = t300 / 0.2e1;
t276 = -t76 / 0.2e1;
t275 = t76 / 0.2e1;
t273 = -t146 / 0.2e1;
t272 = t146 / 0.2e1;
t270 = t148 / 0.2e1;
t269 = -t150 / 0.2e1;
t264 = pkin(8) * t182;
t263 = t300 * Ifges(7,6);
t262 = t76 * Ifges(7,5);
t78 = pkin(3) * t222 + pkin(4) * t149 - qJ(5) * t151 - qJD(5) * t162;
t226 = qJ(4) + t173;
t202 = qJD(3) * t226;
t132 = -t185 * t202 + t221;
t133 = -t187 * t202 - t219;
t89 = t132 * t237 + t180 * t133;
t34 = t179 * t78 + t182 * t89;
t260 = mrSges(5,3) * t140;
t259 = mrSges(5,3) * t141;
t258 = Ifges(4,4) * t185;
t158 = t226 * t185;
t159 = t226 * t187;
t114 = t237 * t158 + t159 * t180;
t54 = t117 * t180 - t237 * t188;
t253 = t114 * t54;
t250 = t146 * Ifges(7,3);
t249 = t148 * mrSges(5,3);
t247 = t150 * Ifges(5,1);
t245 = t189 * t54;
t236 = Ifges(4,5) * qJD(3);
t234 = Ifges(4,6) * qJD(3);
t228 = t151 * t179;
t227 = t162 * t179;
t104 = -pkin(4) * t189 - t162 * qJ(5) + t165;
t115 = -t180 * t158 + t159 * t237;
t48 = t179 * t104 + t182 * t115;
t223 = qJD(1) * t187;
t31 = -mrSges(7,1) * t300 + mrSges(7,2) * t76;
t217 = t31 - t294;
t215 = mrSges(4,3) * t224;
t214 = mrSges(4,3) * t223;
t176 = Ifges(4,4) * t223;
t210 = t237 * pkin(3);
t209 = m(4) * t173 + mrSges(4,3);
t33 = -t179 * t89 + t182 * t78;
t204 = t140 * mrSges(5,1) + t141 * mrSges(5,2);
t47 = t182 * t104 - t115 * t179;
t79 = t128 * t180 + t206;
t88 = t132 * t180 - t237 * t133;
t174 = -t210 - pkin(4);
t168 = t211 * qJD(1);
t195 = Ifges(6,5) * t182 - Ifges(6,6) * t179;
t194 = t179 * t19 + t18 * t182;
t192 = t179 * t29 - t182 * t30;
t32 = -pkin(5) * t189 - t162 * t264 + t47;
t42 = -pkin(8) * t227 + t48;
t11 = -t184 * t42 + t186 * t32;
t12 = t184 * t32 + t186 * t42;
t143 = t166 * t187 + t220;
t134 = -qJD(3) * mrSges(5,2) - t249;
t190 = t134 + t288;
t169 = -qJD(3) * mrSges(4,2) + t214;
t167 = qJD(3) * mrSges(4,1) - t215;
t164 = -t182 * pkin(5) + t174;
t155 = Ifges(4,1) * t224 + t176 + t236;
t154 = t234 + (Ifges(4,2) * t187 + t258) * qJD(1);
t145 = Ifges(5,4) * t148;
t142 = -t166 * t185 + t178;
t139 = Ifges(7,3) * t140;
t137 = t143 * qJD(3);
t116 = mrSges(5,1) * t148 + mrSges(5,2) * t150;
t108 = -t145 + t235 + t247;
t106 = t191 * t162;
t105 = t163 * t162;
t86 = pkin(5) * t227 + t114;
t63 = t140 * Ifges(6,5) + t141 * t197;
t62 = t140 * Ifges(6,6) + t141 * t196;
t57 = t148 * Ifges(6,3) + t251 + t252;
t56 = pkin(5) * t228 + t88;
t53 = -pkin(5) * t230 + t79;
t52 = mrSges(7,1) * t146 - mrSges(7,3) * t76;
t51 = -mrSges(7,2) * t146 + mrSges(7,3) * t300;
t46 = -t151 * t163 + t152 * t162;
t45 = -t151 * t191 - t153 * t162;
t44 = -pkin(5) * t203 + t66;
t43 = pkin(5) * t232 + t54;
t27 = -mrSges(7,2) * t140 + mrSges(7,3) * t39;
t26 = mrSges(7,1) * t140 - mrSges(7,3) * t38;
t25 = -pkin(8) * t228 + t34;
t22 = t250 + t262 + t263;
t17 = pkin(5) * t149 - t151 * t264 + t33;
t4 = -qJD(6) * t12 + t17 * t186 - t184 * t25;
t3 = qJD(6) * t11 + t17 * t184 + t186 * t25;
t7 = [(Ifges(7,5) * t45 + Ifges(7,6) * t46) * t272 + (Ifges(7,1) * t45 + Ifges(7,4) * t46) * t275 + (Ifges(7,4) * t45 + Ifges(7,2) * t46) * t277 + t45 * t281 + t46 * t282 - t106 * t283 - t105 * t284 + (t209 * t137 + (-t173 * t169 - t154 / 0.2e1 + t168 * mrSges(4,1) - t234 / 0.2e1 - t209 * t143 + (t211 * mrSges(4,1) - 0.3e1 / 0.2e1 * t258 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t187) * qJD(1) + (t116 + m(5) * (t147 + t225) + qJD(1) * (-mrSges(5,1) * t189 + mrSges(5,2) * t162)) * pkin(3)) * qJD(3)) * t185 + (t114 * t141 - t115 * t140 - t149 * t73 - t151 * t72 + t162 * t54 + t189 * t55) * mrSges(5,3) + (-t140 * t162 + t141 * t189 + t149 * t269 + t151 * t271) * Ifges(5,4) - (t18 * mrSges(6,1) - t19 * mrSges(6,2) + t139 / 0.2e1 + t195 * t141 + (Ifges(5,2) + Ifges(6,3) + Ifges(7,3) / 0.2e1) * t140 + t287) * t189 + (-t1 * t105 + t106 * t2 - t45 * t5 + t46 * t6) * mrSges(7,3) + (-Ifges(7,5) * t106 - Ifges(7,6) * t105) * t274 + (-Ifges(7,4) * t106 - Ifges(7,2) * t105) * t279 + (-Ifges(7,1) * t106 - Ifges(7,4) * t105) * t280 + t43 * (mrSges(7,1) * t105 - mrSges(7,2) * t106) - t294 * t88 + m(7) * (t1 * t12 + t11 * t2 + t3 * t6 + t4 * t5 + t43 * t86 + t44 * t56) + (t57 / 0.2e1 + t263 / 0.2e1 + t262 / 0.2e1 + t250 / 0.2e1 + t22 / 0.2e1 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t148 - t299) * t149 + t89 * t134 + t114 * t91 + t48 * t98 + t47 * t99 + t165 * t204 + (t247 / 0.2e1 + t195 * t270 + t108 / 0.2e1 + (-t179 * t30 - t182 * t29) * mrSges(6,3) + t301) * t151 + (t209 * t136 + (-t173 * t167 + t155 / 0.2e1 + 0.3e1 / 0.2e1 * t176 + t236 / 0.2e1 - t209 * t142 + 0.2e1 * t168 * mrSges(4,2)) * qJD(3)) * t187 + m(5) * (t115 * t55 - t72 * t88 + t73 * t89 + t253) + m(6) * (t18 * t47 + t19 * t48 + t29 * t33 + t30 * t34 + t66 * t88 + t253) + (t54 * t198 + t195 * t274 + t62 * t268 + t63 * t267 - t194 * mrSges(6,3) + (Ifges(5,1) + Ifges(6,1) * t182 ^ 2 / 0.2e1 + (-t256 + t254 / 0.2e1) * t179) * t141) * t162 + t11 * t26 + t12 * t27 + t44 * (-t46 * mrSges(7,1) + t45 * mrSges(7,2)) + t3 * t51 + t4 * t52 + t56 * t31 + t86 * t13 + t34 * t92 + t33 * t93; -t105 * t26 - t106 * t27 + t45 * t51 + t46 * t52 + (-t260 + t289) * t162 - (t259 + t295) * t189 + t190 * t151 + t217 * t149 + (-t185 * t167 + t187 * t169 + (-t185 ^ 2 - t187 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + m(7) * (-t1 * t106 - t105 * t2 + t149 * t44 - t189 * t43 + t45 * t6 + t46 * t5) + m(6) * (t149 * t66 - t151 * t192 + t162 * t193 - t245) + m(4) * (t136 * t185 - t137 * t187 + (-t142 * t185 + t143 * t187) * qJD(3)) + m(5) * (-t149 * t72 + t151 * t73 + t162 * t55 - t245); (Ifges(7,4) * t163 - Ifges(7,2) * t191) * t279 + (Ifges(7,1) * t163 - Ifges(7,4) * t191) * t280 + (Ifges(6,5) * t179 + Ifges(7,5) * t163 + Ifges(6,6) * t182 - Ifges(7,6) * t191) * t274 + t43 * (mrSges(7,1) * t191 + mrSges(7,2) * t163) + (-t1 * t191 - t163 * t2 - t292 * t6 + t293 * t5) * mrSges(7,3) - t191 * t284 + (mrSges(7,1) * t292 - mrSges(7,2) * t293) * t44 + t294 * t79 + t62 * t267 + t288 * qJD(5) + t289 * t171 + (-t168 * (mrSges(4,1) * t185 + mrSges(4,2) * t187) - (Ifges(4,1) * t187 - t258) * t224 / 0.2e1) * qJD(1) + (-t229 * t29 - t230 * t30 + t193) * mrSges(6,3) + (t167 + t215) * t143 + (-qJD(5) * t192 + t171 * t193 + t174 * t54 - t29 * t40 - t30 * t41 - t66 * t79) * m(6) - t152 * t281 - t153 * t282 + t163 * t283 + t73 * t248 + (Ifges(7,1) * t97 + Ifges(7,4) * t96) * t276 + (-Ifges(7,5) * t152 - Ifges(7,6) * t153) * t272 + (-Ifges(7,1) * t152 - Ifges(7,4) * t153) * t275 + (-Ifges(7,4) * t152 - Ifges(7,2) * t153) * t277 + (-t169 + t214) * t142 + (-Ifges(5,1) * t148 + t22 - t246 + t57) * t269 + (Ifges(7,5) * t97 + Ifges(7,6) * t96) * t273 - (-Ifges(4,2) * t224 + t155 + t176) * t223 / 0.2e1 + t218 * Ifges(4,5) * t187 / 0.2e1 + (Ifges(7,4) * t97 + Ifges(7,2) * t96) * t278 + t296 * t51 + ((t180 * t55 - t237 * t54) * pkin(3) - t147 * t216 + t72 * t79 - t73 * t80) * m(5) + (t108 - t145) * t270 + t179 * t63 / 0.2e1 + t174 * t91 + t164 * t13 + (Ifges(7,5) * t276 - Ifges(5,2) * t270 + Ifges(7,6) * t278 + Ifges(6,3) * t271 + Ifges(7,3) * t273 + t299) * t150 + Ifges(5,5) * t141 - Ifges(5,6) * t140 - t80 * t134 - t136 * mrSges(4,2) - t137 * mrSges(4,1) + t112 * t26 + t113 * t27 + t297 * t52 + (t1 * t113 + t112 * t2 + t164 * t43 + t296 * t6 + t297 * t5 - t44 * t53) * m(7) - Ifges(4,6) * t208 / 0.2e1 - t116 * t216 + (-t195 * t271 + t301) * t148 + t154 * t224 / 0.2e1 - t72 * t249 + (Ifges(6,1) * t179 + t256) * t231 / 0.2e1 - (Ifges(6,2) * t182 + t257) * t232 / 0.2e1 - t210 * t259 - t260 * t265 + (-mrSges(6,1) * t182 + mrSges(6,2) * t179 - mrSges(5,1)) * t54 - t53 * t31 - t55 * mrSges(5,2) - t41 * t92 - t40 * t93 - t96 * t23 / 0.2e1 - t97 * t24 / 0.2e1; -t191 * t26 + t163 * t27 + t179 * t98 + t182 * t99 - t292 * t52 - t293 * t51 - t217 * t150 + t190 * t148 + t204 + (t1 * t163 - t150 * t44 - t191 * t2 - t292 * t5 - t293 * t6) * m(7) + (-t148 * t192 - t150 * t66 + t194) * m(6) + (t148 * t73 + t150 * t72 + t200) * m(5); t126 * t93 - t203 * t92 - t300 * t51 + t76 * t52 + (-t300 * t6 + t5 * t76 + t43) * m(7) + (t126 * t29 - t203 * t30 + t54) * m(6) + t295; t139 - t44 * (mrSges(7,1) * t76 + mrSges(7,2) * t300) + (Ifges(7,1) * t300 - t266) * t276 + t23 * t275 + (Ifges(7,5) * t300 - Ifges(7,6) * t76) * t273 - t5 * t51 + t6 * t52 + (t300 * t5 + t6 * t76) * mrSges(7,3) + (-Ifges(7,2) * t76 + t24 + t71) * t278 + t287;];
tauc  = t7(:);
