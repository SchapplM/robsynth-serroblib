% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:48
% EndTime: 2019-12-05 18:37:59
% DurationCPUTime: 4.87s
% Computational Cost: add. (7403->409), mult. (19875->575), div. (0->0), fcn. (14502->8), ass. (0->200)
t221 = sin(qJ(5));
t224 = cos(qJ(5));
t222 = sin(qJ(3));
t223 = sin(qJ(2));
t225 = cos(qJ(3));
t226 = cos(qJ(2));
t195 = -t222 * t223 + t225 * t226;
t183 = t195 * qJD(1);
t196 = t222 * t226 + t225 * t223;
t184 = t196 * qJD(1);
t219 = sin(pkin(9));
t220 = cos(pkin(9));
t141 = t183 * t219 + t184 * t220;
t133 = pkin(8) * t141;
t218 = qJD(2) + qJD(3);
t271 = -pkin(7) - pkin(6);
t209 = t271 * t226;
t201 = qJD(1) * t209;
t185 = t222 * t201;
t208 = t271 * t223;
t200 = qJD(1) * t208;
t191 = qJD(2) * pkin(2) + t200;
t149 = t225 * t191 + t185;
t179 = t184 * qJ(4);
t119 = t149 - t179;
t108 = pkin(3) * t218 + t119;
t188 = t225 * t201;
t150 = t191 * t222 - t188;
t246 = qJ(4) * t183;
t120 = t150 + t246;
t109 = t219 * t120;
t59 = t220 * t108 - t109;
t40 = pkin(4) * t218 - t133 + t59;
t230 = t220 * t183 - t184 * t219;
t256 = pkin(8) * t230;
t242 = t220 * t120;
t60 = t219 * t108 + t242;
t42 = t60 + t256;
t10 = t221 * t40 + t224 * t42;
t84 = t141 * t224 + t221 * t230;
t255 = t10 * t84;
t153 = -t200 * t222 + t188;
t122 = t153 - t246;
t154 = t225 * t200 + t185;
t123 = -t179 + t154;
t241 = t220 * t222;
t251 = pkin(2) * qJD(3);
t279 = -t220 * t122 + t123 * t219 + (-t219 * t225 - t241) * t251;
t243 = t219 * t222;
t278 = -t219 * t122 - t220 * t123 + (t220 * t225 - t243) * t251;
t292 = t256 + t279;
t291 = t133 + t278;
t290 = t218 / 0.2e1;
t268 = t230 / 0.2e1;
t289 = t141 * Ifges(5,4);
t288 = -t141 * t221 + t224 * t230;
t287 = t289 / 0.2e1 + t268 * Ifges(5,2) + Ifges(5,6) * t290;
t212 = pkin(2) * t225 + pkin(3);
t175 = -pkin(2) * t243 + t220 * t212;
t168 = pkin(4) + t175;
t177 = pkin(2) * t241 + t212 * t219;
t131 = t168 * t221 + t177 * t224;
t285 = -qJD(5) * t131 - t291 * t221 + t224 * t292;
t130 = t168 * t224 - t177 * t221;
t284 = qJD(5) * t130 + t221 * t292 + t291 * t224;
t283 = mrSges(5,3) * t230;
t282 = Ifges(5,4) * t230;
t211 = pkin(3) * t220 + pkin(4);
t257 = pkin(3) * t219;
t178 = t211 * t221 + t224 * t257;
t61 = -t119 * t219 - t242;
t43 = t61 - t256;
t62 = t220 * t119 - t109;
t44 = -t133 + t62;
t281 = -t178 * qJD(5) + t221 * t44 - t224 * t43;
t176 = t211 * t224 - t221 * t257;
t280 = t176 * qJD(5) - t221 * t43 - t224 * t44;
t159 = t222 * t208 - t225 * t209;
t275 = -t288 / 0.2e1;
t274 = t288 / 0.2e1;
t273 = -t84 / 0.2e1;
t272 = t84 / 0.2e1;
t270 = (pkin(1) * mrSges(3,1));
t269 = (pkin(1) * mrSges(3,2));
t267 = t141 / 0.2e1;
t265 = t183 / 0.2e1;
t264 = t184 / 0.2e1;
t217 = qJD(5) + t218;
t263 = -t217 / 0.2e1;
t213 = -pkin(2) * t226 - pkin(1);
t207 = qJD(1) * t213;
t260 = m(4) * t207;
t259 = Ifges(6,4) * t84;
t258 = pkin(3) * t184;
t156 = t218 * t196;
t146 = t156 * qJD(1);
t235 = qJD(2) * t271;
t229 = qJD(1) * t235;
t192 = t223 * t229;
t193 = t226 * t229;
t236 = qJD(3) * t225;
t237 = qJD(3) * t222;
t95 = t191 * t236 + t225 * t192 + t222 * t193 + t201 * t237;
t50 = -qJ(4) * t146 + qJD(4) * t183 + t95;
t155 = t218 * t195;
t145 = t155 * qJD(1);
t96 = -qJD(3) * t150 - t192 * t222 + t225 * t193;
t51 = -qJ(4) * t145 - qJD(4) * t184 + t96;
t21 = t219 * t51 + t220 * t50;
t202 = t223 * t235;
t203 = t226 * t235;
t103 = t225 * t202 + t222 * t203 + t208 * t236 + t209 * t237;
t67 = -qJ(4) * t156 + qJD(4) * t195 + t103;
t104 = -t159 * qJD(3) - t202 * t222 + t225 * t203;
t68 = -qJ(4) * t155 - qJD(4) * t196 + t104;
t34 = t219 * t68 + t220 * t67;
t254 = mrSges(4,3) * t183;
t253 = mrSges(4,3) * t184;
t252 = Ifges(3,4) * t223;
t250 = t184 * Ifges(4,4);
t249 = t60 * t141;
t248 = Ifges(3,5) * qJD(2);
t247 = Ifges(3,6) * qJD(2);
t245 = qJD(2) * mrSges(3,1);
t244 = qJD(2) * mrSges(3,2);
t158 = t225 * t208 + t209 * t222;
t136 = -qJ(4) * t196 + t158;
t137 = qJ(4) * t195 + t159;
t74 = t219 * t136 + t220 * t137;
t240 = qJD(1) * t223;
t239 = qJD(1) * t226;
t238 = qJD(2) * t223;
t215 = pkin(2) * t240;
t234 = t248 / 0.2e1;
t233 = -t247 / 0.2e1;
t90 = -t145 * t219 - t146 * t220;
t91 = t145 * t220 - t146 * t219;
t232 = -t90 * mrSges(5,1) + t91 * mrSges(5,2);
t27 = qJD(5) * t288 + t221 * t90 + t224 * t91;
t28 = -qJD(5) * t84 - t221 * t91 + t224 * t90;
t231 = -t28 * mrSges(6,1) + t27 * mrSges(6,2);
t126 = pkin(3) * t146 + qJD(2) * t215;
t142 = pkin(2) * t238 + pkin(3) * t156;
t20 = -t219 * t50 + t220 * t51;
t33 = -t219 * t67 + t220 * t68;
t73 = t220 * t136 - t137 * t219;
t105 = pkin(4) * t141 + t258;
t7 = -pkin(8) * t91 + t20;
t8 = pkin(8) * t90 + t21;
t9 = -t221 * t42 + t224 * t40;
t2 = qJD(5) * t9 + t221 * t7 + t224 * t8;
t3 = -qJD(5) * t10 - t221 * t8 + t224 * t7;
t228 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t27 + Ifges(6,6) * t28;
t152 = t195 * t219 + t196 * t220;
t55 = -pkin(8) * t152 + t73;
t151 = t195 * t220 - t196 * t219;
t56 = pkin(8) * t151 + t74;
t29 = -t221 * t56 + t224 * t55;
t30 = t221 * t55 + t224 * t56;
t93 = t151 * t224 - t152 * t221;
t94 = t151 * t221 + t152 * t224;
t165 = -t195 * pkin(3) + t213;
t157 = -t183 * pkin(3) + qJD(4) + t207;
t134 = t183 * Ifges(4,2) + t218 * Ifges(4,6) + t250;
t180 = Ifges(4,4) * t183;
t135 = t184 * Ifges(4,1) + t218 * Ifges(4,5) + t180;
t36 = Ifges(6,2) * t288 + Ifges(6,6) * t217 + t259;
t80 = Ifges(6,4) * t288;
t37 = Ifges(6,1) * t84 + Ifges(6,5) * t217 + t80;
t79 = t141 * Ifges(5,1) + t218 * Ifges(5,5) + t282;
t99 = -pkin(4) * t230 + t157;
t227 = t20 * mrSges(5,1) - t21 * mrSges(5,2) + t141 * t287 - t141 * (Ifges(5,1) * t230 - t289) / 0.2e1 - t157 * (mrSges(5,1) * t141 + mrSges(5,2) * t230) + Ifges(4,5) * t145 - Ifges(4,6) * t146 + Ifges(5,5) * t91 - t95 * mrSges(4,2) + t96 * mrSges(4,1) + Ifges(5,6) * t90 - t207 * (mrSges(4,1) * t184 + mrSges(4,2) * t183) + t228 + t59 * t283 - t184 * (Ifges(4,1) * t183 - t250) / 0.2e1 + t149 * t254 + t134 * t264 - (-t36 / 0.2e1 + t99 * mrSges(6,1) + Ifges(6,6) * t263 + Ifges(6,4) * t273 + Ifges(6,2) * t275) * t84 + (-t37 / 0.2e1 - t99 * mrSges(6,2) + t9 * mrSges(6,3) + Ifges(6,5) * t263 + Ifges(6,1) * t273 + Ifges(6,4) * t275) * t288 - (-Ifges(5,2) * t141 + t282 + t79) * t230 / 0.2e1 - (-Ifges(4,2) * t184 + t135 + t180) * t183 / 0.2e1 - (Ifges(4,5) * t183 + Ifges(5,5) * t230 - Ifges(4,6) * t184 - Ifges(5,6) * t141) * t218 / 0.2e1;
t214 = Ifges(3,4) * t239;
t205 = mrSges(3,3) * t239 - t244;
t204 = -mrSges(3,3) * t240 + t245;
t182 = Ifges(3,1) * t240 + t214 + t248;
t181 = t247 + (Ifges(3,2) * t226 + t252) * qJD(1);
t162 = mrSges(4,1) * t218 - t253;
t161 = -mrSges(4,2) * t218 + t254;
t160 = t215 + t258;
t148 = -mrSges(4,1) * t183 + mrSges(4,2) * t184;
t125 = mrSges(5,1) * t218 - mrSges(5,3) * t141;
t124 = -mrSges(5,2) * t218 + t283;
t113 = -t151 * pkin(4) + t165;
t100 = t105 + t215;
t98 = t155 * t220 - t156 * t219;
t97 = -t155 * t219 - t156 * t220;
t86 = -mrSges(5,1) * t230 + mrSges(5,2) * t141;
t72 = mrSges(6,1) * t217 - mrSges(6,3) * t84;
t71 = -mrSges(6,2) * t217 + mrSges(6,3) * t288;
t66 = -pkin(4) * t97 + t142;
t54 = -pkin(4) * t90 + t126;
t39 = -mrSges(6,1) * t288 + mrSges(6,2) * t84;
t32 = -qJD(5) * t94 - t221 * t98 + t224 * t97;
t31 = qJD(5) * t93 + t221 * t97 + t224 * t98;
t16 = pkin(8) * t97 + t34;
t15 = -pkin(8) * t98 + t33;
t5 = -qJD(5) * t30 + t15 * t224 - t16 * t221;
t4 = qJD(5) * t29 + t15 * t221 + t16 * t224;
t1 = [m(5) * (t126 * t165 + t142 * t157 + t20 * t73 + t21 * t74 + t33 * t59 + t34 * t60) + m(6) * (t10 * t4 + t113 * t54 + t2 * t30 + t29 * t3 + t5 * t9 + t66 * t99) + t97 * t287 + t155 * t135 / 0.2e1 - t156 * t134 / 0.2e1 + t157 * (-mrSges(5,1) * t97 + mrSges(5,2) * t98) + t103 * t161 + t104 * t162 + t126 * (-mrSges(5,1) * t151 + mrSges(5,2) * t152) + t142 * t86 + t34 * t124 + t33 * t125 + t98 * t79 / 0.2e1 + t99 * (-mrSges(6,1) * t32 + mrSges(6,2) * t31) + t54 * (-mrSges(6,1) * t93 + mrSges(6,2) * t94) + t4 * t71 + t5 * t72 + t66 * t39 + t32 * t36 / 0.2e1 + t31 * t37 / 0.2e1 + m(4) * (t103 * t150 + t104 * t149 + t158 * t96 + t159 * t95) + t207 * (mrSges(4,1) * t156 + mrSges(4,2) * t155) + (-t195 * t146 - t156 * t265) * Ifges(4,2) + (t195 * t145 - t196 * t146 + t155 * t265 - t156 * t264) * Ifges(4,4) + (-t145 * t158 - t146 * t159 - t149 * t155 - t150 * t156 + t195 * t95 - t196 * t96) * mrSges(4,3) + t213 * (mrSges(4,1) * t146 + mrSges(4,2) * t145) + (t10 * t32 + t2 * t93 - t27 * t29 + t28 * t30 - t3 * t94 - t31 * t9) * mrSges(6,3) + (t151 * t21 - t152 * t20 - t59 * t98 + t60 * t97 - t73 * t91 + t74 * t90) * mrSges(5,3) + (-pkin(6) * t204 + t182 / 0.2e1 + t234 + (-(2 * t269) + 0.3e1 / 0.2e1 * Ifges(3,4) * t226) * qJD(1)) * t226 * qJD(2) + t113 * t231 + t165 * t232 + t217 * (Ifges(6,5) * t31 + Ifges(6,6) * t32) / 0.2e1 + (Ifges(4,5) * t155 + Ifges(5,5) * t98 - Ifges(4,6) * t156 + Ifges(5,6) * t97) * t290 + (t196 * t145 + t155 * t264) * Ifges(4,1) + (t91 * t152 + t267 * t98) * Ifges(5,1) + (t151 * t91 + t90 * t152 + t267 * t97 + t268 * t98) * Ifges(5,4) + (t151 * t90 + t268 * t97) * Ifges(5,2) + (-pkin(6) * t205 - t181 / 0.2e1 + t233 + (-(2 * t270) - 0.3e1 / 0.2e1 * t252 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t226) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t195 + mrSges(4,2) * t196) + t148 + 0.2e1 * t260) * pkin(2)) * t238 + (t94 * t27 + t272 * t31) * Ifges(6,1) + (t274 * t32 + t93 * t28) * Ifges(6,2) + (t93 * t27 + t272 * t32 + t274 * t31 + t94 * t28) * Ifges(6,4); (-t130 * t27 + t131 * t28 + t255) * mrSges(6,3) - t160 * t86 - t154 * t161 - t153 * t162 - t100 * t39 + (-t175 * t91 + t177 * t90 + t249) * mrSges(5,3) + t150 * t253 + t227 + ((t161 * t225 - t162 * t222) * qJD(3) + (-t145 * t225 - t146 * t222) * mrSges(4,3)) * pkin(2) + t285 * t72 + t284 * t71 + t279 * t125 + t278 * t124 + ((t234 - t182 / 0.2e1 - t214 / 0.2e1 + qJD(1) * t269 + (t204 - t245) * pkin(6)) * t226 + (t233 + t181 / 0.2e1 + (t270 + t252 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t226) * qJD(1) + (t205 + t244) * pkin(6) + (-t148 - t260) * pkin(2)) * t223) * qJD(1) + (t10 * t284 - t99 * t100 + t3 * t130 + t2 * t131 + t285 * t9) * m(6) + (-t157 * t160 + t20 * t175 + t21 * t177 + t278 * t60 + t279 * t59) * m(5) + (-t149 * t153 - t150 * t154 + (t222 * t95 + t225 * t96 + (-t149 * t222 + t150 * t225) * qJD(3)) * pkin(2)) * m(4); -t149 * t161 - t62 * t124 - t61 * t125 - t105 * t39 - t86 * t258 + t227 + (t249 + (t219 * t90 - t220 * t91) * pkin(3)) * mrSges(5,3) + t281 * t72 + t280 * t71 + (t162 + t253) * t150 + (-t176 * t27 + t178 * t28 + t255) * mrSges(6,3) + (t10 * t280 - t105 * t99 + t176 * t3 + t178 * t2 + t281 * t9) * m(6) + (-t157 * t258 - t59 * t61 - t60 * t62 + (t20 * t220 + t21 * t219) * pkin(3)) * m(5); -t230 * t124 + t141 * t125 - t288 * t71 + t84 * t72 + t231 + t232 + (-t10 * t288 + t84 * t9 + t54) * m(6) + (t141 * t59 - t230 * t60 + t126) * m(5); -t99 * (mrSges(6,1) * t84 + mrSges(6,2) * t288) + (Ifges(6,1) * t288 - t259) * t273 + t36 * t272 + (Ifges(6,5) * t288 - Ifges(6,6) * t84) * t263 - t9 * t71 + t10 * t72 + (t288 * t9 + t255) * mrSges(6,3) + t228 + (-Ifges(6,2) * t84 + t37 + t80) * t275;];
tauc = t1(:);
