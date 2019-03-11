% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:53:17
% EndTime: 2019-03-09 02:53:35
% DurationCPUTime: 10.01s
% Computational Cost: add. (6166->511), mult. (14248->736), div. (0->0), fcn. (9876->8), ass. (0->229)
t189 = sin(pkin(10));
t191 = cos(pkin(10));
t192 = sin(qJ(6));
t194 = cos(qJ(6));
t162 = t189 * t194 + t191 * t192;
t153 = t162 * qJD(6);
t190 = sin(pkin(9));
t193 = sin(qJ(3));
t251 = cos(pkin(9));
t274 = cos(qJ(3));
t199 = -t190 * t274 - t193 * t251;
t146 = t199 * qJD(1);
t94 = t162 * t146;
t303 = t153 - t94;
t208 = t189 * t192 - t191 * t194;
t152 = t208 * qJD(6);
t95 = t208 * t146;
t313 = t95 - t152;
t145 = qJD(6) - t146;
t215 = t251 * t274;
t207 = qJD(1) * t215;
t240 = qJD(1) * t193;
t150 = -t190 * t240 + t207;
t124 = qJD(3) * t189 + t150 * t191;
t221 = t191 * qJD(3) - t150 * t189;
t73 = t124 * t194 + t192 * t221;
t273 = Ifges(7,4) * t73;
t311 = -t124 * t192 + t194 * t221;
t21 = Ifges(7,2) * t311 + Ifges(7,6) * t145 + t273;
t292 = t21 / 0.2e1;
t69 = Ifges(7,4) * t311;
t22 = Ifges(7,1) * t73 + Ifges(7,5) * t145 + t69;
t291 = t22 / 0.2e1;
t312 = qJ(2) * (m(3) + m(4)) + mrSges(3,3);
t237 = qJD(1) * qJD(3);
t139 = t199 * t237;
t35 = qJD(6) * t311 - t139 * t208;
t290 = t35 / 0.2e1;
t36 = -qJD(6) * t73 - t139 * t162;
t289 = t36 / 0.2e1;
t225 = t193 * t237;
t138 = -qJD(3) * t207 + t190 * t225;
t283 = -t138 / 0.2e1;
t272 = pkin(3) * t190;
t178 = qJ(5) + t272;
t271 = pkin(8) + t178;
t156 = t271 * t189;
t157 = t271 * t191;
t112 = -t156 * t192 + t157 * t194;
t247 = t146 * t191;
t195 = -pkin(1) - pkin(7);
t173 = qJD(1) * t195 + qJD(2);
t142 = -qJ(4) * t240 + t173 * t193;
t132 = t190 * t142;
t164 = t274 * t173;
t228 = qJD(1) * t274;
t143 = -qJ(4) * t228 + t164;
t101 = t143 * t251 - t132;
t219 = pkin(3) * t228;
t102 = t150 * pkin(4) - t146 * qJ(5) + t219;
t40 = -t101 * t189 + t191 * t102;
t26 = pkin(5) * t150 - pkin(8) * t247 + t40;
t248 = t146 * t189;
t41 = t191 * t101 + t189 * t102;
t29 = -pkin(8) * t248 + t41;
t309 = -qJD(5) * t162 - qJD(6) * t112 + t192 * t29 - t194 * t26;
t111 = -t156 * t194 - t157 * t192;
t308 = -qJD(5) * t208 + qJD(6) * t111 - t192 * t26 - t194 * t29;
t307 = t124 * Ifges(6,5) + t73 * Ifges(7,5) + Ifges(6,6) * t221 + Ifges(7,6) * t311 - t146 * Ifges(6,3) + t145 * Ifges(7,3);
t239 = qJD(3) * t193;
t147 = -qJD(3) * t215 + t190 * t239;
t306 = t208 * qJD(1) + t147 * t162 - t152 * t199;
t305 = -t162 * qJD(1) + t147 * t208 + t153 * t199;
t267 = mrSges(5,3) * t150;
t302 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t221 + mrSges(6,2) * t124 + t267;
t269 = mrSges(5,3) * t139;
t249 = t139 * t191;
t250 = t139 * t189;
t85 = mrSges(6,1) * t250 + mrSges(6,2) * t249;
t301 = t85 + t269;
t241 = qJ(4) - t195;
t187 = qJD(1) * qJD(2);
t227 = qJD(3) * t274;
t217 = qJD(1) * t227;
t165 = pkin(3) * t217 + t187;
t58 = -pkin(4) * t138 - qJ(5) * t139 - qJD(5) * t150 + t165;
t218 = t173 * t227;
t238 = t193 * qJD(4);
t121 = t218 + (-qJ(4) * t227 - t238) * qJD(1);
t226 = t274 * qJD(4);
t231 = t173 * t239;
t197 = -t231 + (qJ(4) * t239 - t226) * qJD(1);
t68 = t251 * t121 + t190 * t197;
t63 = qJD(3) * qJD(5) + t68;
t18 = -t189 * t63 + t191 * t58;
t19 = t189 * t58 + t191 * t63;
t210 = -t18 * t189 + t19 * t191;
t98 = mrSges(6,2) * t138 - mrSges(6,3) * t250;
t99 = -mrSges(6,1) * t138 - mrSges(6,3) * t249;
t298 = -t189 * t99 + t191 * t98;
t90 = mrSges(6,2) * t146 + mrSges(6,3) * t221;
t91 = -mrSges(6,1) * t146 - mrSges(6,3) * t124;
t297 = -t189 * t91 + t191 * t90;
t148 = t199 * qJD(3);
t159 = t190 * t193 - t215;
t67 = t121 * t190 - t251 * t197;
t259 = t159 * t67;
t135 = qJD(3) * pkin(3) + t143;
t87 = t135 * t251 - t132;
t223 = t251 * t142;
t88 = t190 * t135 + t223;
t296 = t147 * t88 - t148 * t87 + t199 * t68 - t259;
t295 = qJD(1) ^ 2;
t294 = Ifges(7,4) * t290 + Ifges(7,2) * t289 + Ifges(7,6) * t283;
t293 = Ifges(7,1) * t290 + Ifges(7,4) * t289 + Ifges(7,5) * t283;
t264 = Ifges(6,4) * t191;
t212 = -Ifges(6,2) * t189 + t264;
t288 = Ifges(6,6) * t283 + t139 * t212 / 0.2e1;
t287 = -t311 / 0.2e1;
t286 = t311 / 0.2e1;
t285 = -t73 / 0.2e1;
t284 = t73 / 0.2e1;
t282 = -t145 / 0.2e1;
t281 = t145 / 0.2e1;
t280 = t146 / 0.2e1;
t279 = -t146 / 0.2e1;
t276 = t150 / 0.2e1;
t80 = qJD(3) * qJ(5) + t88;
t166 = pkin(3) * t240 + qJD(1) * qJ(2) + qJD(4);
t89 = -pkin(4) * t146 - qJ(5) * t150 + t166;
t39 = t189 * t89 + t191 * t80;
t174 = pkin(3) * t227 + qJD(2);
t70 = -pkin(4) * t147 - qJ(5) * t148 + qJD(5) * t159 + t174;
t140 = t239 * t241 - t226;
t168 = t241 * t274;
t141 = -qJD(3) * t168 - t238;
t97 = t190 * t140 + t141 * t251;
t34 = t189 * t70 + t191 * t97;
t270 = mrSges(5,3) * t138;
t268 = mrSges(5,3) * t146;
t266 = Ifges(4,4) * t193;
t265 = Ifges(6,4) * t189;
t167 = t241 * t193;
t119 = -t167 * t190 + t251 * t168;
t263 = t119 * t67;
t260 = t150 * Ifges(5,4);
t246 = t148 * t189;
t245 = t148 * t191;
t244 = t159 * t189;
t243 = t159 * t191;
t171 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t228;
t242 = t193 * t171;
t181 = t193 * pkin(3) + qJ(2);
t114 = -pkin(4) * t199 + qJ(5) * t159 + t181;
t120 = -t167 * t251 - t190 * t168;
t52 = t189 * t114 + t191 * t120;
t28 = -mrSges(7,1) * t311 + mrSges(7,2) * t73;
t236 = t28 + t302;
t235 = Ifges(7,5) * t35 + Ifges(7,6) * t36 - Ifges(7,3) * t138;
t233 = Ifges(4,4) * t274;
t232 = t251 * pkin(3);
t230 = -t250 / 0.2e1;
t229 = t249 / 0.2e1;
t11 = -t36 * mrSges(7,1) + t35 * mrSges(7,2);
t33 = -t189 * t97 + t191 * t70;
t38 = -t189 * t80 + t191 * t89;
t222 = -t138 * mrSges(5,1) + t139 * mrSges(5,2);
t51 = t191 * t114 - t120 * t189;
t96 = -t251 * t140 + t141 * t190;
t100 = t143 * t190 + t223;
t180 = -t232 - pkin(4);
t214 = -mrSges(6,1) * t189 - mrSges(6,2) * t191;
t213 = t191 * Ifges(6,1) - t265;
t211 = Ifges(6,5) * t191 - Ifges(6,6) * t189;
t17 = -pkin(5) * t146 - pkin(8) * t124 + t38;
t27 = pkin(8) * t221 + t39;
t5 = t17 * t194 - t192 * t27;
t6 = t17 * t192 + t194 * t27;
t209 = t189 * t38 - t191 * t39;
t37 = -pkin(5) * t199 + pkin(8) * t243 + t51;
t42 = pkin(8) * t244 + t52;
t12 = -t192 * t42 + t194 * t37;
t13 = t192 * t37 + t194 * t42;
t130 = -qJD(3) * mrSges(5,2) + t268;
t206 = -t130 - t297;
t205 = -Ifges(4,5) * t193 - Ifges(4,6) * t274;
t79 = -qJD(3) * pkin(4) + qJD(5) - t87;
t204 = t79 * t214;
t104 = t162 * t159;
t203 = qJ(2) * (mrSges(4,1) * t274 - mrSges(4,2) * t193);
t202 = t193 * (-Ifges(4,2) * t274 - t266);
t201 = (Ifges(4,1) * t274 - t266) * qJD(1);
t200 = (-Ifges(4,2) * t193 + t233) * qJD(1);
t163 = (t193 * mrSges(4,1) + mrSges(4,2) * t274) * qJD(1);
t198 = (-Ifges(4,1) * t193 - t233) * t274;
t170 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t240;
t169 = -t191 * pkin(5) + t180;
t155 = Ifges(4,5) * qJD(3) + t201;
t154 = Ifges(4,6) * qJD(3) + t200;
t144 = Ifges(5,4) * t146;
t113 = -mrSges(5,1) * t146 + mrSges(5,2) * t150;
t110 = t150 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t144;
t109 = t146 * Ifges(5,2) + Ifges(5,6) * qJD(3) + t260;
t106 = t208 * t159;
t105 = t208 * t199;
t103 = t162 * t199;
t84 = -pkin(5) * t244 + t119;
t64 = pkin(5) * t248 + t100;
t62 = pkin(5) * t246 + t96;
t61 = -t138 * Ifges(6,5) + t139 * t213;
t57 = t124 * Ifges(6,1) + Ifges(6,4) * t221 - t146 * Ifges(6,5);
t56 = t124 * Ifges(6,4) + Ifges(6,2) * t221 - Ifges(6,6) * t146;
t50 = -pkin(5) * t221 + t79;
t49 = mrSges(7,1) * t145 - mrSges(7,3) * t73;
t48 = -mrSges(7,2) * t145 + mrSges(7,3) * t311;
t47 = pkin(5) * t250 + t67;
t46 = -t148 * t162 - t152 * t159;
t44 = qJD(6) * t104 - t148 * t208;
t25 = -pkin(8) * t246 + t34;
t24 = mrSges(7,2) * t138 + mrSges(7,3) * t36;
t23 = -mrSges(7,1) * t138 - mrSges(7,3) * t35;
t16 = -pkin(5) * t147 - pkin(8) * t245 + t33;
t15 = -pkin(8) * t250 + t19;
t14 = -pkin(5) * t138 - pkin(8) * t249 + t18;
t4 = -qJD(6) * t13 + t16 * t194 - t192 * t25;
t3 = qJD(6) * t12 + t16 * t192 + t194 * t25;
t2 = -qJD(6) * t6 + t14 * t194 - t15 * t192;
t1 = qJD(6) * t5 + t14 * t192 + t15 * t194;
t7 = [(t198 + 0.2e1 * t203 - t202) * t237 + (t170 * t227 - t171 * t239) * t195 + (Ifges(7,1) * t106 + Ifges(7,4) * t104 - Ifges(7,5) * t199) * t290 + (t138 * t199 + t147 * t280) * Ifges(5,2) + (Ifges(7,5) * t106 + Ifges(7,6) * t104 - t159 * t211 - (Ifges(6,3) + Ifges(7,3)) * t199) * t283 + (-Ifges(6,5) * t199 - t159 * t213) * t229 + (-Ifges(6,6) * t199 - t159 * t212) * t230 + (-t138 * t159 + t139 * t199 + t147 * t276 + t148 * t280) * Ifges(5,4) + t165 * (-mrSges(5,1) * t199 - mrSges(5,2) * t159) + t1 * (mrSges(7,2) * t199 + mrSges(7,3) * t104) + t2 * (-mrSges(7,1) * t199 - mrSges(7,3) * t106) + t19 * (mrSges(6,2) * t199 + mrSges(6,3) * t244) + t18 * (-mrSges(6,1) * t199 + mrSges(6,3) * t243) + (Ifges(7,4) * t106 + Ifges(7,2) * t104 - Ifges(7,6) * t199) * t289 - (-Ifges(6,3) * t138 + t139 * t211 + t235) * t199 / 0.2e1 - t148 * t204 + t296 * mrSges(5,3) + m(6) * (t18 * t51 + t19 * t52 + t33 * t38 + t34 * t39 + t263) + m(5) * (t120 * t68 + t165 * t181 + t166 * t174 + t88 * t97 + t263) + m(7) * (t1 * t13 + t12 * t2 + t3 * t6 + t4 * t5 + t47 * t84 + t50 * t62) + t38 * (-mrSges(6,1) * t147 - mrSges(6,3) * t245) + t39 * (mrSges(6,2) * t147 - mrSges(6,3) * t246) + (-t139 * t159 + t148 * t276) * Ifges(5,1) + t221 * (-Ifges(6,6) * t147 + t148 * t212) / 0.2e1 + t166 * (-mrSges(5,1) * t147 + mrSges(5,2) * t148) + qJD(3) * (Ifges(5,5) * t148 + Ifges(5,6) * t147) / 0.2e1 + (-Ifges(6,3) * t147 + t148 * t211) * t279 + t124 * (-Ifges(6,5) * t147 + t148 * t213) / 0.2e1 - (t200 + t154) * t227 / 0.2e1 - (t201 + t155) * t239 / 0.2e1 + t301 * t119 + (-m(5) * t87 + m(6) * t79 + t302) * t96 + t174 * t113 - t56 * t246 / 0.2e1 + 0.2e1 * qJD(2) * t163 + t57 * t245 / 0.2e1 - t61 * t243 / 0.2e1 + t148 * t110 / 0.2e1 + t147 * t109 / 0.2e1 + t5 * (-mrSges(7,1) * t147 - mrSges(7,3) * t44) + t6 * (mrSges(7,2) * t147 + mrSges(7,3) * t46) + t97 * t130 - t307 * t147 / 0.2e1 + 0.2e1 * t312 * t187 + t214 * t259 + t120 * t270 + (Ifges(7,5) * t44 + Ifges(7,6) * t46 - Ifges(7,3) * t147) * t281 + (Ifges(7,1) * t44 + Ifges(7,4) * t46 - Ifges(7,5) * t147) * t284 + (Ifges(7,4) * t44 + Ifges(7,2) * t46 - Ifges(7,6) * t147) * t286 + t244 * t288 + t44 * t291 + t46 * t292 + t106 * t293 + t104 * t294 + qJD(3) ^ 2 * t205 / 0.2e1 + t12 * t23 + t13 * t24 + t3 * t48 + t4 * t49 + t50 * (-mrSges(7,1) * t46 + mrSges(7,2) * t44) + t62 * t28 + t181 * t222 + t84 * t11 + t34 * t90 + t33 * t91 + t52 * t98 + t51 * t99 + t47 * (-mrSges(7,1) * t104 + mrSges(7,2) * t106); t103 * t23 + t105 * t24 + t306 * t49 + t305 * t48 + (t170 * t274 - t242) * qJD(3) - t312 * t295 - (t270 + t298) * t199 + (t11 + t301) * t159 - t236 * t148 + t206 * t147 - m(5) * t296 + m(6) * (t147 * t209 - t148 * t79 - t199 * t210 + t259) + (-t163 - t113 - m(5) * t166 - t189 * t90 - t191 * t91 - m(6) * (t39 * t189 + t38 * t191)) * qJD(1) + (t1 * t105 + t103 * t2 - t148 * t50 + t159 * t47 + t305 * t6 + t306 * t5) * m(7); (-t198 / 0.2e1 + t202 / 0.2e1 - t203) * t295 + (Ifges(6,5) * t189 + Ifges(7,5) * t162 + Ifges(6,6) * t191 - Ifges(7,6) * t208) * t283 + t47 * (mrSges(7,1) * t208 + mrSges(7,2) * t162) + (Ifges(7,4) * t162 - Ifges(7,2) * t208) * t289 + (Ifges(7,1) * t162 - Ifges(7,4) * t208) * t290 - t208 * t294 + (-Ifges(7,5) * t95 - Ifges(7,6) * t94 + Ifges(7,3) * t150) * t282 + (-Ifges(7,1) * t95 - Ifges(7,4) * t94 + Ifges(7,5) * t150) * t285 + (-Ifges(7,4) * t95 - Ifges(7,2) * t94 + Ifges(7,6) * t150) * t287 + (-Ifges(5,2) * t150 + t110 + t144) * t279 + t313 * t291 - t303 * t292 + t297 * qJD(5) + t298 * t178 + (t247 * t38 + t248 * t39 + t210) * mrSges(6,3) + (-mrSges(6,1) * t191 + mrSges(6,2) * t189 - mrSges(5,1)) * t67 + (-t1 * t208 - t162 * t2 - t303 * t6 - t313 * t5) * mrSges(7,3) + (mrSges(7,1) * t303 + mrSges(7,2) * t313) * t50 + (-Ifges(7,5) * t152 - Ifges(7,6) * t153) * t281 + (-Ifges(7,1) * t152 - Ifges(7,4) * t153) * t284 + (-Ifges(7,4) * t152 - Ifges(7,2) * t153) * t286 - t221 * (Ifges(6,6) * t150 + t146 * t212) / 0.2e1 - t166 * (mrSges(5,1) * t150 + mrSges(5,2) * t146) - qJD(3) * (Ifges(5,5) * t146 - Ifges(5,6) * t150) / 0.2e1 + (Ifges(6,3) * t150 + t146 * t211) * t280 - t124 * (Ifges(6,5) * t150 + t146 * t213) / 0.2e1 - t302 * t100 - t38 * mrSges(6,1) * t150 + t39 * mrSges(6,2) * t150 + t6 * mrSges(7,2) * t150 - t5 * mrSges(7,1) * t150 + t87 * t268 - t232 * t269 + (Ifges(6,1) * t189 + t264) * t229 + (Ifges(6,2) * t191 + t265) * t230 + (-qJD(5) * t209 - t100 * t79 + t178 * t210 + t180 * t67 - t38 * t40 - t39 * t41) * m(6) + t189 * t61 / 0.2e1 + t180 * t85 + t169 * t11 + t56 * t248 / 0.2e1 - t57 * t247 / 0.2e1 + t155 * t240 / 0.2e1 - t205 * t237 / 0.2e1 + ((t190 * t68 - t251 * t67) * pkin(3) + t87 * t100 - t88 * t101 - t166 * t219) * m(5) + Ifges(5,6) * t138 + Ifges(5,5) * t139 - t101 * t130 - mrSges(4,1) * t231 + t111 * t23 - t170 * t164 + t112 * t24 + t173 * t242 - (Ifges(5,1) * t146 - t260 + t307) * t150 / 0.2e1 + t308 * t48 + t309 * t49 + (t1 * t112 + t111 * t2 + t169 * t47 + t308 * t6 + t309 * t5 - t50 * t64) * m(7) + t88 * t267 + t270 * t272 + t109 * t276 + t191 * t288 + t162 * t293 + t146 * t204 - Ifges(4,6) * t217 - mrSges(4,2) * t218 - t113 * t219 - t64 * t28 - t68 * mrSges(5,2) - Ifges(4,5) * t225 + t154 * t228 / 0.2e1 - t41 * t90 - t40 * t91; -t208 * t23 + t162 * t24 + t189 * t98 + t191 * t99 - t303 * t49 + t313 * t48 - t236 * t150 + t206 * t146 + t222 + (t1 * t162 - t150 * t50 - t2 * t208 - t303 * t5 + t313 * t6) * m(7) + (t146 * t209 - t150 * t79 + t18 * t191 + t189 * t19) * m(6) + (-t146 * t88 + t150 * t87 + t165) * m(5); t124 * t91 - t221 * t90 - t311 * t48 + t73 * t49 + t11 + t85 + (-t311 * t6 + t5 * t73 + t47) * m(7) + (t124 * t38 - t221 * t39 + t67) * m(6); -t1 * mrSges(7,2) + t2 * mrSges(7,1) - t50 * (mrSges(7,1) * t73 + mrSges(7,2) * t311) + (Ifges(7,1) * t311 - t273) * t285 + t21 * t284 + (Ifges(7,5) * t311 - Ifges(7,6) * t73) * t282 - t5 * t48 + t6 * t49 + (t311 * t5 + t6 * t73) * mrSges(7,3) + t235 + (-Ifges(7,2) * t73 + t22 + t69) * t287;];
tauc  = t7(:);
