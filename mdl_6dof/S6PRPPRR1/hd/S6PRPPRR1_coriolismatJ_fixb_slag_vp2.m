% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:14:07
% EndTime: 2019-03-08 19:14:13
% DurationCPUTime: 2.72s
% Computational Cost: add. (8574->303), mult. (20293->440), div. (0->0), fcn. (23187->12), ass. (0->190)
t192 = sin(pkin(12));
t193 = cos(pkin(12));
t266 = sin(pkin(11));
t318 = pkin(2) * t266;
t240 = qJ(4) + t318;
t321 = (m(5) * t240 + mrSges(5,3)) * (t192 ^ 2 + t193 ^ 2);
t267 = sin(pkin(6));
t238 = t267 * t266;
t268 = cos(pkin(11));
t239 = t268 * t267;
t295 = sin(qJ(2));
t297 = cos(qJ(2));
t155 = t238 * t295 - t239 * t297;
t294 = sin(qJ(5));
t296 = cos(qJ(5));
t176 = -t192 * t296 - t193 * t294;
t101 = t176 * t155;
t174 = t192 * t294 - t193 * t296;
t102 = t174 * t155;
t195 = cos(qJ(6));
t194 = sin(qJ(6));
t204 = t238 * t297 + t239 * t295;
t82 = t102 * t195 + t194 * t204;
t270 = t82 * t195;
t81 = -t102 * t194 + t195 * t204;
t271 = t81 * t194;
t232 = -t270 + t271;
t309 = -m(7) / 0.2e1;
t320 = (t270 / 0.2e1 - t271 / 0.2e1) * mrSges(7,3) - (-pkin(5) * t101 - pkin(9) * t232) * t309 - t102 * mrSges(6,2) / 0.2e1;
t317 = pkin(2) * t268;
t190 = t194 ^ 2;
t191 = t195 ^ 2;
t241 = mrSges(7,3) * (t191 / 0.2e1 + t190 / 0.2e1);
t186 = Ifges(7,5) * t195;
t285 = Ifges(7,6) * t194;
t223 = -t186 / 0.2e1 + t285 / 0.2e1;
t316 = Ifges(6,4) + t223;
t178 = -mrSges(7,1) * t195 + mrSges(7,2) * t194;
t315 = -t178 + mrSges(6,1);
t275 = t195 * mrSges(7,2);
t279 = t194 * mrSges(7,1);
t179 = t275 + t279;
t133 = t179 * t176;
t314 = t176 * mrSges(6,3) + t133;
t187 = Ifges(7,4) * t195;
t181 = Ifges(7,1) * t194 + t187;
t228 = pkin(8) + t240;
t167 = t228 * t193;
t215 = t192 * t228;
t125 = t167 * t294 + t215 * t296;
t292 = pkin(5) * t176;
t142 = pkin(9) * t174 - t292;
t83 = t125 * t194 + t142 * t195;
t84 = -t125 * t195 + t142 * t194;
t231 = -t194 * t83 + t195 * t84;
t236 = Ifges(7,2) * t194 - t187;
t311 = -m(7) * pkin(5) - t315;
t310 = m(6) / 0.2e1;
t308 = m(7) / 0.2e1;
t269 = cos(pkin(6));
t143 = t192 * t269 + t193 * t204;
t201 = t192 * t204 - t193 * t269;
t89 = t143 * t296 - t201 * t294;
t65 = t155 * t195 - t194 * t89;
t307 = t65 / 0.2e1;
t260 = t176 * t195;
t261 = t176 * t194;
t131 = mrSges(7,1) * t260 - mrSges(7,2) * t261;
t306 = -t131 / 0.2e1;
t140 = t174 * mrSges(7,1) + mrSges(7,3) * t260;
t305 = -t140 / 0.2e1;
t304 = -t178 / 0.2e1;
t303 = t179 / 0.2e1;
t302 = -t194 / 0.2e1;
t301 = t194 / 0.2e1;
t299 = -t195 / 0.2e1;
t298 = t195 / 0.2e1;
t251 = -t190 - t191;
t244 = t251 * t174;
t116 = pkin(9) * t244 + t292;
t293 = m(7) * t116;
t291 = mrSges(7,3) * t174;
t289 = Ifges(7,4) * t194;
t288 = Ifges(7,5) * t174;
t286 = Ifges(7,6) * t174;
t88 = t143 * t294 + t201 * t296;
t284 = t101 * t88;
t282 = t174 * mrSges(6,3);
t280 = t176 * t88;
t278 = t194 * t65;
t276 = t194 * t84;
t66 = t155 * t194 + t195 * t89;
t274 = t195 * t66;
t273 = t195 * t83;
t262 = t174 * t101;
t27 = (t176 * t232 + t262) * t308 + (-t102 * t176 + t262) * t310;
t265 = qJD(1) * t27;
t264 = t125 * t101;
t263 = t125 * t176;
t139 = -mrSges(7,1) * t176 + t195 * t291;
t259 = t194 * t139;
t258 = t194 * t140;
t137 = mrSges(7,2) * t176 + t194 * t291;
t257 = t195 * t137;
t138 = -mrSges(7,2) * t174 + mrSges(7,3) * t261;
t256 = t195 * t138;
t26 = t174 * t306 + (t138 * t301 + t140 * t298 - t176 * t241) * t176;
t255 = t26 * qJD(2);
t168 = t174 * mrSges(6,2);
t219 = t137 * t302 + t139 * t299;
t227 = t174 * t241;
t28 = t168 + (t304 + mrSges(6,1)) * t176 - t227 + 0.2e1 * (t116 / 0.4e1 - t276 / 0.4e1 - t273 / 0.4e1) * m(7) + t219;
t254 = t28 * qJD(2);
t222 = t279 / 0.2e1 + t275 / 0.2e1;
t212 = t222 * t174;
t218 = t258 / 0.2e1 - t256 / 0.2e1;
t33 = t212 + t218;
t253 = t33 * qJD(2);
t243 = (-0.1e1 - t251) * t176;
t50 = m(7) * t174 * t243;
t49 = t50 / 0.2e1;
t252 = t49 * qJD(2);
t126 = t167 * t296 - t215 * t294;
t132 = t174 * t179;
t242 = -pkin(3) - t317;
t177 = -pkin(4) * t193 + t242;
t120 = pkin(5) * t174 + pkin(9) * t176 + t177;
t70 = t120 * t195 - t126 * t194;
t71 = t120 * t194 + t126 * t195;
t233 = t194 * t70 - t195 * t71;
t203 = (t126 + t233) * t308 - t132 / 0.2e1 + t218;
t11 = ((-t125 - t231) * t308 - t257 / 0.2e1 + t259 / 0.2e1 + t133 / 0.2e1) * t176 + t203 * t174;
t250 = qJD(4) * t49 + qJD(5) * t11 + qJD(6) * t26;
t237 = Ifges(7,1) * t195 - t289;
t180 = Ifges(7,2) * t195 + t289;
t235 = Ifges(7,5) * t194 + Ifges(7,6) * t195;
t234 = -t274 + t278;
t121 = t155 * t204;
t199 = t193 * t143 + t192 * t201;
t12 = m(7) * (t65 * t81 + t66 * t82 + t284) + m(6) * (t102 * t89 + t121 + t284) + m(5) * (-t155 * t199 + t121);
t230 = t12 * qJD(1) + t27 * qJD(3);
t224 = t234 + t89;
t13 = (t174 * t224 + t243 * t88) * t308;
t17 = m(7) * t224 * t88;
t229 = t17 * qJD(1) + t13 * qJD(3);
t226 = t81 * mrSges(7,1) / 0.2e1 - t82 * mrSges(7,2) / 0.2e1;
t225 = -t83 * mrSges(7,1) / 0.2e1 + t84 * mrSges(7,2) / 0.2e1;
t113 = t176 * t236 + t286;
t135 = t176 * t181;
t221 = -t113 / 0.4e1 + t135 / 0.4e1 - pkin(9) * t138 / 0.2e1;
t115 = -t176 * t237 + t288;
t134 = t176 * t180;
t220 = pkin(9) * t305 + t115 / 0.4e1 + t134 / 0.4e1;
t217 = t180 * t301 + t181 * t299;
t216 = t176 * t235;
t141 = -t176 * mrSges(6,1) - t168;
t198 = t203 * t88 + (t125 * t89 + t83 * t65 + t84 * t66) * t308 + t155 * t141 / 0.2e1 + t139 * t307 + t66 * t137 / 0.2e1 - t89 * t133 / 0.2e1;
t2 = (t304 + mrSges(6,1) / 0.2e1) * t101 + t198 - t320;
t112 = -Ifges(7,6) * t176 + t174 * t236;
t114 = -Ifges(7,5) * t176 - t174 * t237;
t3 = t83 * t140 + t70 * t139 + t84 * t138 + t71 * t137 + m(7) * (t125 * t126 + t70 * t83 + t71 * t84) - t126 * t133 - t125 * t132 + t177 * t141 + (t112 * t301 + t114 * t299 - t316 * t176) * t176 + (t113 * t301 + t115 * t299 + (-Ifges(7,3) - Ifges(6,2) + Ifges(6,1)) * t176 + t316 * t174) * t174;
t214 = t2 * qJD(1) + t3 * qJD(2) + t11 * qJD(3);
t202 = (t274 / 0.2e1 - t278 / 0.2e1) * t176 * mrSges(7,3) + t138 * t307 + t66 * t305 + t88 * t306;
t5 = t202 - t226;
t7 = -t125 * t131 + t70 * t138 - t71 * t140 + (t174 * t235 / 0.2e1 + t113 * t298 + t135 * t299 - t233 * mrSges(7,3) + (t134 + t115) * t301) * t176;
t213 = t5 * qJD(1) + t7 * qJD(2) + t26 * qJD(3);
t211 = -t179 / 0.2e1 + t222;
t210 = t13 * qJD(1) + t11 * qJD(2) + t50 * qJD(3);
t196 = (t174 * t234 - t280) * t308 + (-t174 * t89 - t280) * t310 + m(5) * t199 / 0.2e1;
t197 = (t194 * t82 + t195 * t81) * t309 + (-m(5) / 0.2e1 - m(6) / 0.2e1) * t204;
t16 = t196 + t197;
t21 = t314 * t176 + (-t256 + t258 + t282) * t174 + m(7) * (t174 * t233 - t263) + m(6) * (-t174 * t126 - t263) + t321;
t209 = -qJD(1) * t16 - qJD(2) * t21 - qJD(3) * t49;
t208 = pkin(5) * t131 / 0.2e1 + t125 * t303 + t174 * t186 / 0.4e1;
t19 = t211 * t88;
t200 = pkin(9) * t241 + (-t236 + t181) * t194 / 0.4e1 + (t180 / 0.4e1 - t237 / 0.4e1) * t195;
t8 = (t288 / 0.2e1 + t220) * t195 + (-0.3e1 / 0.4e1 * t286 + t221) * t194 + (Ifges(7,3) / 0.2e1 + t200) * t176 + t208 + t225;
t90 = pkin(5) * t179 - t236 * t299 + t237 * t302 + t217;
t92 = t211 * t174;
t207 = t19 * qJD(1) - t8 * qJD(2) + t92 * qJD(3) + t90 * qJD(5);
t93 = t174 * t303 + t212;
t34 = t212 - t218;
t30 = (t273 + t276) * t308 + t293 / 0.2e1 + t176 * t304 - t227 - t219;
t20 = (t222 + t303) * t88;
t15 = t196 - t197;
t9 = (-t286 / 0.4e1 + t221) * t194 + t220 * t195 + t208 - t225 + (-Ifges(7,3) / 0.2e1 + t200) * t176 + t223 * t174;
t6 = t202 + t226;
t4 = qJD(2) * t27 + qJD(5) * t13;
t1 = t198 + (t178 / 0.2e1 - mrSges(6,1) / 0.2e1) * t101 + t320;
t10 = [qJD(2) * t12 + qJD(5) * t17 (m(7) * (t70 * t81 + t71 * t82 + t264) + t82 * t138 + t81 * t140 + m(6) * (t126 * t102 + t264) - t102 * t282 + (-mrSges(3,1) * t295 - mrSges(3,2) * t297) * t267 + (-m(4) * t317 + m(5) * t242 + m(6) * t177 - mrSges(5,1) * t193 + mrSges(6,1) * t174 + mrSges(5,2) * t192 - mrSges(6,2) * t176 - mrSges(4,1)) * t204 - t314 * t101 + (-m(4) * t318 + mrSges(4,2) - t321) * t155) * qJD(2) + t15 * qJD(4) + t1 * qJD(5) + t6 * qJD(6) + t230, t4, qJD(2) * t15, t1 * qJD(2) + (t311 * t89 + (mrSges(6,2) + (m(7) * pkin(9) + mrSges(7,3)) * t251) * t88) * qJD(5) + t20 * qJD(6) + t229, t6 * qJD(2) + t20 * qJD(5) + (-mrSges(7,1) * t66 - mrSges(7,2) * t65) * qJD(6); qJD(4) * t16 + qJD(5) * t2 + qJD(6) * t5 - t230, qJD(4) * t21 + qJD(5) * t3 + qJD(6) * t7, t250 - t265, qJD(5) * t30 + qJD(6) * t34 - t209, t30 * qJD(4) + t9 * qJD(6) + t214 + (-t216 / 0.2e1 + t114 * t301 + t112 * t298 + pkin(5) * t132 + Ifges(6,6) * t176 + t125 * mrSges(6,2) + (-Ifges(6,5) + t217) * t174 + t311 * t126 + (m(7) * t231 + t257 - t259) * pkin(9) + t231 * mrSges(7,3)) * qJD(5), t34 * qJD(4) + t9 * qJD(5) + (-t71 * mrSges(7,1) - t70 * mrSges(7,2) + t216) * qJD(6) + t213; t4, t250 + t265, t50 * qJD(5), t252 (mrSges(7,3) * t244 + t176 * t315 + t168 + t293) * qJD(5) + t93 * qJD(6) + t210, qJD(5) * t93 + qJD(6) * t131 + t255; -qJD(2) * t16, -qJD(5) * t28 - qJD(6) * t33 + t209, -t252, 0, -t254, -qJD(6) * t179 - t253; -qJD(2) * t2 - qJD(6) * t19 - t229, qJD(4) * t28 + qJD(6) * t8 - t214, -qJD(6) * t92 - t210, t254, -t90 * qJD(6) (pkin(9) * t178 + t186 - t285) * qJD(6) - t207; -t5 * qJD(2) + t19 * qJD(5), qJD(4) * t33 - qJD(5) * t8 - t213, qJD(5) * t92 - t255, t253, t207, 0;];
Cq  = t10;
