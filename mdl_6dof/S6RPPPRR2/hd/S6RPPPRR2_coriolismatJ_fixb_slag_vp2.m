% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPPRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:42
% EndTime: 2019-03-09 01:31:47
% DurationCPUTime: 3.04s
% Computational Cost: add. (6278->245), mult. (11360->345), div. (0->0), fcn. (11692->8), ass. (0->156)
t171 = sin(qJ(6));
t167 = t171 ^ 2;
t172 = cos(qJ(6));
t168 = t172 ^ 2;
t251 = t167 + t168;
t306 = 0.1e1 - t251;
t236 = m(7) * t306;
t169 = sin(pkin(10));
t170 = cos(pkin(10));
t282 = sin(qJ(5));
t283 = cos(qJ(5));
t148 = t169 * t282 - t170 * t283;
t149 = t169 * t283 + t170 * t282;
t261 = t148 * t149;
t310 = t236 * t261;
t311 = t310 * qJD(5);
t229 = -Ifges(7,5) * t172 + Ifges(7,6) * t171;
t309 = t149 * t229;
t290 = t148 / 0.2e1;
t226 = cos(pkin(9)) * pkin(1) + pkin(2) + qJ(4);
t218 = pkin(7) + t226;
t197 = t218 * t282;
t198 = t218 * t283;
t189 = -t169 * t197 + t170 * t198;
t281 = pkin(5) * t148;
t233 = t149 * pkin(8) - t281;
t182 = t171 * t233 - t172 * t189;
t183 = -t171 * t189 - t172 * t233;
t301 = t171 * t183 + t172 * t182;
t308 = -m(5) / 0.2e1;
t307 = -m(6) / 0.2e1;
t266 = t172 * mrSges(7,1);
t267 = t171 * mrSges(7,2);
t153 = -t266 + t267;
t305 = -t153 + mrSges(6,1);
t265 = t172 * mrSges(7,2);
t268 = t171 * mrSges(7,1);
t154 = t265 + t268;
t304 = t149 * t154;
t123 = t148 ^ 2;
t296 = t149 ^ 2;
t303 = t123 + t296;
t164 = Ifges(7,4) * t172;
t156 = Ifges(7,1) * t171 + t164;
t302 = t169 ^ 2 + t170 ^ 2;
t190 = -t169 * t198 - t170 * t197;
t300 = t148 * t190 + t189 * t149;
t258 = t171 * t148;
t221 = -t149 * mrSges(7,2) + mrSges(7,3) * t258;
t207 = t172 * t221;
t256 = t172 * t148;
t223 = t149 * mrSges(7,1) + mrSges(7,3) * t256;
t211 = t171 * t223;
t299 = t207 / 0.2e1 - t211 / 0.2e1;
t274 = Ifges(7,4) * t171;
t155 = Ifges(7,2) * t172 + t274;
t285 = -t172 / 0.2e1;
t286 = t171 / 0.2e1;
t298 = t155 * t286 + t156 * t285;
t230 = Ifges(7,2) * t171 - t164;
t209 = t172 * t223;
t210 = t171 * t221;
t297 = -t210 / 0.2e1 - t209 / 0.2e1 + t251 * mrSges(7,3) * t290;
t294 = -m(7) / 0.2e1;
t293 = m(7) / 0.2e1;
t237 = t251 * t149;
t292 = m(7) * (-pkin(8) * t237 + t281);
t291 = -t148 / 0.2e1;
t289 = -t149 / 0.2e1;
t288 = t149 / 0.2e1;
t287 = t154 / 0.2e1;
t284 = t172 / 0.2e1;
t280 = pkin(5) * t149;
t278 = mrSges(7,3) * t149;
t276 = Ifges(6,4) * t148;
t273 = Ifges(7,5) * t149;
t271 = Ifges(7,6) * t149;
t269 = Ifges(7,3) * t148;
t113 = t148 * t153;
t219 = -t266 / 0.2e1 + t267 / 0.2e1;
t15 = t113 * t289 + t219 * t261;
t264 = qJD(1) * t15;
t157 = sin(pkin(9)) * pkin(1) + qJ(3);
t248 = t169 * pkin(4) + t157;
t103 = pkin(8) * t148 + t248 + t280;
t186 = -t172 * t103 + t171 * t190;
t241 = -t149 * mrSges(6,1) + t148 * mrSges(6,2);
t50 = t171 * t103 + t172 * t190;
t19 = t210 + t209 + m(7) * (t171 * t50 - t172 * t186) + m(6) * t248 + mrSges(4,3) + t170 * mrSges(5,2) + t169 * mrSges(5,1) - t241 + (m(5) + m(4)) * t157;
t263 = qJD(1) * t19;
t179 = t113 * t290 + t149 * t297;
t12 = t179 + t219;
t262 = t12 * qJD(1);
t260 = t149 * t153;
t257 = t171 * t149;
t222 = mrSges(7,2) * t148 + mrSges(7,3) * t257;
t255 = t172 * t149;
t224 = -mrSges(7,1) * t148 + mrSges(7,3) * t255;
t176 = (t171 * t182 - t172 * t183) * t293 + t222 * t286 + t224 * t284;
t225 = t292 / 0.2e1 + t153 * t291;
t142 = t149 * mrSges(6,2);
t242 = -t148 * mrSges(6,1) - t142;
t16 = -t176 + t225 - t242 - t251 * t278 / 0.2e1;
t259 = t16 * qJD(1);
t254 = t304 * qJD(1);
t191 = (-t251 * t296 - t123) * t293 + t303 * t307 + t302 * t308;
t215 = t251 * t294 + t307 + t308;
t31 = t191 + t215;
t253 = t31 * qJD(1);
t36 = -t310 / 0.2e1;
t252 = t36 * qJD(1);
t245 = t256 / 0.2e1;
t238 = t251 * t148;
t231 = Ifges(7,1) * t172 - t274;
t185 = t171 * t186;
t184 = t148 * t185;
t3 = -t189 * t113 + t186 * t221 + t50 * t223 - mrSges(7,3) * t184 + ((-Ifges(7,4) * t258 - t273) * t171 + (-t50 * mrSges(7,3) - t271 + (t164 + (Ifges(7,1) - Ifges(7,2)) * t171) * t148) * t172) * t148;
t228 = -t3 * qJD(1) - t15 * qJD(2);
t178 = -t148 * t189 - t149 * t185 - t255 * t50;
t217 = t123 * t154;
t9 = t217 - t149 * t207 + t149 * t211 + m(7) * t178 + m(6) * ((-t148 * t198 + t149 * t197) * t170 + (t148 * t197 + t149 * t198) * t169) + t303 * mrSges(6,3) + (m(5) * t226 + mrSges(5,3)) * t302;
t227 = -qJD(1) * t9 - qJD(2) * t36;
t220 = t268 / 0.2e1 + t265 / 0.2e1;
t192 = t172 * (-t148 * t231 + t273);
t193 = t171 * (t148 * t230 + t271);
t195 = -Ifges(7,6) * t148 + t149 * t230;
t196 = -Ifges(7,5) * t148 - t149 * t231;
t1 = t196 * t245 - t195 * t258 / 0.2e1 - t182 * t221 - t50 * t222 + t183 * t223 + t186 * t224 - m(7) * (t182 * t50 + t183 * t186 + t189 * t190) + (-Ifges(6,2) * t149 - t276) * t291 - t248 * t242 + (t229 * t148 + t276 + (-Ifges(6,1) + Ifges(7,3)) * t149) * t290 + (t193 - t269 + t309) * t289 + t300 * t154 + (t192 - 0.2e1 * Ifges(6,4) * t149 + (-Ifges(6,1) + Ifges(6,2)) * t148) * t288;
t208 = t172 * t222;
t212 = t171 * t224;
t7 = t208 * t289 + t212 * t288 + (t149 * t301 - t50 * t256 - t184 + t300) * t294 + (t304 + t299) * t148;
t8 = t212 * t291 + t208 * t290 - t217 / 0.2e1 + t296 * t287 + (-t301 * t148 + t178) * t294 + (t190 * t294 + t299) * t149;
t214 = -t1 * qJD(1) - t8 * qJD(2) - t7 * qJD(3);
t206 = t220 * t149;
t205 = -t154 / 0.2e1 + t220;
t23 = t296 * t236 / 0.2e1 - t306 * t293 * t123;
t204 = -t7 * qJD(1) + t23 * qJD(2) + qJD(3) * t310;
t203 = -t8 * qJD(1) - qJD(2) * t310 + t23 * qJD(3);
t175 = -pkin(5) * t113 / 0.2e1 - t309 / 0.4e1 - t193 / 0.4e1 + t192 / 0.4e1 + t189 * t287 + t155 * t245 - t231 * t256 / 0.4e1 + (0.2e1 * t156 - t230) * t258 / 0.4e1 + t297 * pkin(8);
t177 = -t269 / 0.2e1 - t183 * mrSges(7,1) / 0.2e1 - t182 * mrSges(7,2) / 0.2e1 - Ifges(7,5) * t255 / 0.2e1 + Ifges(7,6) * t257 / 0.2e1;
t5 = t175 - t177;
t64 = pkin(5) * t154 - t171 * t231 / 0.2e1 - t230 * t285 + t298;
t65 = t205 * t149;
t66 = t205 * t148;
t194 = t5 * qJD(1) - t65 * qJD(2) - t66 * qJD(3) - t64 * qJD(5);
t68 = t149 * t287 + t206;
t67 = (t220 + t287) * t148;
t30 = t191 - t215;
t22 = t23 * qJD(5);
t21 = t206 + t299;
t18 = (-t167 / 0.2e1 - t168 / 0.2e1) * t278 + t176 + t225;
t13 = t179 - t219;
t6 = t7 * qJD(5);
t4 = t175 + t177;
t2 = qJD(4) * t36 - qJD(5) * t8 - qJD(6) * t15;
t10 = [qJD(3) * t19 + qJD(4) * t9 - qJD(5) * t1 - qJD(6) * t3, t2, qJD(4) * t30 + qJD(6) * t13 + t263 - t6, qJD(3) * t30 + qJD(5) * t18 + qJD(6) * t21 - t227, t18 * qJD(4) + t4 * qJD(6) + t214 + (pkin(5) * t304 + t196 * t286 + t195 * t284 + (Ifges(7,5) * t171 + Ifges(7,6) * t172) * t291 + Ifges(6,6) * t148 + t189 * mrSges(6,2) + (-m(7) * pkin(5) - t305) * t190 + (-Ifges(6,5) + t298) * t149 + (m(7) * t301 + t208 - t212) * pkin(8) + t301 * mrSges(7,3)) * qJD(5), t13 * qJD(3) + t21 * qJD(4) + t4 * qJD(5) + (-t50 * mrSges(7,1) + mrSges(7,2) * t186 + Ifges(7,5) * t258 + Ifges(7,6) * t256) * qJD(6) + t228; t2, -t311, t22, t252 (-mrSges(7,3) * t237 + t148 * t305 + t142 + t292) * qJD(5) + t68 * qJD(6) + t203, qJD(5) * t68 - qJD(6) * t113 - t264; qJD(4) * t31 + qJD(6) * t12 - t263 - t6, t22, t311, t253 (t260 + m(7) * (-pkin(8) * t238 - t280) - mrSges(7,3) * t238 + t241) * qJD(5) + t67 * qJD(6) + t204, t67 * qJD(5) + qJD(6) * t260 + t262; -qJD(3) * t31 - qJD(5) * t16 - qJD(6) * t304 + t227, -t252, -t253, 0, -t259, -qJD(6) * t154 - t254; qJD(4) * t16 + qJD(6) * t5 - t214, -qJD(6) * t65 - t203, -qJD(6) * t66 - t204, t259, -t64 * qJD(6) (pkin(8) * t153 - t229) * qJD(6) + t194; -qJD(3) * t12 + qJD(4) * t304 - qJD(5) * t5 - t228, qJD(5) * t65 + t264, qJD(5) * t66 - t262, t254, -t194, 0;];
Cq  = t10;
