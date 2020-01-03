% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:36
% EndTime: 2019-12-31 18:44:42
% DurationCPUTime: 2.54s
% Computational Cost: add. (3498->352), mult. (7508->486), div. (0->0), fcn. (5837->6), ass. (0->182)
t306 = Ifges(6,4) + Ifges(5,5);
t181 = sin(qJ(4));
t178 = t181 ^ 2;
t183 = cos(qJ(4));
t179 = t183 ^ 2;
t309 = t178 + t179;
t184 = cos(qJ(3));
t218 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t205 = t218 * t184;
t182 = sin(qJ(3));
t173 = Ifges(6,5) * t181;
t301 = Ifges(6,1) * t183 + t173;
t87 = -Ifges(6,4) * t184 + t182 * t301;
t261 = Ifges(5,4) * t181;
t144 = Ifges(5,1) * t183 - t261;
t89 = -Ifges(5,5) * t184 + t182 * t144;
t308 = -t205 + t87 / 0.2e1 + t89 / 0.2e1;
t277 = t182 / 0.2e1;
t307 = -mrSges(5,1) - mrSges(6,1);
t305 = Ifges(6,2) + Ifges(5,3);
t166 = sin(pkin(8)) * pkin(1) + pkin(6);
t234 = t166 * t183;
t206 = -qJ(5) + t234;
t214 = -cos(pkin(8)) * pkin(1) - pkin(2);
t116 = -pkin(3) * t184 - t182 * pkin(7) + t214;
t233 = t181 * t116;
t34 = t184 * t206 + t233;
t226 = t184 * t166;
t40 = t183 * t226 + t233;
t304 = -t34 + t40;
t207 = t166 * t181 + pkin(4);
t229 = t183 * t116;
t35 = t184 * t207 - t229;
t225 = t184 * t181;
t39 = -t166 * t225 + t229;
t303 = t35 + t39;
t200 = pkin(4) * t183 + qJ(5) * t181;
t130 = -pkin(3) - t200;
t203 = t183 * mrSges(6,1) + t181 * mrSges(6,3);
t302 = m(6) * t130 - t203;
t300 = Ifges(6,6) * t181 + t306 * t183;
t176 = Ifges(5,4) * t183;
t299 = -Ifges(5,2) * t181 + t176;
t143 = Ifges(5,1) * t181 + t176;
t204 = t183 * mrSges(5,1) - t181 * mrSges(5,2);
t298 = t204 + t203;
t231 = t182 * t166;
t265 = t182 * pkin(3);
t266 = pkin(7) * t184;
t149 = t265 - t266;
t235 = t149 * t183;
t49 = t181 * t231 + t235;
t121 = t181 * t149;
t230 = t182 * t183;
t50 = -t166 * t230 + t121;
t297 = -t49 * t181 + t50 * t183;
t41 = -t182 * t206 + t121;
t42 = -t182 * t207 - t235;
t197 = t181 * t42 + t183 * t41;
t296 = Ifges(6,3) * t183 - t173;
t217 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t170 = m(6) * qJ(5) + mrSges(6,3);
t295 = m(5) / 0.2e1;
t294 = -m(6) / 0.2e1;
t293 = m(6) / 0.2e1;
t292 = -mrSges(5,1) / 0.2e1;
t291 = -mrSges(6,1) / 0.2e1;
t290 = mrSges(5,2) / 0.2e1;
t287 = t40 / 0.2e1;
t286 = t42 / 0.2e1;
t285 = -qJ(5) / 0.2e1;
t108 = t200 * t182;
t284 = t108 / 0.2e1;
t227 = t183 * t184;
t219 = mrSges(6,2) * t227;
t248 = t182 * mrSges(6,1);
t127 = t219 - t248;
t283 = t127 / 0.2e1;
t282 = t130 / 0.2e1;
t238 = qJ(5) * t183;
t268 = pkin(4) * t181;
t134 = -t238 + t268;
t281 = t134 / 0.2e1;
t245 = t183 * mrSges(6,3);
t251 = t181 * mrSges(6,1);
t135 = -t245 + t251;
t280 = t135 / 0.2e1;
t246 = t183 * mrSges(5,2);
t252 = t181 * mrSges(5,1);
t136 = t246 + t252;
t279 = t136 / 0.2e1;
t278 = t181 / 0.2e1;
t276 = -t183 / 0.2e1;
t275 = t183 / 0.2e1;
t274 = -t184 / 0.2e1;
t270 = m(6) * t134;
t269 = m(6) * t184;
t264 = m(6) * qJD(5);
t260 = Ifges(6,5) * t183;
t258 = Ifges(5,6) * t184;
t256 = t178 * mrSges(6,2);
t255 = t178 * mrSges(5,3);
t254 = t179 * mrSges(6,2);
t253 = t179 * mrSges(5,3);
t249 = t182 * mrSges(5,1);
t247 = t182 * mrSges(5,2);
t232 = t181 * t182;
t122 = t184 * mrSges(5,2) - mrSges(5,3) * t232;
t124 = -t184 * mrSges(5,1) - mrSges(5,3) * t230;
t125 = t184 * mrSges(6,1) + mrSges(6,2) * t230;
t171 = t184 * mrSges(6,3);
t129 = -mrSges(6,2) * t232 - t171;
t7 = t269 * t284 + (t124 * t275 + t125 * t276 + (t304 * t181 + t303 * t183) * t294 + (t253 / 0.2e1 + t255 / 0.2e1 + t256 / 0.2e1 + t254 / 0.2e1) * t182 + (t122 + t129) * t278 + t298 * t184 / 0.2e1) * t182;
t241 = t7 * qJD(1);
t240 = -t204 - mrSges(4,1);
t109 = t135 * t182;
t193 = t134 + t166;
t51 = t193 * t182;
t15 = m(6) * (-t184 * t34 - t230 * t51) - t109 * t230 - t184 * t129;
t237 = qJD(1) * t15;
t236 = t130 * t182;
t228 = t183 * t203;
t224 = t309 * t266;
t222 = qJD(3) * t184;
t221 = qJD(4) * t182;
t220 = m(6) * t286;
t163 = Ifges(6,5) * t230;
t83 = -Ifges(6,6) * t184 + Ifges(6,3) * t232 + t163;
t85 = t182 * t299 - t258;
t216 = t83 / 0.2e1 - t85 / 0.2e1;
t212 = t122 / 0.2e1 + t129 / 0.2e1;
t211 = -t124 / 0.2e1 + t125 / 0.2e1;
t139 = Ifges(5,2) * t183 + t261;
t210 = -t296 / 0.2e1 - t139 / 0.2e1;
t141 = Ifges(6,1) * t181 - t260;
t209 = t141 / 0.2e1 + t143 / 0.2e1;
t138 = Ifges(6,3) * t181 + t260;
t110 = t184 * t135;
t111 = t184 * t136;
t123 = -mrSges(5,3) * t225 - t247;
t126 = -mrSges(5,3) * t227 + t249;
t128 = -mrSges(6,2) * t225 + t182 * mrSges(6,3);
t52 = t193 * t184;
t84 = Ifges(6,6) * t182 + t138 * t184;
t86 = Ifges(5,6) * t182 + t184 * t299;
t88 = Ifges(6,4) * t182 + t184 * t301;
t90 = Ifges(5,5) * t182 + t144 * t184;
t3 = t52 * t109 + t51 * t110 + t50 * t122 + t40 * t123 + t49 * t124 + t42 * t125 + t39 * t126 + t35 * t127 + t34 * t128 + t41 * t129 + m(6) * (t34 * t41 + t35 * t42 + t51 * t52) + m(5) * (t39 * t49 + t40 * t50) + (Ifges(4,4) * t184 + t214 * mrSges(4,2) + t308 * t183 + (-t184 * t217 + t216) * t181) * t184 + (t166 * t111 + t214 * mrSges(4,1) - Ifges(4,4) * t182 + (t88 / 0.2e1 + t90 / 0.2e1 + t218 * t182) * t183 + (t84 / 0.2e1 - t86 / 0.2e1 + t217 * t182) * t181 + (Ifges(4,1) - Ifges(4,2) + (m(5) * t166 + t136) * t166 - t305) * t184) * t182;
t192 = m(5) * t297;
t6 = (-t111 / 0.2e1 - t110 / 0.2e1 + t212 * t183 + t211 * t181 + (-t181 * t39 + t183 * t40) * t295 + (t181 * t35 + t183 * t34 - t52) * t293 - m(5) * t226 / 0.2e1) * t184 + (t109 / 0.2e1 + (t123 / 0.2e1 + t128 / 0.2e1 + t247 / 0.2e1) * t183 + (-t126 / 0.2e1 + t283 + t249 / 0.2e1) * t181 + t192 / 0.2e1 + (t197 + t51) * t293 + t231 * t295) * t182;
t199 = t3 * qJD(1) + t6 * qJD(2);
t112 = t296 * t182;
t113 = t182 * t139;
t114 = -Ifges(6,1) * t232 + t163;
t115 = t182 * t143;
t162 = Ifges(6,6) * t230;
t5 = m(6) * (t108 * t51 + t34 * t39 + t35 * t40) + t108 * t109 + t162 * t274 + t39 * t122 - t40 * t124 + t40 * t125 + t39 * t129 + ((t51 * mrSges(6,1) + t258 / 0.2e1 - t40 * mrSges(5,3) - t34 * mrSges(6,2) + t114 / 0.2e1 - t115 / 0.2e1 + mrSges(5,1) * t231 + t216) * t183 + (t51 * mrSges(6,3) + t39 * mrSges(5,3) - t35 * mrSges(6,2) + t112 / 0.2e1 + t113 / 0.2e1 - mrSges(5,2) * t231 - t308) * t181) * t182;
t198 = t5 * qJD(1) - t7 * qJD(2);
t23 = 0.4e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * (-0.1e1 + t309) * t184 * t182;
t196 = t6 * qJD(1) + t23 * qJD(2);
t190 = (-t181 * t51 + (-t236 - t266) * t183) * t293 - t181 * t109 / 0.2e1;
t13 = t219 + (t291 - t228 / 0.2e1) * t182 + t220 - t190;
t45 = t302 * t181;
t195 = qJD(1) * t13 + qJD(3) * t45;
t19 = -t171 + 0.2e1 * (t233 / 0.4e1 - t40 / 0.4e1 + (t234 / 0.4e1 + t285) * t184) * m(6);
t194 = qJD(1) * t19 + qJD(4) * t170;
t10 = -pkin(3) * t136 - t134 * t203 + (t135 + t270) * t130 + (-t138 / 0.2e1 + t299 / 0.2e1 + t209) * t183 + (t301 / 0.2e1 + t144 / 0.2e1 + t210) * t181;
t16 = (t279 + t280 + (-mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1) * t183 + (t292 + t291) * t181 + (-t268 / 0.2e1 + t238 / 0.2e1 + t281) * m(6)) * t184;
t185 = t166 * t279 + (t144 / 0.4e1 + t301 / 0.4e1 - t139 / 0.4e1 - t296 / 0.4e1 + pkin(3) * t292 + mrSges(6,1) * t282) * t183 + (-t143 / 0.4e1 - t141 / 0.4e1 - t299 / 0.4e1 + t138 / 0.4e1 + pkin(3) * t290 + mrSges(6,3) * t282) * t181 + (mrSges(6,2) + mrSges(5,3)) * pkin(7) * (-t179 / 0.2e1 - t178 / 0.2e1);
t186 = (t108 * t130 + t134 * t51) * t293 - t203 * t284 + t109 * t281 + t51 * t280 - t300 * t184 / 0.4e1;
t187 = (-pkin(4) * t42 + qJ(5) * t41) * t294 + pkin(4) * t283 + t128 * t285 - t41 * mrSges(6,3) / 0.2e1 + mrSges(6,1) * t286 + t49 * t292 + t50 * t290;
t188 = (t287 - t34 / 0.2e1) * mrSges(6,2) + (t304 * t293 - t212) * pkin(7) + t114 / 0.4e1 - t115 / 0.4e1 + t83 / 0.4e1 - t85 / 0.4e1;
t189 = -t113 / 0.4e1 - t112 / 0.4e1 + t89 / 0.4e1 + t87 / 0.4e1 + (t39 / 0.2e1 + t35 / 0.2e1) * mrSges(6,2) + (t303 * t293 + t211) * pkin(7);
t2 = (-t205 + t189) * t183 + ((0.3e1 / 0.4e1 * Ifges(5,6) - Ifges(6,6) / 0.2e1) * t184 + t188) * t181 + (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1 + t185) * t182 + t186 + t187;
t191 = t2 * qJD(1) - t16 * qJD(2) + t10 * qJD(3);
t131 = (m(6) * pkin(7) + mrSges(6,2)) * t183;
t18 = (t233 + (-0.2e1 * qJ(5) + t234) * t184) * t293 + m(6) * t287 + t129;
t17 = -t134 * t269 / 0.2e1 + (-t246 / 0.2e1 - t252 / 0.2e1 - t251 / 0.2e1 + t245 / 0.2e1 - t270 / 0.2e1) * t184 + (t136 + t135) * t274;
t14 = t228 * t277 + t220 - t248 / 0.2e1 + t190;
t4 = qJD(3) * t6 - qJD(4) * t7;
t1 = (t258 / 0.4e1 + t188) * t181 + t189 * t183 + t186 - t187 + t185 * t182 + t305 * t277 + t217 * t225 + t306 * t227 / 0.2e1;
t8 = [qJD(3) * t3 + qJD(4) * t5 + qJD(5) * t15, t4, t1 * qJD(4) + t14 * qJD(5) + (Ifges(4,5) + t209 * t183 + t210 * t181 + (-m(5) * pkin(3) + t240) * t166) * t222 + t199 + (mrSges(4,2) * t231 - pkin(3) * t111 + t130 * t110 - Ifges(4,6) * t182 + t84 * t276 + t86 * t275 + ((t123 + t128) * t183 + (-t126 + t127) * t181 + m(6) * t197 + t192) * pkin(7) + t302 * t52 + (t88 + t90) * t278 + ((Ifges(5,6) - Ifges(6,6)) * t183 + t306 * t181) * t277 + t297 * mrSges(5,3) + t197 * mrSges(6,2)) * qJD(3), t1 * qJD(3) + (t162 + (-m(6) * pkin(4) + t307) * t40 + (-mrSges(5,2) + t170) * t39) * qJD(4) + t18 * qJD(5) + ((-qJ(5) * mrSges(6,2) - Ifges(5,6)) * t183 + (pkin(4) * mrSges(6,2) - t306) * t181) * t221 + t198, qJD(3) * t14 + qJD(4) * t18 + t237; t4, t23 * qJD(3), t17 * qJD(4) + (-t203 + t240) * qJD(3) * t182 + 0.2e1 * ((t224 - t265) * t295 + (t224 + t236) * t293) * qJD(3) + ((-mrSges(4,2) + t253 + t254 + t255 + t256) * qJD(3) + t181 * t264) * t184 + t196, -t241 + t17 * qJD(3) - m(6) * t108 * qJD(4) + ((mrSges(5,2) - mrSges(6,3)) * qJD(4) * t181 + (t307 * qJD(4) + t264) * t183) * t182, (t181 * t222 + t183 * t221) * m(6); qJD(4) * t2 - qJD(5) * t13 - t199, -qJD(4) * t16 - t196, qJD(4) * t10 - qJD(5) * t45, t131 * qJD(5) + t191 + (-Ifges(5,6) * t181 + (-m(6) * t200 - t298) * pkin(7) - t200 * mrSges(6,2) + t300) * qJD(4), qJD(4) * t131 - t195; -qJD(3) * t2 + qJD(5) * t19 - t198, qJD(3) * t16 + t241, -t191, t170 * qJD(5), t194; qJD(3) * t13 - qJD(4) * t19 - t237, 0, t195, -t194, 0;];
Cq = t8;
