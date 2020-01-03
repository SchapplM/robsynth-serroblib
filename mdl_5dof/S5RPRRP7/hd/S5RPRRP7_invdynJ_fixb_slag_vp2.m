% Calculate vector of inverse dynamics joint torques for
% S5RPRRP7
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:33
% EndTime: 2019-12-31 18:44:53
% DurationCPUTime: 11.18s
% Computational Cost: add. (2819->465), mult. (5876->602), div. (0->0), fcn. (3383->10), ass. (0->207)
t336 = Ifges(5,1) + Ifges(6,1);
t335 = Ifges(6,4) + Ifges(5,5);
t333 = Ifges(5,6) - Ifges(6,6);
t135 = sin(qJ(4));
t138 = cos(qJ(4));
t219 = t138 * qJD(3);
t136 = sin(qJ(3));
t228 = qJD(1) * t136;
t94 = t135 * t228 - t219;
t280 = t94 / 0.2e1;
t225 = qJD(3) * t135;
t95 = t138 * t228 + t225;
t339 = -t95 / 0.2e1;
t139 = cos(qJ(3));
t227 = qJD(1) * t139;
t115 = qJD(4) - t227;
t338 = -t115 / 0.2e1;
t337 = -mrSges(5,3) - mrSges(6,2);
t334 = Ifges(6,2) + Ifges(5,3);
t304 = -t135 * t333 + t138 * t335;
t248 = Ifges(6,5) * t135;
t250 = Ifges(5,4) * t135;
t302 = t138 * t336 + t248 - t250;
t178 = t138 * mrSges(6,1) + t135 * mrSges(6,3);
t180 = mrSges(5,1) * t138 - mrSges(5,2) * t135;
t332 = -t178 - t180;
t218 = qJD(1) * qJD(3);
t100 = qJDD(1) * t136 + t139 * t218;
t237 = qJD(4) * t94;
t48 = qJDD(3) * t135 + t100 * t138 - t237;
t285 = t48 / 0.2e1;
t49 = qJD(4) * t95 - t138 * qJDD(3) + t100 * t135;
t284 = -t49 / 0.2e1;
t99 = qJDD(1) * t139 - t136 * t218;
t90 = qJDD(4) - t99;
t282 = t90 / 0.2e1;
t331 = t99 / 0.2e1;
t330 = t100 / 0.2e1;
t329 = t218 / 0.2e1;
t133 = sin(pkin(8));
t119 = pkin(1) * t133 + pkin(6);
t328 = qJD(2) * qJD(3) + t119 * qJDD(1);
t220 = qJD(4) * t138;
t222 = qJD(4) * t135;
t104 = t119 * qJD(1);
t224 = qJD(3) * t136;
t32 = t136 * qJDD(2) - t104 * t224 + t139 * t328;
t29 = qJDD(3) * pkin(7) + t32;
t134 = cos(pkin(8));
t267 = pkin(1) * t134;
t120 = -pkin(2) - t267;
t103 = t120 * qJDD(1);
t47 = -pkin(3) * t99 - pkin(7) * t100 + t103;
t226 = qJD(2) * t136;
t74 = t104 * t139 + t226;
t66 = qJD(3) * pkin(7) + t74;
t229 = t139 * pkin(3) + t136 * pkin(7);
t158 = -pkin(2) - t229;
t67 = (t158 - t267) * qJD(1);
t3 = t135 * t47 + t138 * t29 + t67 * t220 - t222 * t66;
t25 = t135 * t67 + t138 * t66;
t4 = -qJD(4) * t25 - t135 * t29 + t138 * t47;
t183 = -t135 * t4 + t138 * t3;
t24 = -t135 * t66 + t138 * t67;
t327 = -t24 * t220 - t25 * t222 + t183;
t16 = -pkin(4) * t115 + qJD(5) - t24;
t17 = qJ(5) * t115 + t25;
t1 = qJ(5) * t90 + qJD(5) * t115 + t3;
t2 = -pkin(4) * t90 + qJDD(5) - t4;
t184 = t1 * t138 + t135 * t2;
t326 = t16 * t220 - t17 * t222 + t184;
t18 = mrSges(5,1) * t90 - mrSges(5,3) * t48;
t19 = -t90 * mrSges(6,1) + t48 * mrSges(6,2);
t20 = -mrSges(5,2) * t90 - mrSges(5,3) * t49;
t21 = -mrSges(6,2) * t49 + mrSges(6,3) * t90;
t325 = (t20 + t21) * t138 + (-t18 + t19) * t135;
t123 = Ifges(4,4) * t227;
t268 = Ifges(6,5) * t94;
t88 = Ifges(5,4) * t94;
t313 = t335 * t115 + t336 * t95 + t268 - t88;
t87 = Ifges(6,5) * t95;
t34 = t115 * Ifges(6,6) + t94 * Ifges(6,3) + t87;
t324 = Ifges(4,1) * t228 + Ifges(4,5) * qJD(3) + t135 * t34 + t138 * t313 + t123;
t252 = Ifges(4,4) * t136;
t172 = t139 * Ifges(4,2) + t252;
t323 = Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t172 / 0.2e1 + t335 * t339 + t333 * t280 + t334 * t338;
t322 = -t4 * mrSges(5,1) + t2 * mrSges(6,1) + t3 * mrSges(5,2) - t1 * mrSges(6,3);
t320 = -m(5) - m(6);
t319 = t335 * t90 + (-Ifges(5,4) + Ifges(6,5)) * t49 + t336 * t48;
t258 = -qJD(1) / 0.2e1;
t316 = -mrSges(4,3) + mrSges(3,2);
t271 = mrSges(5,3) * t94;
t60 = -mrSges(5,2) * t115 - t271;
t273 = mrSges(6,2) * t94;
t63 = mrSges(6,3) * t115 - t273;
t256 = t60 + t63;
t270 = mrSges(5,3) * t95;
t61 = mrSges(5,1) * t115 - t270;
t272 = mrSges(6,2) * t95;
t62 = -mrSges(6,1) * t115 + t272;
t255 = -t61 + t62;
t11 = mrSges(5,1) * t49 + mrSges(5,2) * t48;
t312 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t100 + t11;
t161 = pkin(4) * t135 - qJ(5) * t138;
t311 = qJD(4) * t161 - qJD(5) * t135 - t226 - (qJD(1) * t161 + t104) * t139;
t212 = mrSges(4,3) * t228;
t310 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t94 + mrSges(5,2) * t95 + t212;
t182 = mrSges(4,1) * t139 - mrSges(4,2) * t136;
t309 = -t182 - mrSges(3,1);
t307 = t334 * t136 + t304 * t139;
t306 = t335 * t136 + t302 * t139;
t177 = t135 * mrSges(6,1) - t138 * mrSges(6,3);
t179 = mrSges(5,1) * t135 + mrSges(5,2) * t138;
t73 = qJD(2) * t139 - t136 * t104;
t65 = -qJD(3) * pkin(3) - t73;
t23 = pkin(4) * t94 - qJ(5) * t95 + t65;
t305 = t23 * t177 + t65 * t179;
t303 = t335 * t135 + t333 * t138;
t247 = Ifges(6,5) * t138;
t249 = Ifges(5,4) * t138;
t301 = t336 * t135 - t247 + t249;
t298 = -t333 * t49 + t334 * t90 + t335 * t48;
t223 = qJD(3) * t139;
t33 = qJDD(2) * t139 - t104 * t223 - t136 * t328;
t297 = -t136 * t33 + t139 * t32;
t296 = t337 * t136;
t295 = -m(4) + t320;
t294 = -t48 * Ifges(6,5) / 0.2e1 - t90 * Ifges(6,6) / 0.2e1 + Ifges(5,4) * t285 + Ifges(5,6) * t282 + (Ifges(6,3) + Ifges(5,2)) * t284;
t162 = pkin(4) * t138 + qJ(5) * t135;
t101 = -pkin(3) - t162;
t181 = mrSges(4,1) * t136 + mrSges(4,2) * t139;
t293 = t181 + t337 * t139 + (m(5) * pkin(3) - m(6) * t101 - t332) * t136;
t230 = t138 * t139;
t86 = t120 - t229;
t253 = t119 * t230 + t135 * t86;
t263 = pkin(7) * t139;
t185 = pkin(3) * t136 - t263;
t98 = t185 * qJD(3);
t289 = -qJD(4) * t253 + t138 * t98;
t288 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t287 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t269 = Ifges(5,4) * t95;
t37 = -t94 * Ifges(5,2) + t115 * Ifges(5,6) + t269;
t286 = -t37 / 0.2e1;
t283 = t49 / 0.2e1;
t281 = -t94 / 0.2e1;
t278 = t95 / 0.2e1;
t275 = t135 / 0.2e1;
t137 = sin(qJ(1));
t266 = pkin(1) * t137;
t140 = cos(qJ(1));
t131 = t140 * pkin(1);
t257 = qJD(4) / 0.2e1;
t97 = t185 * qJD(1);
t41 = t135 * t97 + t138 * t73;
t254 = t135 * t98 + t86 * t220;
t251 = Ifges(4,4) * t139;
t242 = t138 * t86;
t132 = qJ(1) + pkin(8);
t125 = cos(t132);
t235 = t125 * t136;
t234 = t125 * t139;
t233 = t135 * t136;
t232 = t135 * t139;
t221 = qJD(4) * t136;
t211 = mrSges(4,3) * t227;
t124 = sin(t132);
t206 = t125 * pkin(2) + t124 * pkin(6) + t131;
t205 = t119 * t224;
t203 = t135 * t223;
t192 = t220 / 0.2e1;
t190 = t119 * t135 + pkin(4);
t171 = -Ifges(5,2) * t135 + t249;
t170 = Ifges(5,2) * t138 + t250;
t167 = Ifges(4,5) * t139 - Ifges(4,6) * t136;
t164 = Ifges(6,3) * t135 + t247;
t163 = -Ifges(6,3) * t138 + t248;
t40 = -t135 * t73 + t138 * t97;
t157 = t119 + t161;
t153 = t120 * qJD(1) * t181;
t152 = t136 * (Ifges(4,1) * t139 - t252);
t30 = -qJDD(3) * pkin(3) - t33;
t146 = Ifges(5,6) * t136 + t139 * t171;
t145 = Ifges(6,6) * t136 + t139 * t164;
t107 = -qJD(3) * mrSges(4,2) + t211;
t85 = t179 * t136;
t76 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t99;
t71 = t124 * t135 + t125 * t230;
t70 = -t124 * t138 + t125 * t232;
t69 = t124 * t230 - t125 * t135;
t68 = t124 * t232 + t125 * t138;
t59 = t157 * t136;
t56 = mrSges(6,1) * t94 - mrSges(6,3) * t95;
t55 = pkin(4) * t95 + qJ(5) * t94;
t53 = -t119 * t232 + t242;
t52 = t139 * t190 - t242;
t51 = -qJ(5) * t139 + t253;
t28 = -pkin(4) * t228 - t40;
t27 = qJ(5) * t228 + t41;
t22 = (qJD(4) * t162 - qJD(5) * t138) * t136 + t157 * t223;
t15 = t135 * t205 + t289;
t14 = (-t136 * t219 - t139 * t222) * t119 + t254;
t13 = -t190 * t224 - t289;
t12 = (-t119 * t222 - qJD(5)) * t139 + (-t119 * t138 + qJ(5)) * t224 + t254;
t10 = mrSges(6,1) * t49 - mrSges(6,3) * t48;
t5 = pkin(4) * t49 - qJ(5) * t48 - qJD(5) * t95 + t30;
t6 = [(-t1 * mrSges(6,2) - t3 * mrSges(5,3) - t294) * t233 - t107 * t205 + m(5) * (t14 * t25 + t15 * t24 + t3 * t253 + t4 * t53) + (t139 * t76 + m(5) * (t136 * t30 + t223 * t65) + m(4) * ((-t136 * t74 - t139 * t73) * qJD(3) + t297) + t310 * t223 + t312 * t136) * t119 + t120 * (-mrSges(4,1) * t99 + mrSges(4,2) * t100) + t324 * t223 / 0.2e1 + t253 * t20 + t30 * t85 + t22 * t56 + t59 * t10 + t14 * t60 + t15 * t61 + t13 * t62 + t12 * t63 + t51 * t21 + t52 * t19 + t53 * t18 + (-m(3) * t131 - m(4) * t206 - mrSges(2,1) * t140 + mrSges(2,2) * t137 + t320 * (pkin(3) * t234 + pkin(7) * t235 + t206) - t288 * t71 - t287 * t70 + t337 * t235 + t309 * t125 + t316 * t124) * g(2) + (Ifges(4,6) * qJDD(3) + (-Ifges(4,2) * t136 + t251) * t329 + Ifges(4,4) * t330 + Ifges(4,2) * t331 - Ifges(6,6) * t283 - Ifges(5,6) * t284 - t334 * t282 - t298 / 0.2e1 - t335 * t285 + t322) * t139 + (-t4 * mrSges(5,3) + t2 * mrSges(6,2) + t319 / 0.2e1) * t136 * t138 + t203 * t286 + (qJD(3) * t306 - t221 * t301) * t278 + (qJD(3) * t307 - t221 * t303) * t115 / 0.2e1 + qJD(3) * t153 + (t24 * mrSges(5,1) - t16 * mrSges(6,1) - t25 * mrSges(5,2) - t74 * mrSges(4,3) + t17 * mrSges(6,3) - t323) * t224 + (qJD(3) * t145 - t163 * t221) * t280 + (qJD(3) * t146 - t170 * t221) * t281 + m(6) * (t1 * t51 + t12 * t17 + t13 * t16 + t2 * t52 + t22 * t23 + t5 * t59) + qJD(3) ^ 2 * t167 / 0.2e1 + (Ifges(4,1) * t100 + Ifges(4,5) * qJDD(3) + t164 * t283 + t171 * t284 + t5 * t177 + t34 * t192 + t304 * t282 + t302 * t285) * t136 + (m(3) * t266 + mrSges(2,1) * t137 + mrSges(2,2) * t140 + t288 * t69 + t287 * t68 + t295 * (t125 * pkin(6) - t266) + t316 * t125 + (m(4) * pkin(2) + t158 * t320 - t296 - t309) * t124) * g(1) + t152 * t329 + t251 * t330 + (t252 + t172) * t331 + (-t223 * t73 + t297) * mrSges(4,3) + (mrSges(5,1) * t65 + mrSges(6,1) * t23 - mrSges(6,2) * t17 - mrSges(5,3) * t25) * (t136 * t220 + t203) + (mrSges(5,2) * t65 + mrSges(6,2) * t16 - mrSges(5,3) * t24 - mrSges(6,3) * t23) * (-t135 * t221 + t139 * t219) - (t135 * t313 + t138 * t37) * t221 / 0.2e1 + (m(4) * t120 - t182) * t103 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t134 - 0.2e1 * mrSges(3,2) * t133 + m(3) * (t133 ^ 2 + t134 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1); m(3) * qJDD(2) + (-m(3) + t295) * g(3) + (-t10 + (t135 * t255 + t138 * t256 + t107) * qJD(3) + m(5) * (t219 * t25 - t225 * t24 - t30) + m(6) * (t16 * t225 + t17 * t219 - t5) + m(4) * (qJD(3) * t74 + t33) - t312) * t139 + (t76 + (t56 + t310) * qJD(3) + (-t135 * t256 + t138 * t255) * qJD(4) + m(5) * (qJD(3) * t65 + t327) + m(6) * (qJD(3) * t23 + t326) + m(4) * (-qJD(3) * t73 + t32) + t325) * t136; -t167 * t218 / 0.2e1 + (t257 * t302 + t258 * t306) * t95 + ((t146 / 0.2e1 - t145 / 0.2e1) * t94 - t24 * (mrSges(5,1) * t136 - mrSges(5,3) * t230) - t16 * (-mrSges(6,1) * t136 + mrSges(6,2) * t230) - t25 * (-mrSges(5,2) * t136 - mrSges(5,3) * t232) - t17 * (-mrSges(6,2) * t232 + mrSges(6,3) * t136) - t153 + t152 * t258) * qJD(1) + t294 * t138 + (-t182 + t320 * t229 + (-m(6) * t162 + t332) * t139 + t296) * g(3) + (t263 * t320 + t293) * g(2) * t124 + (t257 * t304 + t258 * t307) * t115 + t293 * t125 * g(1) - (-Ifges(4,2) * t228 + t123 + t324) * t227 / 0.2e1 + (m(5) * ((-t25 * t135 - t24 * t138) * qJD(4) + t183) + m(6) * ((-t17 * t135 + t16 * t138) * qJD(4) + t184) + t255 * t220 - t256 * t222 + t234 * t320 * g(1) + t325) * pkin(7) + t326 * mrSges(6,2) + t327 * mrSges(5,3) + (-t107 + t211) * t73 + Ifges(4,6) * t99 + Ifges(4,5) * t100 + t101 * t10 - t41 * t60 - t40 * t61 - t28 * t62 - t27 * t63 - t32 * mrSges(4,2) + t33 * mrSges(4,1) - pkin(3) * t11 + (-m(5) * t65 + t212 - t310) * t74 + t311 * t56 + t313 * t192 + t163 * t283 + t170 * t284 + (-t171 / 0.2e1 + t164 / 0.2e1) * t237 + (t34 / 0.2e1 + t286) * t222 + t301 * t285 + t303 * t282 + (t275 * t37 - t305) * t227 + t305 * qJD(4) + t323 * t228 - t5 * t178 - t30 * t180 + t319 * t275 + (-pkin(3) * t30 - t24 * t40 - t25 * t41) * m(5) + (t101 * t5 - t16 * t28 - t17 * t27 + t311 * t23) * m(6) + Ifges(4,3) * qJDD(3); (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t17 - t23 * t55) * m(6) + t298 + (-t336 * t94 - t269 + t34 + t87) * t339 + (t85 - (-m(6) * t161 - t177) * t136) * g(3) - t23 * (mrSges(6,1) * t95 + mrSges(6,3) * t94) - t65 * (mrSges(5,1) * t95 - mrSges(5,2) * t94) - t55 * t56 + qJD(5) * t63 - pkin(4) * t19 + qJ(5) * t21 - t322 + (-m(6) * t16 - t255 + t270) * t25 + (-m(6) * t17 - t256 - t271) * t24 + (-Ifges(5,2) * t95 + t313 - t88) * t280 + (-t287 * t71 + t288 * t70) * g(1) + (-t287 * t69 + t288 * t68) * g(2) + t17 * t272 + t16 * t273 + t37 * t278 + (Ifges(6,3) * t95 - t268) * t281 + (-t333 * t95 - t335 * t94) * t338; -t115 * t63 + t95 * t56 + (-g(1) * t70 - g(2) * t68 - g(3) * t233 - t17 * t115 + t23 * t95 + t2) * m(6) + t19;];
tau = t6;
