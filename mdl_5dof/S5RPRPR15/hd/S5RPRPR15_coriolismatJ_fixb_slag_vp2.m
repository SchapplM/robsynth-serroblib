% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR15_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR15_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:40
% EndTime: 2019-12-31 18:36:45
% DurationCPUTime: 2.06s
% Computational Cost: add. (5072->328), mult. (10610->469), div. (0->0), fcn. (10280->6), ass. (0->175)
t171 = sin(pkin(8));
t172 = cos(pkin(8));
t254 = sin(qJ(5));
t255 = cos(qJ(5));
t143 = -t254 * t171 + t255 * t172;
t138 = Ifges(6,4) * t143;
t185 = t255 * t171 + t254 * t172;
t94 = Ifges(6,1) * t185 + t138;
t282 = -Ifges(6,2) * t185 + t138 + t94;
t173 = sin(qJ(3));
t174 = cos(qJ(3));
t281 = t173 * t174;
t235 = t172 * mrSges(5,2);
t238 = t171 * mrSges(5,1);
t280 = t173 * (t235 + t238);
t155 = pkin(3) * t174 + qJ(4) * t173;
t146 = t172 * t155;
t175 = -pkin(1) - pkin(6);
t217 = t174 * t175;
t104 = -t171 * t217 + t146;
t105 = t171 * t155 + t172 * t217;
t194 = -t104 * t171 + t105 * t172;
t253 = pkin(3) * t173;
t151 = -qJ(4) * t174 + qJ(2) + t253;
t140 = t172 * t151;
t218 = t173 * t175;
t102 = -t171 * t218 + t140;
t103 = t171 * t151 + t172 * t218;
t222 = t171 * t174;
t148 = -t173 * mrSges(5,2) - mrSges(5,3) * t222;
t219 = t172 * t174;
t150 = t173 * mrSges(5,1) - mrSges(5,3) * t219;
t279 = -m(5) * (t102 * t172 + t103 * t171) - t171 * t148 - t172 * t150;
t278 = -m(5) / 0.2e1;
t277 = m(5) / 0.2e1;
t276 = -m(6) / 0.2e1;
t275 = m(6) / 0.2e1;
t125 = t185 * t174;
t99 = -t173 * mrSges(6,2) - t125 * mrSges(6,3);
t274 = -t99 / 0.2e1;
t124 = t185 * t173;
t126 = t143 * t173;
t127 = t143 * t174;
t273 = m(6) * (t124 * t127 - t125 * t126);
t101 = t173 * mrSges(6,1) - t127 * mrSges(6,3);
t272 = -t101 / 0.2e1;
t271 = t124 / 0.2e1;
t270 = -t125 / 0.2e1;
t268 = -t126 / 0.2e1;
t267 = t127 / 0.2e1;
t265 = t143 / 0.2e1;
t263 = t185 / 0.2e1;
t262 = -t185 / 0.2e1;
t261 = t171 / 0.2e1;
t260 = -t172 / 0.2e1;
t259 = t172 / 0.2e1;
t258 = t173 / 0.2e1;
t257 = -t174 / 0.2e1;
t256 = t174 / 0.2e1;
t252 = pkin(7) + qJ(4);
t251 = Ifges(5,4) * t171;
t250 = Ifges(5,4) * t172;
t249 = Ifges(6,4) * t185;
t100 = mrSges(6,1) * t174 + mrSges(6,3) * t126;
t237 = t171 * Ifges(5,2);
t122 = Ifges(5,6) * t174 + (t237 - t250) * t173;
t234 = t172 * Ifges(5,1);
t123 = Ifges(5,5) * t174 + (-t234 + t251) * t173;
t205 = pkin(4) * t171 - t175;
t141 = t205 * t173;
t142 = t205 * t174;
t223 = t171 * t173;
t147 = -mrSges(5,2) * t174 + mrSges(5,3) * t223;
t220 = t172 * t173;
t149 = mrSges(5,1) * t174 + mrSges(5,3) * t220;
t192 = Ifges(6,5) * t268 + Ifges(6,6) * t271;
t243 = t127 * mrSges(6,2);
t246 = t125 * mrSges(6,1);
t200 = t243 + t246;
t233 = t172 * Ifges(5,5);
t236 = t171 * Ifges(5,6);
t206 = -t171 * t175 + pkin(4);
t84 = -pkin(7) * t219 + t206 * t173 + t140;
t90 = -pkin(7) * t222 + t103;
t35 = -t254 * t90 + t255 * t84;
t36 = t254 * t84 + t255 * t90;
t89 = pkin(7) * t220 + t206 * t174 + t146;
t95 = pkin(7) * t223 + t105;
t37 = -t254 * t95 + t255 * t89;
t38 = t254 * t89 + t255 * t95;
t68 = -Ifges(6,4) * t126 + Ifges(6,2) * t124 + Ifges(6,6) * t174;
t242 = t127 * Ifges(6,4);
t69 = -t125 * Ifges(6,2) + t173 * Ifges(6,6) + t242;
t70 = -Ifges(6,1) * t126 + Ifges(6,4) * t124 + Ifges(6,5) * t174;
t115 = Ifges(6,4) * t125;
t71 = t127 * Ifges(6,1) + t173 * Ifges(6,5) - t115;
t244 = t126 * mrSges(6,2);
t247 = t124 * mrSges(6,1);
t74 = -t244 - t247;
t98 = -mrSges(6,2) * t174 + mrSges(6,3) * t124;
t1 = t103 * t147 + t105 * t148 + t102 * t149 + t104 * t150 - t141 * t200 + t142 * t74 + t71 * t268 + t70 * t267 + t68 * t270 + t69 * t271 + t36 * t98 + t38 * t99 + t35 * t100 + t37 * t101 + m(6) * (-t141 * t142 + t35 * t37 + t36 * t38) + m(5) * (t102 * t104 + t103 * t105) + (qJ(2) * mrSges(4,1) + Ifges(6,5) * t267 + Ifges(6,6) * t270 - t171 * t122 / 0.2e1 + t123 * t259 + t175 * t280 + (-Ifges(4,4) + t233 / 0.2e1 - t236 / 0.2e1) * t174) * t174 + (-qJ(2) * mrSges(4,2) + (Ifges(4,4) - t233 + t236) * t173 + (-Ifges(4,1) + Ifges(4,2) + Ifges(5,3) + Ifges(6,3) - m(5) * t175 ^ 2 + (t175 * mrSges(5,2) - t234 / 0.2e1) * t172 + (t175 * mrSges(5,1) + t250 - t237 / 0.2e1) * t171) * t174 + t192) * t173;
t248 = t1 * qJD(1);
t152 = t252 * t171;
t154 = t252 * t172;
t96 = -t255 * t152 - t254 * t154;
t245 = t125 * t96;
t97 = -t254 * t152 + t255 * t154;
t241 = t127 * t97;
t240 = t143 * mrSges(6,1);
t239 = t185 * mrSges(6,2);
t216 = -Ifges(6,5) * t125 - Ifges(6,6) * t127;
t73 = t127 * mrSges(6,1) - mrSges(6,2) * t125;
t75 = -Ifges(6,2) * t127 - t115;
t76 = -Ifges(6,1) * t125 - t242;
t4 = t35 * t99 - t36 * t101 + t216 * t258 + t142 * t73 + (-t36 * mrSges(6,3) + t76 / 0.2e1 - t69 / 0.2e1) * t127 - (-t35 * mrSges(6,3) + t71 / 0.2e1 + t75 / 0.2e1) * t125;
t232 = t4 * qJD(1);
t226 = t126 * t127;
t227 = t124 * t125;
t179 = (-t226 / 0.2e1 - t227 / 0.2e1) * mrSges(6,3) + t124 * t274 + t101 * t268 + t73 * t257;
t190 = -t240 / 0.2e1 + t239 / 0.2e1;
t7 = t179 + t190;
t231 = t7 * qJD(1);
t230 = -mrSges(5,1) * t172 + mrSges(5,2) * t171 - mrSges(4,1);
t16 = t173 * mrSges(4,1) + t174 * mrSges(4,2) + t143 * t101 + t185 * t99 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(6) * (t143 * t35 + t185 * t36) - t279;
t225 = t16 * qJD(1);
t224 = t171 * t149;
t221 = t172 * t147;
t215 = Ifges(6,5) * t143 - Ifges(6,6) * t185;
t169 = t171 ^ 2;
t170 = t172 ^ 2;
t214 = t169 + t170;
t213 = qJD(3) * t173;
t212 = t273 / 0.2e1;
t211 = m(5) * t258;
t91 = mrSges(6,1) * t185 + mrSges(6,2) * t143;
t210 = t91 * t257;
t207 = -t223 / 0.2e1;
t204 = t214 * mrSges(5,3);
t202 = t214 * qJ(4);
t199 = Ifges(6,1) * t143 - t249;
t24 = m(5) * (-0.1e1 + t214) * t281 + (t226 + t227 - t281) * m(6);
t195 = -t102 * t171 + t103 * t172;
t177 = (t195 * t174 + (t194 - 0.2e1 * t217) * t173) * t278 + (-t124 * t37 - t35 * t125 + t126 * t38 + t36 * t127 + t141 * t174 + t142 * t173) * t276 + t100 * t271 - t125 * t272 + t98 * t268 + t127 * t274;
t186 = (t143 * t96 + t185 * t97) * t275;
t189 = t148 * t260 + t150 * t261;
t191 = -t246 / 0.2e1 - t243 / 0.2e1;
t5 = (t224 / 0.2e1 - t221 / 0.2e1 + t191) * t173 + (t74 / 0.2e1 - t280 + t189) * t174 + t186 + t177;
t198 = -t5 * qJD(1) + t24 * qJD(2);
t197 = qJD(1) * t73 + qJD(3) * t91;
t14 = -t125 * t99 - t127 * t101 + m(6) * (-t125 * t36 - t127 * t35) + t279 * t174;
t187 = t14 * qJD(1) + qJD(2) * t212;
t184 = (t124 * t185 + t126 * t143) * t275;
t165 = -pkin(4) * t172 - pkin(3);
t93 = t143 * Ifges(6,2) + t249;
t11 = -t165 * t91 + t199 * t262 + t93 * t263 - t282 * t143 / 0.2e1;
t17 = t210 - t191;
t176 = (-t241 / 0.2e1 + t245 / 0.2e1) * mrSges(6,3) + t142 * t91 / 0.2e1 + t165 * t73 / 0.2e1 + t173 * t215 / 0.4e1 + t96 * t99 / 0.2e1 + t97 * t272 - t282 * t125 / 0.4e1 + (t71 + t75) * t143 / 0.4e1 - (t69 / 0.4e1 - t76 / 0.4e1) * t185 + (-t93 / 0.4e1 + t199 / 0.4e1) * t127;
t180 = Ifges(6,3) * t256 + t37 * mrSges(6,1) / 0.2e1 - t38 * mrSges(6,2) / 0.2e1 + t192;
t2 = t176 - t180;
t183 = t2 * qJD(1) + t17 * qJD(2) - t11 * qJD(3);
t178 = (-t125 * t265 + t127 * t263) * mrSges(6,3) + t195 * t277 + (-t125 * t97 - t127 * t96 + t143 * t36 - t185 * t35) * t275 + t99 * t265 + t101 * t262 - t189;
t181 = -t141 * t276 + t247 / 0.2e1 + t244 / 0.2e1;
t10 = (t175 * t278 + t235 / 0.2e1 + t238 / 0.2e1) * t173 + t178 + t181;
t21 = (t143 ^ 2 + t185 ^ 2) * mrSges(6,3) + t204 + m(6) * (t143 * t97 - t185 * t96) + m(5) * t202;
t26 = t184 + (t276 + (t169 / 0.2e1 + t170 / 0.2e1 - 0.1e1 / 0.2e1) * m(5)) * t173;
t182 = qJD(1) * t10 + qJD(2) * t26 + qJD(3) * t21;
t92 = t239 - t240;
t30 = qJD(4) * t212;
t25 = m(6) * t258 + t214 * t211 + t184 + t211;
t18 = t210 + t191;
t9 = t175 * t211 - mrSges(5,2) * t220 / 0.2e1 + mrSges(5,1) * t207 + t178 - t181;
t8 = t179 - t190;
t6 = t200 * t258 - t150 * t222 / 0.2e1 + t149 * t207 + t148 * t219 / 0.2e1 + t147 * t220 / 0.2e1 + t280 * t256 + t186 - t177 + (t74 - t280) * t257;
t3 = t176 + t180;
t12 = [qJD(2) * t16 + qJD(3) * t1 + qJD(4) * t14 + qJD(5) * t4, t225 + m(6) * (-t124 * t143 + t126 * t185) * qJD(2) + t6 * qJD(3) + t30 + t8 * qJD(5), t248 + t6 * qJD(2) + t9 * qJD(4) + t3 * qJD(5) + (-Ifges(4,5) + (Ifges(5,1) * t171 + t250) * t260 + (Ifges(5,2) * t172 + t251) * t261 + (-m(5) * pkin(3) + t230) * t175) * t213 + (-Ifges(4,6) * t174 + t123 * t261 + t122 * t259 + t165 * t74 + t70 * t263 - t141 * t92 + t68 * t265 + pkin(3) * t280 + t94 * t268 + t93 * t271 + t97 * t98 + t96 * t100 + m(6) * (-t141 * t165 + t37 * t96 + t38 * t97) - mrSges(4,2) * t217 + (m(5) * t194 + t221 - t224) * qJ(4) + (Ifges(5,5) * t171 + Ifges(6,5) * t185 + Ifges(5,6) * t172 + Ifges(6,6) * t143) * t256 + (t38 * t143 - t185 * t37) * mrSges(6,3) + t194 * mrSges(5,3)) * qJD(3), t9 * qJD(3) + t187, t232 + t8 * qJD(2) + t3 * qJD(3) + (-mrSges(6,1) * t36 - mrSges(6,2) * t35 + t216) * qJD(5); -qJD(3) * t5 + qJD(5) * t7 - t225 + t30, t24 * qJD(3), t25 * qJD(4) + t18 * qJD(5) + (t92 + t230) * t213 + t198 + ((t125 * t185 + t127 * t143) * mrSges(6,3) + (-mrSges(4,2) + t204) * t174 + 0.2e1 * (t165 * t173 + t241 - t245) * t275 + 0.2e1 * (t174 * t202 - t253) * t277) * qJD(3), qJD(1) * t212 + t25 * qJD(3), t231 + t18 * qJD(3) + (-mrSges(6,1) * t126 + mrSges(6,2) * t124) * qJD(5); qJD(2) * t5 + qJD(4) * t10 + qJD(5) * t2 - t248, qJD(4) * t26 + qJD(5) * t17 - t198, qJD(4) * t21 - qJD(5) * t11, t182, (-mrSges(6,1) * t97 - mrSges(6,2) * t96 + t215) * qJD(5) + t183; -t10 * qJD(3) + t73 * qJD(5) - t187, -qJD(1) * t273 / 0.2e1 - t26 * qJD(3), qJD(5) * t91 - t182, 0, t197; -qJD(2) * t7 - qJD(3) * t2 - qJD(4) * t73 - t232, -t17 * qJD(3) - t231, -qJD(4) * t91 - t183, -t197, 0;];
Cq = t12;
