% Calculate vector of inverse dynamics joint torques for
% S5RRRPP5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:36
% EndTime: 2019-12-31 20:57:54
% DurationCPUTime: 9.64s
% Computational Cost: add. (3023->453), mult. (6764->543), div. (0->0), fcn. (4210->8), ass. (0->205)
t316 = Ifges(6,4) + Ifges(5,5);
t317 = Ifges(5,4) + Ifges(4,5);
t328 = -Ifges(6,5) + t317;
t330 = mrSges(4,1) + mrSges(5,1);
t329 = -mrSges(6,2) - mrSges(5,3);
t315 = Ifges(6,2) + Ifges(5,3);
t324 = Ifges(5,6) - Ifges(6,6);
t298 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t183 = sin(qJ(3));
t184 = sin(qJ(2));
t186 = cos(qJ(2));
t276 = cos(qJ(3));
t129 = t183 * t186 + t184 * t276;
t114 = t129 * qJD(1);
t327 = t316 * t114;
t226 = t276 * t186;
t243 = qJD(1) * t184;
t113 = -qJD(1) * t226 + t183 * t243;
t326 = t316 * t113;
t314 = -Ifges(4,6) + Ifges(5,6);
t181 = qJD(2) + qJD(3);
t323 = t315 * t113 + t324 * t181 + t327;
t179 = t186 * pkin(2);
t313 = t179 + pkin(1);
t142 = t313 * qJD(1);
t322 = -qJ(4) * t114 - t142;
t239 = qJD(1) * qJD(2);
t137 = qJDD(1) * t186 - t184 * t239;
t108 = Ifges(4,4) * t113;
t321 = t298 * t114 + t328 * t181 - t108 + t326;
t182 = qJ(2) + qJ(3);
t177 = sin(t182);
t178 = cos(t182);
t320 = -t330 * t178 + (mrSges(4,2) + t329) * t177;
t318 = t137 / 0.2e1;
t277 = t186 / 0.2e1;
t180 = qJDD(2) + qJDD(3);
t138 = qJDD(1) * t184 + t186 * t239;
t196 = t129 * qJD(3);
t56 = qJD(1) * t196 - t137 * t276 + t183 * t138;
t33 = mrSges(6,2) * t180 + mrSges(6,3) * t56;
t35 = -mrSges(5,2) * t56 + mrSges(5,3) * t180;
t312 = t33 + t35;
t263 = mrSges(6,3) * t113;
t93 = mrSges(6,2) * t181 + t263;
t267 = mrSges(5,2) * t113;
t98 = mrSges(5,3) * t181 - t267;
t311 = t93 + t98;
t265 = mrSges(4,3) * t113;
t94 = -mrSges(4,2) * t181 - t265;
t310 = -t94 - t98;
t264 = mrSges(4,3) * t114;
t266 = mrSges(5,2) * t114;
t309 = t330 * t181 - t264 - t266;
t308 = t186 * Ifges(3,2);
t188 = -pkin(7) - pkin(6);
t143 = t188 * t184;
t133 = qJD(1) * t143;
t124 = qJD(2) * pkin(2) + t133;
t144 = t188 * t186;
t134 = qJD(1) * t144;
t247 = t183 * t134;
t80 = t276 * t124 + t247;
t225 = qJD(4) - t80;
t92 = t183 * t143 - t276 * t144;
t187 = cos(qJ(1));
t248 = t178 * t187;
t250 = t177 * t187;
t307 = pkin(3) * t248 + qJ(4) * t250;
t185 = sin(qJ(1));
t249 = t178 * t185;
t145 = qJ(4) * t249;
t306 = -m(6) * t145 + t329 * t249;
t147 = qJ(4) * t248;
t305 = -m(6) * t147 + t329 * t248;
t189 = -pkin(3) - pkin(4);
t213 = m(6) * t189 - mrSges(6,1);
t203 = t213 * t177;
t274 = pkin(2) * t184;
t302 = m(6) * t274 - m(5) * (-pkin(3) * t177 - t274) + t177 * mrSges(5,1) - t203;
t242 = qJD(1) * t186;
t271 = pkin(6) * t186;
t272 = pkin(6) * t184;
t301 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t243) * t271 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t242) * t272;
t126 = t137 * pkin(6);
t127 = t138 * pkin(6);
t300 = t126 * t186 + t127 * t184;
t299 = g(1) * t187 + g(2) * t185;
t297 = -Ifges(4,4) + t316;
t201 = -t183 * t184 + t226;
t195 = t201 * qJD(3);
t55 = qJD(1) * t195 + t183 * t137 + t138 * t276;
t32 = -t180 * mrSges(5,1) + t55 * mrSges(5,2);
t223 = qJD(3) * t276;
t240 = qJD(3) * t183;
t89 = qJDD(2) * pkin(2) - pkin(7) * t138 - t127;
t90 = pkin(7) * t137 + t126;
t14 = -t124 * t240 + t134 * t223 - t183 * t90 + t276 * t89;
t198 = qJDD(4) - t14;
t9 = -t180 * pkin(3) + t198;
t296 = m(5) * t9 + t32;
t295 = -t178 * mrSges(6,1) + t320;
t58 = -pkin(3) * t181 + t225;
t294 = m(5) * t58 - t309;
t293 = mrSges(2,2) - mrSges(5,2) - mrSges(4,3) - m(3) * pkin(6) - mrSges(3,3) - m(6) * (-qJ(5) - t188) + mrSges(6,3);
t141 = -mrSges(3,1) * t186 + mrSges(3,2) * t184;
t292 = -m(3) * pkin(1) - mrSges(2,1) + t141 + t320;
t286 = -t113 / 0.2e1;
t285 = t113 / 0.2e1;
t283 = t114 / 0.2e1;
t282 = -t201 / 0.2e1;
t279 = -t181 / 0.2e1;
t278 = t181 / 0.2e1;
t275 = pkin(2) * t183;
t273 = pkin(4) * t114;
t165 = t178 * pkin(3);
t262 = Ifges(3,4) * t184;
t261 = Ifges(3,4) * t186;
t258 = t114 * mrSges(6,3);
t257 = t114 * Ifges(4,4);
t255 = qJ(5) * t113;
t254 = qJDD(1) * pkin(1);
t157 = t177 * qJ(4);
t251 = t177 * t185;
t72 = t114 * pkin(3) + t113 * qJ(4);
t84 = t276 * t133 + t247;
t227 = t276 * t134;
t81 = t183 * t124 - t227;
t244 = t165 + t157;
t241 = qJD(2) * t184;
t235 = t276 * pkin(2);
t234 = pkin(2) * t243;
t233 = pkin(2) * t241;
t230 = t179 + t244;
t228 = qJD(2) * t188;
t224 = -t56 * mrSges(6,1) + t55 * mrSges(6,2);
t30 = -t180 * mrSges(6,1) - t55 * mrSges(6,3);
t218 = t239 / 0.2e1;
t216 = -t313 - t157;
t83 = t133 * t183 - t227;
t91 = -t276 * t143 - t144 * t183;
t148 = t187 * t313;
t215 = -t185 * t188 + t148;
t214 = pkin(2) * t223;
t167 = -t235 - pkin(3);
t212 = qJ(4) * t129 + t313;
t44 = t81 + t255;
t209 = mrSges(3,1) * t184 + mrSges(3,2) * t186;
t207 = mrSges(4,1) * t177 + mrSges(4,2) * t178;
t205 = t262 + t308;
t204 = Ifges(3,5) * t186 - Ifges(3,6) * t184;
t105 = -pkin(2) * t137 - t254;
t59 = t234 + t72;
t202 = pkin(1) * t209;
t200 = t184 * (Ifges(3,1) * t186 - t262);
t13 = t124 * t223 + t134 * t240 + t183 * t89 + t276 * t90;
t135 = t184 * t228;
t136 = t186 * t228;
t19 = t276 * t135 + t183 * t136 + t143 * t223 + t144 * t240;
t8 = t180 * qJ(4) + t181 * qJD(4) + t13;
t85 = qJD(2) * t201 + t195;
t194 = qJ(4) * t85 + qJD(4) * t129 - t233;
t191 = qJ(4) * t55 + qJD(4) * t114 - t105;
t20 = qJD(3) * t92 + t183 * t135 - t276 * t136;
t104 = t114 * qJ(5);
t17 = t181 * t189 - t104 + t225;
t18 = t113 * t189 + qJD(5) - t322;
t176 = t181 * qJ(4);
t28 = t176 + t44;
t3 = -t55 * qJ(5) - t114 * qJD(5) + t180 * t189 + t198;
t4 = qJ(5) * t56 + qJD(5) * t113 + t8;
t57 = pkin(3) * t113 + t322;
t65 = -t113 * Ifges(4,2) + t181 * Ifges(4,6) + t257;
t71 = t176 + t81;
t190 = -t9 * mrSges(5,1) - t13 * mrSges(4,2) + t14 * mrSges(4,1) + t4 * mrSges(6,2) + t8 * mrSges(5,3) - t3 * mrSges(6,1) - t57 * (mrSges(5,1) * t114 + mrSges(5,3) * t113) - t18 * (-mrSges(6,1) * t114 - mrSges(6,2) * t113) + t142 * (mrSges(4,1) * t114 - mrSges(4,2) * t113) - t17 * t263 + t81 * t264 - t80 * t265 + t71 * t266 + t58 * t267 + (-Ifges(6,5) * t113 + Ifges(6,6) * t114) * t278 + t65 * t283 - t28 * t258 + (t114 * t315 - t326) * t286 + (-t113 * t317 + t114 * t314) * t279 + (-Ifges(4,6) + t324) * t56 + t328 * t55 + (Ifges(5,2) + Ifges(4,3) + Ifges(6,3)) * t180 + (-Ifges(4,2) * t114 - t108 + t321) * t285 - (-t298 * t113 - t257 + t323 + t327) * t114 / 0.2e1;
t170 = Ifges(3,4) * t242;
t164 = t178 * pkin(4);
t163 = qJ(4) + t275;
t156 = -pkin(4) + t167;
t154 = t214 + qJD(4);
t112 = Ifges(3,1) * t243 + Ifges(3,5) * qJD(2) + t170;
t111 = Ifges(3,6) * qJD(2) + qJD(1) * t205;
t95 = -mrSges(6,1) * t181 - t258;
t86 = qJD(2) * t129 + t196;
t79 = -pkin(3) * t201 - t212;
t78 = mrSges(4,1) * t113 + mrSges(4,2) * t114;
t77 = -mrSges(6,1) * t113 + mrSges(6,2) * t114;
t76 = mrSges(5,1) * t113 - mrSges(5,3) * t114;
t70 = -qJ(5) * t201 + t92;
t69 = -qJ(5) * t129 + t91;
t54 = t104 + t84;
t53 = t83 + t255;
t52 = -t189 * t201 + t212;
t43 = t104 + t80;
t42 = -t72 - t273;
t34 = -mrSges(4,2) * t180 - mrSges(4,3) * t56;
t31 = mrSges(4,1) * t180 - mrSges(4,3) * t55;
t27 = -t59 - t273;
t15 = pkin(3) * t86 - t194;
t11 = -t85 * qJ(5) - t129 * qJD(5) + t20;
t10 = qJ(5) * t86 - qJD(5) * t201 + t19;
t6 = t189 * t86 + t194;
t5 = pkin(3) * t56 - t191;
t1 = t189 * t56 + qJDD(5) + t191;
t2 = [m(6) * (t1 * t52 + t10 * t28 + t11 * t17 + t18 * t6 + t3 * t69 + t4 * t70) + t78 * t233 + (-t142 * mrSges(4,2) - t57 * mrSges(5,3) + t18 * mrSges(6,2) - t17 * mrSges(6,3) + t58 * mrSges(5,2) - t80 * mrSges(4,3) + t298 * t283 + t321 / 0.2e1 + t316 * t285 + t317 * t278 + Ifges(6,5) * t279 + Ifges(4,4) * t286) * t85 + (-t129 * t14 + t13 * t201) * mrSges(4,3) + (t137 * t271 + t138 * t272 + t300) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t300) + (t112 * t277 + t204 * qJD(2) / 0.2e1 - t301) * qJD(2) + (-m(5) * (t215 + t307) - m(4) * t215 - m(6) * (t148 + t307) + (-(m(6) * pkin(4) + mrSges(6,1)) * t178 + t292) * t187 + t293 * t185) * g(2) - pkin(1) * (-mrSges(3,1) * t137 + mrSges(3,2) * t138) + t10 * t93 + t11 * t95 + (m(4) * t13 + m(5) * t8 + t34 + t35) * t92 - t180 * (Ifges(6,5) * t129 - Ifges(6,6) * t201) / 0.2e1 + t5 * (-mrSges(5,1) * t201 - mrSges(5,3) * t129) + t1 * (mrSges(6,1) * t201 + mrSges(6,2) * t129) + t105 * (-mrSges(4,1) * t201 + mrSges(4,2) * t129) + (t129 * t316 - t201 * t315) * t56 / 0.2e1 + (-t313 * mrSges(4,2) - t79 * mrSges(5,3) + t129 * t298 + t316 * t282 + (-t297 / 0.2e1 + Ifges(4,4) / 0.2e1) * t201) * t55 + (-Ifges(4,2) * t201 - t313 * mrSges(4,1) - Ifges(4,4) * t129 / 0.2e1 + t79 * mrSges(5,1) + t315 * t282) * t56 + (((m(4) + m(5)) * t188 + t293) * t187 + (m(4) * t313 - m(6) * t216 - t178 * t213 - m(5) * (t216 - t165) - t292) * t185) * g(1) + m(4) * (-t105 * t313 - t142 * t233) + (t129 * t317 + (Ifges(4,6) - t314) * t201) * t180 / 0.2e1 + t69 * t30 + t70 * t33 + t15 * t76 + t6 * t77 + (-mrSges(3,1) * t272 - mrSges(3,2) * t271 + 0.2e1 * Ifges(3,6) * t277) * qJDD(2) + (Ifges(3,1) * t138 + Ifges(3,4) * t318 + Ifges(3,5) * qJDD(2) - t218 * t308) * t184 + (-t142 * mrSges(4,1) + t57 * mrSges(5,1) - t18 * mrSges(6,1) - t65 / 0.2e1 + t28 * mrSges(6,3) - t71 * mrSges(5,2) - t81 * mrSges(4,3) + t323 / 0.2e1 + t297 * t283 + t315 * t285 + t314 * t278 + Ifges(6,6) * t279 - Ifges(4,2) * t286) * t86 + (-m(4) * t14 + t296 - t31) * t91 + (-m(4) * t80 + t294) * t20 + (-t129 * t3 - t201 * t4) * mrSges(6,3) + (t328 * t180 + t297 * t56) * t129 / 0.2e1 + (Ifges(3,4) * t138 + Ifges(3,2) * t137) * t277 + t205 * t318 + (t129 * t9 + t201 * t8) * mrSges(5,2) + (t186 * t261 + t200) * t218 + t324 * t180 * t282 + m(5) * (t15 * t57 + t5 * t79) + (m(4) * t81 + m(5) * t71 - t310) * t19 + t138 * t261 / 0.2e1 + Ifges(2,3) * qJDD(1) + t52 * t224 - t202 * t239 - t111 * t241 / 0.2e1 - t141 * t254; (t187 * t302 + t305) * g(1) + (t185 * t302 + t306) * g(2) + (t301 + (-t200 / 0.2e1 + t202) * qJD(1)) * qJD(1) + t167 * t32 + t156 * t30 + Ifges(3,6) * t137 + Ifges(3,5) * t138 - t126 * mrSges(3,2) - t127 * mrSges(3,1) - t54 * t93 - t53 * t95 + (m(6) * t17 + t294 + t95) * pkin(2) * t240 - t59 * t76 - t27 * t77 - (-Ifges(3,2) * t243 + t112 + t170) * t242 / 0.2e1 + t94 * t214 + (t141 - m(6) * (t164 + t230) - m(5) * t230 - m(4) * t179 + t295) * g(3) + ((t276 * t14 + t13 * t183 + (-t183 * t80 + t276 * t81) * qJD(3)) * pkin(2) + t142 * t234 + t80 * t83 - t81 * t84) * m(4) + (t156 * t3 + t163 * t4 - t17 * t53 - t18 * t27 + (t154 - t54) * t28) * m(6) + t190 + t309 * t83 + t310 * t84 + t311 * t154 + t312 * t163 + (-g(1) * t147 - g(2) * t145 + t163 * t8 + t167 * t9 - t57 * t59 - t58 * t83 + (t154 - t84) * t71) * m(5) + (m(4) * t274 + t207 + t209) * t299 + t31 * t235 + t34 * t275 - t78 * t234 - t204 * t239 / 0.2e1 + t111 * t243 / 0.2e1 + Ifges(3,3) * qJDD(2); -pkin(3) * t32 + t189 * t30 - t42 * t77 - t43 * t93 - t44 * t95 - t72 * t76 + t190 + t309 * t81 + t310 * t80 + t299 * t207 + t311 * qJD(4) + t312 * qJ(4) + (mrSges(5,1) * t251 - t185 * t203 + t306) * g(2) + (mrSges(5,1) * t250 - t187 * t203 + t305) * g(1) + (t4 * qJ(4) - t17 * t44 - t18 * t42 + t189 * t3 + (-t43 + qJD(4)) * t28) * m(6) + (-pkin(3) * t9 + qJ(4) * t8 - t57 * t72 - t58 * t81 - g(1) * (-pkin(3) * t250 + t147) - g(2) * (-pkin(3) * t251 + t145) + t225 * t71) * m(5) + (-m(6) * (t164 + t244) - m(5) * t244 + t295) * g(3); -t311 * t181 + (t76 - t77) * t114 - m(5) * (-t114 * t57 + t181 * t71) + t30 + t296 + (g(3) * t178 - t177 * t299) * (m(5) + m(6)) + (-t114 * t18 - t181 * t28 + t3) * m(6); -t113 * t93 + t114 * t95 + (g(1) * t185 - g(2) * t187 - t28 * t113 + t17 * t114 + t1) * m(6) + t224;];
tau = t2;
