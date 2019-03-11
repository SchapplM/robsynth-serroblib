% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:07:49
% EndTime: 2019-03-09 03:08:05
% DurationCPUTime: 8.61s
% Computational Cost: add. (4957->523), mult. (12199->695), div. (0->0), fcn. (7972->8), ass. (0->229)
t303 = Ifges(6,1) + Ifges(7,1);
t302 = -Ifges(6,4) + Ifges(7,5);
t301 = Ifges(7,4) + Ifges(6,5);
t185 = sin(pkin(10));
t187 = cos(pkin(10));
t189 = sin(qJ(5));
t191 = cos(qJ(5));
t156 = t185 * t191 + t187 * t189;
t192 = cos(qJ(3));
t199 = t156 * t192;
t123 = qJD(1) * t199;
t146 = t156 * qJD(5);
t243 = t123 - t146;
t249 = t185 * t189;
t155 = -t187 * t191 + t249;
t198 = t155 * t192;
t124 = qJD(1) * t198;
t145 = t155 * qJD(5);
t242 = -t124 + t145;
t300 = -Ifges(6,6) + Ifges(7,6);
t316 = Ifges(6,3) + Ifges(7,2);
t315 = -Ifges(4,6) / 0.2e1;
t234 = t185 * qJD(3);
t190 = sin(qJ(3));
t241 = qJD(1) * t190;
t153 = t187 * t241 + t234;
t314 = t153 / 0.2e1;
t240 = qJD(1) * t192;
t313 = -t240 / 0.2e1;
t232 = qJD(1) * qJD(3);
t216 = t190 * t232;
t195 = qJD(3) * t198;
t239 = qJD(3) * t187;
t200 = t185 * t241 - t239;
t99 = t153 * t189 + t191 * t200;
t58 = -qJD(1) * t195 - qJD(5) * t99;
t194 = t153 * t191 - t189 * t200;
t196 = qJD(3) * t199;
t59 = qJD(1) * t196 + qJD(5) * t194;
t312 = t216 * t301 + t302 * t59 + t303 * t58;
t178 = qJD(5) - t240;
t270 = Ifges(7,5) * t99;
t95 = Ifges(6,4) * t99;
t299 = t178 * t301 + t194 * t303 + t270 - t95;
t180 = sin(pkin(9)) * pkin(1) + pkin(7);
t168 = t180 * qJD(1);
t138 = qJD(2) * t190 + t168 * t192;
t220 = t185 * t240;
t111 = pkin(4) * t220 + t138;
t311 = -pkin(5) * t243 + qJ(6) * t242 - qJD(6) * t156 - t111;
t288 = t59 / 0.2e1;
t310 = Ifges(7,3) * t288;
t309 = -Ifges(4,1) / 0.2e1;
t308 = -Ifges(6,6) / 0.2e1;
t307 = Ifges(7,6) / 0.2e1;
t290 = t58 / 0.2e1;
t286 = -t99 / 0.2e1;
t285 = t99 / 0.2e1;
t275 = t178 / 0.2e1;
t306 = Ifges(4,4) * t313;
t279 = t194 / 0.2e1;
t305 = qJD(3) * t315;
t304 = Ifges(5,5) * t314;
t120 = qJD(3) * qJ(4) + t138;
t221 = -cos(pkin(9)) * pkin(1) - pkin(2);
t151 = -pkin(3) * t192 - qJ(4) * t190 + t221;
t125 = t151 * qJD(1);
t67 = -t120 * t185 + t125 * t187;
t42 = -pkin(4) * t240 - pkin(8) * t153 + t67;
t68 = t120 * t187 + t125 * t185;
t49 = -pkin(8) * t200 + t68;
t12 = -t189 * t49 + t191 * t42;
t298 = qJD(6) - t12;
t297 = t216 * t316 + t300 * t59 + t301 * t58;
t215 = t192 * t232;
t214 = t185 * t215;
t262 = mrSges(5,2) * t187;
t129 = mrSges(5,1) * t214 + t215 * t262;
t19 = mrSges(7,1) * t59 - t58 * mrSges(7,3);
t20 = t59 * mrSges(6,1) + mrSges(6,2) * t58;
t296 = -t129 - t19 - t20;
t13 = t189 * t42 + t191 * t49;
t245 = t187 * t192;
t202 = pkin(4) * t190 - pkin(8) * t245;
t197 = t202 * qJD(3);
t157 = t190 * t168;
t233 = t192 * qJD(2);
t182 = qJD(3) * t233;
t117 = t182 + (qJD(4) - t157) * qJD(3);
t209 = pkin(3) * t190 - qJ(4) * t192;
t142 = qJD(3) * t209 - qJD(4) * t190;
t126 = t142 * qJD(1);
t62 = -t117 * t185 + t126 * t187;
t43 = qJD(1) * t197 + t62;
t63 = t117 * t187 + t126 * t185;
t50 = -pkin(8) * t214 + t63;
t4 = -qJD(5) * t13 - t189 * t50 + t191 * t43;
t136 = t187 * t151;
t246 = t187 * t190;
t82 = -pkin(8) * t246 + t136 + (-t180 * t185 - pkin(4)) * t192;
t103 = t151 * t185 + t180 * t245;
t248 = t185 * t190;
t86 = -pkin(8) * t248 + t103;
t263 = t189 * t82 + t191 * t86;
t238 = qJD(3) * t190;
t219 = t180 * t238;
t96 = t142 * t187 + t185 * t219;
t70 = t197 + t96;
t130 = t185 * t142;
t172 = t190 * t180;
t247 = t185 * t192;
t81 = t130 + (-pkin(8) * t247 - t172 * t187) * qJD(3);
t9 = -qJD(5) * t263 - t189 * t81 + t191 * t70;
t235 = qJD(5) * t191;
t237 = qJD(5) * t189;
t3 = t189 * t43 + t191 * t50 + t235 * t42 - t237 * t49;
t1 = qJ(6) * t216 + qJD(6) * t178 + t3;
t2 = -pkin(5) * t216 - t4;
t295 = -t4 * mrSges(6,1) + t2 * mrSges(7,1) + t3 * mrSges(6,2) - t1 * mrSges(7,3);
t10 = -pkin(5) * t178 + t298;
t11 = qJ(6) * t178 + t13;
t170 = t221 * qJD(1);
t227 = t307 + t308;
t228 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t229 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t294 = t228 * t178 + t229 * t194 + t227 * t99 + t11 * mrSges(7,3) + t12 * mrSges(6,1) + t170 * mrSges(4,1) + t67 * mrSges(5,1) + t305 - (Ifges(4,4) * t190 + t192 * Ifges(4,2)) * qJD(1) / 0.2e1 + Ifges(6,6) * t286 + Ifges(7,6) * t285 - Ifges(5,6) * t200 / 0.2e1 + Ifges(5,3) * t313 + t304 - t10 * mrSges(7,1) - t13 * mrSges(6,2) - t68 * mrSges(5,2) + t301 * t279 + t316 * t275;
t293 = Ifges(5,1) / 0.2e1;
t292 = Ifges(7,5) * t290 + t216 * t307 + t310;
t291 = -t58 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t288 + t216 * t308;
t289 = -t59 / 0.2e1;
t280 = -t194 / 0.2e1;
t276 = -t178 / 0.2e1;
t274 = -t185 / 0.2e1;
t273 = t187 / 0.2e1;
t271 = mrSges(6,3) * t99;
t268 = pkin(8) + qJ(4);
t44 = -mrSges(7,2) * t59 + mrSges(7,3) * t216;
t47 = -mrSges(6,2) * t216 - mrSges(6,3) * t59;
t267 = t44 + t47;
t45 = mrSges(6,1) * t216 - mrSges(6,3) * t58;
t46 = -mrSges(7,1) * t216 + mrSges(7,2) * t58;
t266 = t46 - t45;
t137 = -t157 + t233;
t160 = t209 * qJD(1);
t87 = -t137 * t185 + t160 * t187;
t64 = qJD(1) * t202 + t87;
t88 = t137 * t187 + t160 * t185;
t72 = -pkin(8) * t220 + t88;
t25 = t189 * t64 + t191 * t72;
t77 = -mrSges(6,2) * t178 - t271;
t80 = -mrSges(7,2) * t99 + mrSges(7,3) * t178;
t265 = t77 + t80;
t261 = mrSges(6,3) * t194;
t78 = mrSges(6,1) * t178 - t261;
t79 = -mrSges(7,1) * t178 + mrSges(7,2) * t194;
t264 = t79 - t78;
t260 = Ifges(5,1) * t153;
t259 = Ifges(5,4) * t187;
t258 = Ifges(6,4) * t194;
t256 = Ifges(5,5) * t187;
t255 = Ifges(5,2) * t185;
t254 = Ifges(5,6) * t185;
t253 = Ifges(4,5) * qJD(3);
t128 = t138 * qJD(3);
t251 = t128 * t192;
t250 = t180 * t192;
t226 = mrSges(4,3) * t241;
t244 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t200 + mrSges(5,2) * t153 + t226;
t132 = pkin(4) * t192 * t234 + qJD(3) * t250;
t141 = pkin(4) * t248 + t172;
t236 = qJD(5) * t190;
t225 = mrSges(4,3) * t240;
t224 = Ifges(5,5) * t240;
t223 = Ifges(5,6) * t240;
t222 = t128 * t172;
t181 = -pkin(4) * t187 - pkin(3);
t213 = t234 / 0.2e1 - t153 / 0.2e1;
t212 = mrSges(5,1) * t185 + t262;
t211 = Ifges(5,1) * t187 - Ifges(5,4) * t185;
t210 = -t255 + t259;
t208 = -t185 * t62 + t187 * t63;
t207 = -t185 * t67 + t187 * t68;
t24 = -t189 * t72 + t191 * t64;
t28 = -t189 * t86 + t191 * t82;
t163 = t268 * t185;
t164 = t268 * t187;
t203 = -t163 * t191 - t164 * t189;
t110 = -t163 * t189 + t164 * t191;
t8 = t189 * t70 + t191 * t81 + t235 * t82 - t237 * t86;
t201 = t137 * mrSges(4,3) + t241 * t309 + t306 - t253 / 0.2e1 - t170 * mrSges(4,2);
t118 = -qJD(3) * pkin(3) + qJD(4) - t137;
t104 = pkin(4) * t214 + t128;
t89 = pkin(4) * t200 + t118;
t171 = -qJD(3) * mrSges(4,2) + t225;
t140 = (mrSges(5,1) * t190 - mrSges(5,3) * t245) * t232;
t139 = (-mrSges(5,2) * t190 - mrSges(5,3) * t247) * t232;
t134 = t155 * t190;
t133 = t156 * t190;
t127 = -t168 * t238 + t182;
t122 = -mrSges(5,1) * t240 - mrSges(5,3) * t153;
t121 = mrSges(5,2) * t240 - mrSges(5,3) * t200;
t108 = (Ifges(5,5) * t190 + t192 * t211) * t232;
t107 = (Ifges(5,6) * t190 + t192 * t210) * t232;
t102 = -t180 * t247 + t136;
t97 = -t187 * t219 + t130;
t94 = Ifges(7,5) * t194;
t93 = pkin(5) * t155 - qJ(6) * t156 + t181;
t92 = -Ifges(5,4) * t200 - t224 + t260;
t91 = Ifges(5,4) * t153 - Ifges(5,2) * t200 - t223;
t85 = t235 * t246 - t236 * t249 + t196;
t84 = -t156 * t236 - t195;
t76 = qJD(4) * t156 + qJD(5) * t110;
t75 = -qJD(4) * t155 + qJD(5) * t203;
t60 = pkin(5) * t133 + qJ(6) * t134 + t141;
t41 = mrSges(6,1) * t99 + mrSges(6,2) * t194;
t40 = mrSges(7,1) * t99 - mrSges(7,3) * t194;
t39 = pkin(5) * t194 + qJ(6) * t99;
t33 = -Ifges(6,2) * t99 + Ifges(6,6) * t178 + t258;
t30 = Ifges(7,6) * t178 + Ifges(7,3) * t99 + t94;
t27 = pkin(5) * t192 - t28;
t26 = -qJ(6) * t192 + t263;
t23 = pkin(5) * t99 - qJ(6) * t194 + t89;
t22 = -pkin(5) * t241 - t24;
t21 = qJ(6) * t241 + t25;
t18 = pkin(5) * t85 - qJ(6) * t84 + qJD(6) * t134 + t132;
t7 = pkin(5) * t59 - qJ(6) * t58 - qJD(6) * t194 + t104;
t6 = -pkin(5) * t238 - t9;
t5 = qJ(6) * t238 - qJD(6) * t192 + t8;
t14 = [m(4) * (t127 * t250 + t222) - t297 * t192 / 0.2e1 + (mrSges(6,1) * t104 + mrSges(7,1) * t7 - mrSges(7,2) * t1 - mrSges(6,3) * t3 - Ifges(6,2) * t289 + t290 * t302 + t291 + t292 + t310) * t133 + (t180 * t129 + t107 * t274 + t108 * t273 + (mrSges(4,3) + t212) * t128 + (-t185 * t63 - t187 * t62) * mrSges(5,3) + ((Ifges(5,6) * t273 + t315) * qJD(3) - t180 * t171 + (-m(4) * t180 - mrSges(4,3)) * t138 + t304 + t294) * qJD(3) + (t221 * mrSges(4,1) + (-0.3e1 / 0.2e1 * Ifges(4,4) + t256 / 0.2e1 - t254) * t190 - t229 * t134 + t227 * t133) * t232) * t190 + m(5) * (t102 * t62 + t103 * t63 + t67 * t96 + t68 * t97 + t222) + m(6) * (t104 * t141 + t12 * t9 + t13 * t8 + t132 * t89 + t263 * t3 + t28 * t4) + t263 * t47 + m(7) * (t1 * t26 + t10 * t6 + t11 * t5 + t18 * t23 + t2 * t27 + t60 * t7) - (mrSges(6,2) * t104 + mrSges(7,2) * t2 - mrSges(6,3) * t4 - mrSges(7,3) * t7 + Ifges(6,4) * t289 + Ifges(7,5) * t288 + t290 * t303) * t134 + (mrSges(6,2) * t89 + t10 * mrSges(7,2) - t12 * mrSges(6,3) - mrSges(7,3) * t23 + Ifges(6,4) * t286 + Ifges(7,5) * t285 + t299 / 0.2e1 + t301 * t275 + t303 * t279) * t84 + (t63 * mrSges(5,2) - t62 * mrSges(5,1) + (t118 * t212 + t211 * t314 + (Ifges(4,5) / 0.2e1 + t210 * t273) * qJD(3) + t91 * t274 + t92 * t273 + (-t185 * t68 - t187 * t67) * mrSges(5,3) + (-m(4) * t137 + m(5) * t118 + t244) * t180 - t201) * qJD(3) + t127 * mrSges(4,3) - Ifges(7,6) * t288 - Ifges(6,6) * t289 - t301 * t290 + (t221 * mrSges(4,2) + (0.3e1 / 0.2e1 * t254 - 0.3e1 / 0.2e1 * t256 + 0.3e1 / 0.2e1 * Ifges(4,4)) * t192 + (t187 ^ 2 * t293 - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(5,3) + (-0.3e1 / 0.2e1 * t259 + t255) * t185 - t228) * t190) * t232 + t295) * t192 + t18 * t40 + t26 * t44 + t28 * t45 + t27 * t46 + t60 * t19 + t8 * t77 + t9 * t78 - t312 * t134 / 0.2e1 + t6 * t79 + t5 * t80 + (-t11 * mrSges(7,2) - t13 * mrSges(6,3) + Ifges(7,3) * t285 - Ifges(6,2) * t286 + t30 / 0.2e1 + t23 * mrSges(7,1) - t33 / 0.2e1 + t89 * mrSges(6,1) + t300 * t275 + t302 * t279) * t85 + t97 * t121 + t96 * t122 + t132 * t41 + t103 * t139 + t102 * t140 + t141 * t20; t264 * t85 + t265 * t84 + (t139 * t187 - t140 * t185) * t190 - t267 * t134 + t266 * t133 + t296 * t192 + ((t121 * t187 - t122 * t185 + t171 - t225) * t192 + (t40 + t41 - t226 + t244) * t190) * qJD(3) + m(7) * (-t1 * t134 + t10 * t85 + t11 * t84 + t133 * t2 - t192 * t7 + t23 * t238) + m(6) * (-t104 * t192 - t12 * t85 + t13 * t84 - t133 * t4 - t134 * t3 + t238 * t89) + m(4) * (t127 * t190 - t251 + (-t137 * t190 + t138 * t192) * qJD(3)) + m(5) * (-t251 + t208 * t190 + (t118 * t190 + t192 * t207) * qJD(3)); (t63 * mrSges(5,3) + qJD(4) * t121 + qJ(4) * t139 + t107 / 0.2e1 - t128 * mrSges(5,1)) * t187 + (t155 * t302 + t156 * t303) * t290 + (t123 * t302 - t124 * t303) * t280 + (-t145 * t303 + t146 * t302) * t279 + t264 * t76 + t265 * t75 + (-pkin(3) * t128 + qJ(4) * t208 + qJD(4) * t207 - t118 * t138 - t67 * t87 - t68 * t88) * m(5) + (t33 - t30) * (t123 / 0.2e1 - t146 / 0.2e1) + (t104 * t181 + t203 * t4 + t110 * t3 - t111 * t89 + (-t25 + t75) * t13 + (-t24 - t76) * t12) * m(6) - t266 * t203 + (-Ifges(6,4) * t124 - Ifges(7,5) * t145 - Ifges(6,2) * t123 + Ifges(7,3) * t146) * t285 + (-Ifges(6,4) * t145 - Ifges(7,5) * t124 - Ifges(6,2) * t146 + Ifges(7,3) * t123) * t286 + t267 * t110 + (t12 * t242 + t13 * t243 - t155 * t3 - t156 * t4) * mrSges(6,3) + (-t1 * t155 - t10 * t242 + t11 * t243 + t156 * t2) * mrSges(7,2) + (-mrSges(7,1) * t243 + mrSges(7,3) * t242) * t23 + (-mrSges(6,1) * t243 - mrSges(6,2) * t242) * t89 - t244 * t138 + (-t62 * mrSges(5,3) - qJD(4) * t122 - qJ(4) * t140 + t108 / 0.2e1 + t128 * mrSges(5,2)) * t185 + ((t306 + t253 / 0.2e1 + (-t118 * mrSges(5,2) + t67 * mrSges(5,3) - t260 / 0.2e1 - t92 / 0.2e1 + t224 / 0.2e1) * t187 + (-t118 * mrSges(5,1) + t68 * mrSges(5,3) + t239 * t293 + t91 / 0.2e1 - t223 / 0.2e1 - t213 * Ifges(5,4)) * t185 + t201) * t192 + ((t309 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t185 * t210 / 0.2e1) * t240 + (Ifges(4,4) / 0.2e1 + t254 / 0.2e1) * t241 + t213 * Ifges(5,5) + t138 * mrSges(4,3) + t305 + (t155 * t300 + t156 * t301) * qJD(3) / 0.2e1 - t294) * t190) * qJD(1) + (Ifges(7,5) * t156 + Ifges(7,3) * t155) * t288 + (Ifges(6,4) * t156 - Ifges(6,2) * t155) * t289 + t155 * t291 + t155 * t292 - t25 * t77 - t24 * t78 + (t1 * t110 - t203 * t2 + t7 * t93 + t311 * t23 + (-t21 + t75) * t11 + (-t22 + t76) * t10) * m(7) + t311 * t40 + t299 * (t124 / 0.2e1 - t145 / 0.2e1) + t312 * t156 / 0.2e1 + (t123 * t300 - t124 * t301) * t276 + (-t145 * t301 + t146 * t300) * t275 - t22 * t79 - t21 * t80 + t93 * t19 - t111 * t41 - t88 * t121 - t87 * t122 - t127 * mrSges(4,2) - t128 * mrSges(4,1) - pkin(3) * t129 + t7 * (mrSges(7,1) * t155 - mrSges(7,3) * t156) + t104 * (mrSges(6,1) * t155 + mrSges(6,2) * t156) - t137 * t171 + t181 * t20; t265 * t99 + t200 * t121 + t153 * t122 - t264 * t194 + (-t10 * t194 + t11 * t99 + t7) * m(7) + (t12 * t194 + t13 * t99 + t104) * m(6) + (t67 * t153 + t200 * t68 + t128) * m(5) - t296; (-pkin(5) * t2 + qJ(6) * t1 - t10 * t13 + t11 * t298 - t23 * t39) * m(7) + (-t265 - t271) * t12 - t295 + (t261 - t264) * t13 + (t10 * t99 + t11 * t194) * mrSges(7,2) + (Ifges(7,3) * t194 - t270) * t286 - t23 * (mrSges(7,1) * t194 + mrSges(7,3) * t99) - t89 * (mrSges(6,1) * t194 - mrSges(6,2) * t99) + (t194 * t300 - t301 * t99) * t276 + (-Ifges(6,2) * t194 + t299 - t95) * t285 + t33 * t279 - t39 * t40 + qJ(6) * t44 - pkin(5) * t46 + (-t303 * t99 - t258 + t30 + t94) * t280 + qJD(6) * t80 + t297; t194 * t40 - t178 * t80 + 0.2e1 * (t2 / 0.2e1 + t23 * t279 + t11 * t276) * m(7) + t46;];
tauc  = t14(:);
