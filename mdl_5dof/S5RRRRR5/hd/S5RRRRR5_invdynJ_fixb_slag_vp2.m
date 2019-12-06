% Calculate vector of inverse dynamics joint torques for
% S5RRRRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:09
% EndTime: 2019-12-05 18:58:16
% DurationCPUTime: 3.16s
% Computational Cost: add. (5576->374), mult. (7695->516), div. (0->0), fcn. (4376->16), ass. (0->200)
t188 = sin(qJ(4));
t193 = cos(qJ(4));
t145 = -t193 * mrSges(5,1) + mrSges(5,2) * t188;
t274 = t188 / 0.2e1;
t190 = sin(qJ(2));
t257 = qJD(1) * pkin(1);
t233 = t190 * t257;
t195 = cos(qJ(2));
t272 = pkin(1) * t195;
t126 = -qJD(2) * t233 + qJDD(1) * t272;
t180 = qJDD(1) + qJDD(2);
t287 = pkin(2) * t180 - qJD(3) * t233 + t126;
t185 = qJ(4) + qJ(5);
t175 = cos(t185);
t160 = t175 * mrSges(6,1);
t269 = pkin(4) * t193;
t164 = pkin(3) + t269;
t286 = m(5) * pkin(3) + m(6) * t164 + mrSges(4,1) - t145 + t160;
t219 = mrSges(5,1) * t188 + mrSges(5,2) * t193;
t182 = qJD(1) + qJD(2);
t172 = qJD(3) + t182;
t232 = t195 * t257;
t140 = pkin(2) * t182 + t232;
t189 = sin(qJ(3));
t151 = t189 * t233;
t194 = cos(qJ(3));
t98 = t140 * t194 - t151;
t83 = -pkin(3) * t172 - t98;
t285 = t83 * t219 + qJD(4) * (Ifges(5,5) * t193 - Ifges(5,6) * t188) / 0.2e1;
t197 = -pkin(9) - pkin(8);
t284 = -m(5) * pkin(8) + m(6) * t197 + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t187 = sin(qJ(5));
t192 = cos(qJ(5));
t132 = t187 * t193 + t188 * t192;
t230 = qJD(4) * t197;
t137 = t188 * t230;
t138 = t193 * t230;
t147 = t197 * t188;
t178 = t193 * pkin(9);
t148 = pkin(8) * t193 + t178;
t92 = t147 * t187 + t148 * t192;
t281 = -qJD(5) * t92 + t132 * t98 - t137 * t187 + t138 * t192;
t131 = -t187 * t188 + t192 * t193;
t91 = t147 * t192 - t148 * t187;
t280 = qJD(5) * t91 - t131 * t98 + t137 * t192 + t138 * t187;
t115 = t194 * t232 - t151;
t271 = pkin(2) * t189;
t163 = pkin(8) + t271;
t265 = -pkin(9) - t163;
t124 = t265 * t188;
t246 = t163 * t193;
t125 = t178 + t246;
t77 = t124 * t187 + t125 * t192;
t226 = qJD(4) * t265;
t237 = qJD(3) * t194;
t231 = pkin(2) * t237;
t93 = t188 * t226 + t193 * t231;
t94 = -t188 * t231 + t193 * t226;
t279 = -qJD(5) * t77 + t132 * t115 - t187 * t93 + t192 * t94;
t76 = t124 * t192 - t125 * t187;
t278 = qJD(5) * t76 - t131 * t115 + t187 * t94 + t192 * t93;
t99 = t140 * t189 + t194 * t233;
t84 = pkin(8) * t172 + t99;
t225 = (t188 ^ 2 + t193 ^ 2) * t84;
t277 = m(4) / 0.2e1;
t103 = t132 * t172;
t275 = t103 / 0.2e1;
t270 = pkin(2) * t194;
t186 = qJ(1) + qJ(2);
t177 = qJ(3) + t186;
t162 = cos(t177);
t268 = g(3) * t162;
t235 = qJD(4) * t193;
t171 = qJDD(3) + t180;
t239 = qJD(2) * t195;
t127 = (qJD(1) * t239 + qJDD(1) * t190) * pkin(1);
t50 = t194 * t127 + t140 * t237 + t189 * t287;
t40 = pkin(8) * t171 + t50;
t17 = -t188 * t40 - t235 * t84;
t267 = t17 * mrSges(5,3);
t166 = pkin(2) + t272;
t242 = t190 * t194;
t119 = pkin(1) * t242 + t166 * t189;
t113 = pkin(8) + t119;
t266 = -pkin(9) - t113;
t173 = sin(t185);
t263 = mrSges(6,2) * t173;
t262 = mrSges(6,2) * t175;
t102 = t131 * t172;
t261 = mrSges(6,3) * t102;
t260 = Ifges(5,4) * t188;
t259 = Ifges(5,4) * t193;
t258 = Ifges(6,4) * t103;
t256 = qJD(4) * pkin(4);
t255 = t103 * mrSges(6,3);
t236 = qJD(4) * t188;
t16 = t193 * t40 - t236 * t84;
t254 = t16 * t193;
t253 = t17 * t188;
t229 = pkin(9) * t172 + t84;
t69 = t229 * t193;
t252 = t187 * t69;
t251 = t192 * t69;
t249 = t193 * Ifges(5,2);
t161 = sin(t177);
t247 = t161 * t173;
t245 = t172 * t188;
t244 = t172 * t193;
t243 = t189 * t190;
t129 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t244;
t241 = t193 * t129;
t240 = mrSges(6,1) * t247 + t161 * t262;
t238 = qJD(3) * t189;
t234 = m(4) + m(5) + m(6);
t168 = pkin(4) * t236;
t227 = qJD(4) * t266;
t108 = t145 * t172;
t224 = -mrSges(4,1) * t172 + t108;
t223 = t160 - t263;
t118 = -pkin(1) * t243 + t166 * t194;
t64 = -mrSges(6,1) * t102 + mrSges(6,2) * t103;
t221 = -t224 - t64;
t112 = -pkin(3) - t118;
t68 = t229 * t188;
t220 = pkin(2) * t234 + mrSges(3,1);
t218 = -mrSges(6,1) * t173 - t262;
t217 = t249 + t260;
t215 = -t253 + t254;
t67 = -t68 + t256;
t21 = t192 * t67 - t252;
t22 = t187 * t67 + t251;
t89 = t266 * t188;
t90 = t113 * t193 + t178;
t56 = -t187 * t90 + t192 * t89;
t57 = t187 * t89 + t192 * t90;
t128 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t245;
t214 = t128 * t193 + t129 * t188;
t213 = t128 * t188 - t241;
t212 = t189 * t195 + t242;
t211 = mrSges(2,1) + (m(3) + t234) * pkin(1);
t51 = -t189 * t127 - t140 * t238 + t194 * t287;
t209 = t188 * (Ifges(5,1) * t193 - t260);
t208 = t132 * qJD(5);
t207 = t131 * qJD(5);
t41 = -pkin(3) * t171 - t51;
t109 = t171 * t193 - t172 * t236;
t95 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t109;
t110 = t171 * t188 + t172 * t235;
t96 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t110;
t205 = -qJD(4) * t214 - t188 * t96 + t193 * t95;
t179 = qJDD(4) + qJDD(5);
t181 = qJD(4) + qJD(5);
t10 = qJDD(4) * pkin(4) - pkin(9) * t110 + t17;
t11 = pkin(9) * t109 + t16;
t3 = qJD(5) * t21 + t10 * t187 + t11 * t192;
t4 = -qJD(5) * t22 + t10 * t192 - t11 * t187;
t44 = t109 * t187 + t110 * t192 + t172 * t207;
t45 = t109 * t192 - t110 * t187 - t172 * t208;
t59 = Ifges(6,2) * t102 + Ifges(6,6) * t181 + t258;
t97 = Ifges(6,4) * t102;
t60 = Ifges(6,1) * t103 + Ifges(6,5) * t181 + t97;
t73 = -t164 * t172 - t98;
t204 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t21 * t261 + t59 * t275 - t73 * (mrSges(6,1) * t103 + mrSges(6,2) * t102) + Ifges(6,3) * t179 - t103 * (Ifges(6,1) * t102 - t258) / 0.2e1 + Ifges(6,6) * t45 + Ifges(6,5) * t44 - t181 * (Ifges(6,5) * t102 - Ifges(6,6) * t103) / 0.2e1 - (-Ifges(6,2) * t103 + t60 + t97) * t102 / 0.2e1;
t75 = t166 * t238 + (qJD(2) * t212 + t190 * t237) * pkin(1);
t203 = (-t263 + t286) * t162 - t284 * t161;
t202 = -mrSges(6,2) * t247 + t161 * t286 + t162 * t284;
t174 = sin(t186);
t176 = cos(t186);
t201 = -mrSges(3,2) * t174 + t176 * t220 + t203;
t100 = Ifges(5,6) * qJD(4) + t172 * t217;
t146 = Ifges(5,4) * t244;
t101 = Ifges(5,1) * t245 + Ifges(5,5) * qJD(4) + t146;
t27 = -pkin(4) * t109 + t41;
t81 = qJD(4) * t131 + t207;
t82 = -qJD(4) * t132 - t208;
t200 = t193 * (Ifges(5,4) * t110 + Ifges(5,2) * t109) / 0.2e1 + t181 * (Ifges(6,5) * t81 + Ifges(6,6) * t82) / 0.2e1 + Ifges(4,3) * t171 + t41 * t145 + t102 * (Ifges(6,4) * t81 + Ifges(6,2) * t82) / 0.2e1 + t82 * t59 / 0.2e1 + t73 * (-mrSges(6,1) * t82 + mrSges(6,2) * t81) + t81 * t60 / 0.2e1 + t51 * mrSges(4,1) - t100 * t236 / 0.2e1 + t109 * t217 / 0.2e1 + t110 * (Ifges(5,1) * t188 + t259) / 0.2e1 + mrSges(5,3) * t254 + (Ifges(5,1) * t110 + Ifges(5,4) * t109) * t274 + (Ifges(6,1) * t81 + Ifges(6,4) * t82) * t275 + (t101 + t172 * (-Ifges(5,2) * t188 + t259)) * t235 / 0.2e1 + (-t21 * t81 + t22 * t82) * mrSges(6,3) + (0.2e1 * Ifges(5,5) * t274 + Ifges(5,6) * t193) * qJDD(4) + (t209 * t172 / 0.2e1 + t285) * qJD(4) + (mrSges(6,2) * t27 - mrSges(6,3) * t4 + Ifges(6,1) * t44 + Ifges(6,4) * t45 + Ifges(6,5) * t179) * t132 + (-mrSges(6,1) * t27 + mrSges(6,3) * t3 + Ifges(6,4) * t44 + Ifges(6,2) * t45 + Ifges(6,6) * t179) * t131;
t199 = mrSges(3,2) * t176 + t174 * t220 + t202;
t198 = mrSges(3,1) * t126 + Ifges(3,3) * t180 + t200;
t196 = cos(qJ(1));
t191 = sin(qJ(1));
t165 = -pkin(3) - t270;
t143 = -t164 - t270;
t139 = pkin(2) * t238 + t168;
t114 = t212 * t257;
t107 = t112 - t269;
t86 = mrSges(6,1) * t181 - t255;
t85 = -mrSges(6,2) * t181 + t261;
t74 = t166 * t237 + (-t190 * t238 + (t194 * t195 - t243) * qJD(2)) * pkin(1);
t72 = t168 + t75;
t71 = -mrSges(5,1) * t109 + mrSges(5,2) * t110;
t48 = -t188 * t74 + t193 * t227;
t47 = t188 * t227 + t193 * t74;
t36 = -mrSges(6,2) * t179 + mrSges(6,3) * t45;
t35 = mrSges(6,1) * t179 - mrSges(6,3) * t44;
t26 = -t192 * t68 - t252;
t25 = t187 * t68 - t251;
t9 = -mrSges(6,1) * t45 + mrSges(6,2) * t44;
t6 = -qJD(5) * t57 - t187 * t47 + t192 * t48;
t5 = qJD(5) * t56 + t187 * t48 + t192 * t47;
t1 = [t74 * t241 + (mrSges(2,2) * t196 + t191 * t211 + t199) * g(3) + t205 * t113 - t127 * mrSges(3,2) + t112 * t71 + t107 * t9 + t75 * t108 + t5 * t85 + t6 * t86 + (-t128 * t74 - t267) * t188 + m(4) * (t118 * t51 + t119 * t50 + t74 * t99 - t75 * t98) + m(6) * (t107 * t27 + t21 * t6 + t22 * t5 + t3 * t57 + t4 * t56 + t72 * t73) + t72 * t64 + t56 * t35 + t57 * t36 + (-mrSges(2,2) * t191 + t196 * t211 + t201) * g(2) + (-t119 * t171 - t172 * t74 - t50) * mrSges(4,2) + (t118 * t171 - t172 * t75) * mrSges(4,1) + m(5) * (t112 * t41 + t113 * t215 + t225 * t74 + t75 * t83) + t198 + Ifges(2,3) * qJDD(1) + ((-t180 * t190 - t182 * t239) * mrSges(3,2) + (-qJD(2) * t182 * t190 + t180 * t195) * mrSges(3,1) + m(3) * (t126 * t195 + t127 * t190)) * pkin(1); t165 * t71 + t143 * t9 + t139 * t64 + t171 * mrSges(4,1) * t270 - m(4) * (-t114 * t98 + t115 * t99) + t76 * t35 + t77 * t36 + (t182 * t232 - t127) * mrSges(3,2) + 0.2e1 * ((t189 * t50 + t194 * t51) * t277 + (m(5) * (t189 * t83 + t194 * t225) / 0.2e1 + (-t189 * t98 + t194 * t99) * t277) * qJD(3)) * pkin(2) - m(5) * (t114 * t83 + t115 * t225) + (t224 * t189 + (-mrSges(4,2) * t172 - t213) * t194) * pkin(2) * qJD(3) + (t115 * t172 - t171 * t271 - t50) * mrSges(4,2) + t279 * t86 + t221 * t114 + t213 * t115 + t205 * t163 + t199 * g(3) + t201 * g(2) + t278 * t85 + t182 * mrSges(3,1) * t233 + m(5) * (t16 * t246 - t163 * t253 + t165 * t41) - mrSges(5,3) * t253 + t198 + (t143 * t27 + t3 * t77 + t4 * t76 + (t139 - t114) * t73 + t278 * t22 + t279 * t21) * m(6); -t164 * t9 + t91 * t35 + t92 * t36 + (t64 * t256 - t267 + t98 * t128 + (-qJD(4) * t129 - t96) * pkin(8)) * t188 + t203 * g(2) + t202 * g(3) + (-t98 * t129 + (-qJD(4) * t128 + t95) * pkin(8)) * t193 + (t172 * t98 - t50) * mrSges(4,2) - pkin(3) * t71 + t281 * t86 + t280 * t85 + t221 * t99 + t200 + (-t164 * t27 + t3 * t92 + t4 * t91 + (t168 - t99) * t73 + t280 * t22 + t281 * t21) * m(6) + (-pkin(3) * t41 + pkin(8) * t215 - t225 * t98 - t83 * t99) * m(5); Ifges(5,5) * t110 + Ifges(5,6) * t109 - t26 * t85 - t25 * t86 + (-t64 * t245 + t187 * t36 + t192 * t35 + (t187 * t3 + t192 * t4 - g(1) * t193 + (-g(2) * t161 - t172 * t73 + t268) * t188) * m(6) + (-t187 * t86 + t192 * t85 + (-t187 * t21 + t192 * t22) * m(6)) * qJD(5)) * pkin(4) - m(6) * (t21 * t25 + t22 * t26) + (t100 * t274 + (t249 * t274 - t209 / 0.2e1) * t172 - (t146 + t101) * t193 / 0.2e1 - t285) * t172 - t16 * mrSges(5,2) + t17 * mrSges(5,1) + t22 * t255 + (-t161 * t219 - t240) * g(2) + t204 + Ifges(5,3) * qJDD(4) + (-t218 + t219) * t268 + t214 * t84 + (t145 - t223) * g(1); -t218 * t268 - g(1) * t223 - g(2) * t240 - t21 * t85 + t204 + (t86 + t255) * t22;];
tau = t1;
