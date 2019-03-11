% Calculate vector of inverse dynamics joint torques for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPPRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:44
% EndTime: 2019-03-08 19:17:58
% DurationCPUTime: 7.85s
% Computational Cost: add. (3371->446), mult. (7400->625), div. (0->0), fcn. (5601->12), ass. (0->221)
t306 = m(6) + m(7);
t124 = sin(pkin(11));
t126 = sin(pkin(6));
t132 = sin(qJ(2));
t232 = t126 * t132;
t196 = qJD(1) * t232;
t107 = t124 * t196;
t127 = cos(pkin(11));
t135 = cos(qJ(2));
t220 = qJD(1) * t135;
t195 = t126 * t220;
t79 = t127 * t195 - t107;
t302 = -t79 + qJD(4);
t305 = m(6) + m(5);
t131 = sin(qJ(5));
t134 = cos(qJ(5));
t181 = pkin(5) * t134 + pkin(9) * t131;
t304 = qJD(5) * t181 + t302;
t130 = sin(qJ(6));
t133 = cos(qJ(6));
t216 = qJD(5) * t133;
t219 = qJD(2) * t134;
t98 = -t130 * t219 + t216;
t301 = -t98 / 0.2e1;
t218 = qJD(5) * t130;
t99 = t133 * t219 + t218;
t300 = -t99 / 0.2e1;
t210 = t131 * qJD(2);
t120 = qJD(6) + t210;
t299 = -t120 / 0.2e1;
t125 = sin(pkin(10));
t128 = cos(pkin(10));
t129 = cos(pkin(6));
t227 = t129 * t135;
t156 = -t125 * t227 - t128 * t132;
t146 = t156 * pkin(2);
t282 = -m(7) - t305;
t202 = mrSges(6,3) * t219;
t298 = qJD(5) * mrSges(6,1) + mrSges(7,1) * t98 - mrSges(7,2) * t99 - t202;
t174 = -mrSges(7,1) * t133 + mrSges(7,2) * t130;
t149 = m(7) * pkin(5) - t174;
t175 = mrSges(6,1) * t131 + mrSges(6,2) * t134;
t186 = pkin(9) * t134 - qJ(4);
t254 = -mrSges(4,2) + mrSges(5,3);
t297 = m(7) * t186 + t134 * mrSges(7,3) - qJ(4) * t305 - t131 * t149 - t175 - t254;
t117 = qJD(1) * t129 + qJD(3);
t104 = qJD(2) * pkin(2) + t195;
t68 = t104 * t127 - t107;
t176 = qJD(4) - t68;
t266 = -pkin(3) - pkin(8);
t59 = qJD(2) * t266 + t176;
t38 = t117 * t134 + t131 * t59;
t36 = qJD(5) * pkin(9) + t38;
t158 = pkin(5) * t131 - t186;
t69 = t104 * t124 + t127 * t196;
t44 = qJD(2) * t158 + t69;
t10 = t130 * t44 + t133 * t36;
t208 = qJD(2) * qJD(5);
t102 = qJDD(2) * t134 - t131 * t208;
t103 = -qJDD(2) * t131 - t134 * t208;
t209 = qJD(2) * qJD(4);
t236 = qJDD(2) * pkin(2);
t194 = qJD(2) * t232;
t183 = qJD(1) * t194;
t230 = t126 * t135;
t85 = qJDD(1) * t230 - t183;
t81 = t85 + t236;
t190 = qJD(2) * t220;
t86 = (qJDD(1) * t132 + t190) * t126;
t40 = t124 * t81 + t127 * t86;
t31 = qJDD(2) * qJ(4) + t209 + t40;
t13 = -pkin(5) * t103 - pkin(9) * t102 + t31;
t116 = qJDD(1) * t129 + qJDD(3);
t215 = qJD(5) * t134;
t217 = qJD(5) * t131;
t39 = -t124 * t86 + t127 * t81;
t178 = qJDD(4) - t39;
t24 = qJDD(2) * t266 + t178;
t7 = t134 * t116 - t117 * t217 + t131 * t24 + t59 * t215;
t5 = qJDD(5) * pkin(9) + t7;
t9 = -t130 * t36 + t133 * t44;
t1 = qJD(6) * t9 + t13 * t130 + t133 * t5;
t2 = -qJD(6) * t10 + t13 * t133 - t130 * t5;
t180 = t1 * t133 - t130 * t2;
t213 = qJD(6) * t133;
t214 = qJD(6) * t130;
t296 = -t10 * t214 - t9 * t213 + t180;
t295 = -t125 * t132 + t128 * t227;
t57 = qJD(6) * t98 + qJDD(5) * t130 + t102 * t133;
t272 = t57 / 0.2e1;
t58 = -qJD(6) * t99 + qJDD(5) * t133 - t102 * t130;
t271 = t58 / 0.2e1;
t94 = qJDD(6) - t103;
t270 = t94 / 0.2e1;
t294 = t102 / 0.2e1;
t293 = t103 / 0.2e1;
t255 = -mrSges(4,1) + mrSges(5,2);
t14 = -mrSges(7,1) * t58 + mrSges(7,2) * t57;
t253 = -qJDD(5) * mrSges(6,1) + mrSges(6,3) * t102 + t14;
t123 = -pkin(2) * t127 - pkin(3);
t119 = -pkin(8) + t123;
t192 = t119 * t215;
t226 = t130 * t131;
t224 = t131 * t133;
t265 = pkin(2) * t124;
t92 = t158 + t265;
t61 = t119 * t224 + t130 * t92;
t160 = t124 * t135 + t127 * t132;
t84 = t160 * t126;
t76 = qJD(1) * t84;
t292 = -qJD(6) * t61 - t130 * t192 + t133 * t304 + t226 * t76;
t60 = -t119 * t226 + t133 * t92;
t291 = qJD(6) * t60 + t130 * t304 + t133 * t192 - t224 * t76;
t290 = mrSges(6,1) + t149;
t289 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t252 = Ifges(6,4) * t131;
t172 = t134 * Ifges(6,1) - t252;
t93 = Ifges(7,4) * t98;
t54 = t99 * Ifges(7,1) + t120 * Ifges(7,5) + t93;
t288 = Ifges(6,5) * qJD(5) + qJD(2) * t172 + t133 * t54;
t237 = qJD(5) * t38;
t8 = -t116 * t131 + t134 * t24 - t237;
t121 = qJ(4) + t265;
t63 = qJD(2) * qJ(4) + t69;
t287 = t121 * t31 + t302 * t63;
t32 = mrSges(7,1) * t94 - mrSges(7,3) * t57;
t33 = -mrSges(7,2) * t94 + mrSges(7,3) * t58;
t284 = -t130 * t32 + t133 * t33;
t283 = t131 * t7 + t134 * t8;
t251 = Ifges(6,4) * t134;
t169 = -t131 * Ifges(6,2) + t251;
t281 = Ifges(6,6) * qJD(5) / 0.2e1 + qJD(2) * t169 / 0.2e1 + Ifges(7,5) * t300 + Ifges(7,6) * t301 + Ifges(7,3) * t299;
t173 = t130 * mrSges(7,1) + t133 * mrSges(7,2);
t37 = -t117 * t131 + t134 * t59;
t35 = -qJD(5) * pkin(5) - t37;
t258 = t99 * Ifges(7,4);
t53 = t98 * Ifges(7,2) + t120 * Ifges(7,6) + t258;
t280 = t173 * t35 - t130 * t53 / 0.2e1;
t228 = t129 * t132;
t221 = -t124 * t227 - t127 * t228;
t222 = t135 * t127;
t96 = t124 * t132 - t222;
t50 = t125 * t221 - t128 * t96;
t279 = t125 * t96 + t128 * t221;
t278 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t276 = pkin(8) * t306 + mrSges(6,3) + t173 - t255;
t203 = mrSges(6,3) * t210;
t110 = -qJD(5) * mrSges(6,2) - t203;
t6 = -qJDD(5) * pkin(5) - t8;
t71 = -mrSges(7,2) * t120 + mrSges(7,3) * t98;
t72 = mrSges(7,1) * t120 - mrSges(7,3) * t99;
t274 = m(6) * (t8 + t237) + m(7) * (t10 * t216 - t218 * t9 - t6) - qJD(5) * (t130 * t72 - t133 * t71 - t110) - t253;
t136 = qJD(2) ^ 2;
t273 = Ifges(7,1) * t272 + Ifges(7,4) * t271 + Ifges(7,5) * t270;
t268 = t99 / 0.2e1;
t261 = t134 * t6;
t250 = Ifges(7,4) * t130;
t249 = Ifges(7,4) * t133;
t242 = t134 * t76;
t239 = qJD(2) * t76;
t238 = qJD(2) * t79;
t233 = t126 * t131;
t231 = t126 * t134;
t225 = t130 * t134;
t223 = t133 * t134;
t212 = qJD(6) * t134;
t211 = qJDD(2) * mrSges(5,2);
t207 = Ifges(7,5) * t57 + Ifges(7,6) * t58 + Ifges(7,3) * t94;
t118 = pkin(2) * t230;
t205 = m(4) - t282;
t199 = t126 * t222;
t191 = t130 * t217;
t184 = -t208 / 0.2e1;
t182 = t295 * pkin(2);
t78 = qJD(2) * t199 - t124 * t194;
t179 = t31 * t84 + t63 * t78;
t177 = t10 * t130 + t9 * t133;
t171 = Ifges(7,1) * t133 - t250;
t170 = Ifges(7,1) * t130 + t249;
t168 = -Ifges(7,2) * t130 + t249;
t167 = Ifges(7,2) * t133 + t250;
t166 = -Ifges(6,5) * t131 - Ifges(6,6) * t134;
t165 = Ifges(7,5) * t133 - Ifges(7,6) * t130;
t164 = Ifges(7,5) * t130 + Ifges(7,6) * t133;
t83 = t124 * t232 - t199;
t66 = t129 * t134 + t131 * t83;
t23 = t130 * t84 + t133 * t66;
t22 = -t130 * t66 + t133 * t84;
t163 = -t130 * t71 - t133 * t72;
t65 = t129 * t131 - t83 * t134;
t155 = t63 * (mrSges(6,1) * t134 - mrSges(6,2) * t131);
t154 = t131 * (-Ifges(6,2) * t134 - t252);
t153 = t134 * (-Ifges(6,1) * t131 - t251);
t152 = (mrSges(3,1) * t135 - mrSges(3,2) * t132) * t126;
t150 = t96 * t129;
t145 = t130 * t212 + t131 * t216;
t144 = -t133 * t212 + t191;
t141 = Ifges(7,5) * t134 - t131 * t171;
t140 = Ifges(7,6) * t134 - t131 * t168;
t139 = Ifges(7,3) * t134 - t131 * t165;
t88 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t103;
t137 = t88 + t163 * qJD(6) - t298 * qJD(5) + m(7) * (qJD(5) * t35 + t296) + m(6) * (-qJD(5) * t37 + t7) + t284;
t105 = t116 * t129;
t101 = t181 * qJD(2);
t100 = t175 * qJD(2);
t77 = qJD(2) * t84;
t67 = -mrSges(6,1) * t103 + mrSges(6,2) * t102;
t62 = -qJD(2) * pkin(3) + t176;
t49 = t125 * t150 - t128 * t160;
t46 = -t125 * t160 - t128 * t150;
t34 = -qJDD(2) * pkin(3) + t178;
t30 = t128 * t231 + t46 * t131;
t28 = t125 * t231 - t131 * t49;
t20 = qJD(5) * t66 - t77 * t134;
t19 = -qJD(5) * t65 + t131 * t77;
t16 = t101 * t130 + t133 * t37;
t15 = t101 * t133 - t130 * t37;
t11 = t57 * Ifges(7,4) + t58 * Ifges(7,2) + t94 * Ifges(7,6);
t4 = qJD(6) * t22 + t130 * t78 + t133 * t19;
t3 = -qJD(6) * t23 - t130 * t19 + t133 * t78;
t12 = [m(2) * qJDD(1) + t78 * t100 + t19 * t110 + t22 * t32 + t23 * t33 + t3 * t72 + t4 * t71 + t66 * t88 + t84 * t67 + t253 * t65 - t298 * t20 + (-mrSges(3,1) * t132 - mrSges(3,2) * t135) * t136 * t126 + (t254 * t78 + t255 * t77) * qJD(2) + (t254 * t84 + t255 * t83 + t152) * qJDD(2) + (-m(2) - m(3) - t205) * g(3) + m(3) * (qJDD(1) * t129 ^ 2 + (t132 * t86 + t135 * t85) * t126) + m(4) * (-t39 * t83 + t40 * t84 - t68 * t77 + t69 * t78 + t105) + m(5) * (t34 * t83 + t62 * t77 + t105 + t179) + m(6) * (t19 * t38 - t20 * t37 - t65 * t8 + t66 * t7 + t179) + m(7) * (t1 * t23 + t10 * t4 + t2 * t22 + t20 * t35 + t3 * t9 + t6 * t65); (t127 * t236 + t239 + t39) * mrSges(4,1) + qJDD(5) * (Ifges(6,5) * t134 - Ifges(6,6) * t131) + (t85 + t183) * mrSges(3,1) + t121 * t67 + t60 * t32 + t61 * t33 + (t68 * t76 - t69 * t79 + (t124 * t40 + t127 * t39) * pkin(2)) * m(4) - t253 * t119 * t134 + t35 * (-mrSges(7,1) * t144 - mrSges(7,2) * t145) + (-m(4) * t182 - t295 * mrSges(3,1) - (-t125 * t135 - t128 * t228) * mrSges(3,2) + t282 * (t46 * pkin(3) + t182) - t276 * t46 - t297 * t279) * g(2) + t298 * (-t119 * t217 - t242) + (-t1 * t225 + t10 * t144 + t145 * t9 - t2 * t223) * mrSges(7,3) + (t119 * t88 - t76 * t110 - Ifges(6,4) * t102 / 0.2e1 - Ifges(6,2) * t103 / 0.2e1 + t207 / 0.2e1 + Ifges(7,3) * t270 + Ifges(7,6) * t271 + Ifges(7,5) * t272 + t278) * t131 + t154 * t184 + (t9 * mrSges(7,1) - t10 * mrSges(7,2) - t281) * t215 + (-m(4) * t118 - t152 + t282 * (-t83 * pkin(3) + t118) + t297 * t84 + t276 * t83) * g(3) + t291 * t71 + (t1 * t61 + t2 * t60 + (t217 * t35 - t261) * t119 + t242 * t35 + t292 * t9 + t291 * t10) * m(7) + t292 * t72 + (Ifges(5,1) + Ifges(4,3) + Ifges(3,3)) * qJDD(2) + t123 * t211 + qJD(5) * t155 + t169 * t293 + t172 * t294 + (-m(4) * t146 - t156 * mrSges(3,1) - (t125 * t228 - t128 * t135) * mrSges(3,2) - t276 * t49 + t297 * t50 + (-m(5) - t306) * (t49 * pkin(3) + t146)) * g(1) - (t130 * t54 + t133 * t53) * t212 / 0.2e1 + (qJDD(2) * t121 + t209 - t238 + t31) * mrSges(5,3) - t288 * t217 / 0.2e1 + (-t215 * t38 + t217 * t37 - t283) * mrSges(6,3) + (((-t37 * t131 + t38 * t134) * qJD(5) + t283) * t119 - (t131 * t38 + t134 * t37) * t76 + t287) * m(6) + (t123 * t34 - t62 * t76 + t287) * m(5) + qJD(5) ^ 2 * t166 / 0.2e1 + t31 * t175 + t110 * t192 + (t126 * t190 - t86) * mrSges(3,2) + t302 * t100 + (t34 - t239) * mrSges(5,2) - t11 * t225 / 0.2e1 + t98 * (qJD(5) * t140 - t167 * t212) / 0.2e1 + t120 * (qJD(5) * t139 - t164 * t212) / 0.2e1 + t153 * t208 / 0.2e1 + (Ifges(6,1) * t294 + Ifges(6,4) * t293 + t165 * t270 + t168 * t271 + t171 * t272) * t134 + t173 * t261 + (qJD(5) * t141 - t170 * t212) * t268 + t223 * t273 + (-t124 * t236 + t238 - t40) * mrSges(4,2) + t53 * t191 / 0.2e1; 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t116 - t274 * t131 + t137 * t134 + (-t129 * g(3) + (-g(1) * t125 + g(2) * t128) * t126) * t205; t211 - t136 * mrSges(5,3) + t34 * m(5) + t274 * t134 + t137 * t131 + (-m(7) * t177 - t305 * t63 - t100 + t163) * qJD(2) - t282 * (g(1) * t49 + g(2) * t46 - g(3) * t83); t133 * t11 / 0.2e1 + (-t110 - t203) * t37 + Ifges(6,5) * t102 + Ifges(6,6) * t103 - t16 * t71 - t15 * t72 - pkin(5) * t14 - t7 * mrSges(6,2) + t8 * mrSges(6,1) + (-t10 * (-mrSges(7,2) * t134 + mrSges(7,3) * t226) - t9 * (mrSges(7,1) * t134 + mrSges(7,3) * t224) - t155) * qJD(2) + (-m(7) * t35 + t202 + t298) * t38 + (-pkin(5) * t6 - t10 * t16 - t15 * t9) * m(7) + t166 * t184 + t281 * t219 + t296 * mrSges(7,3) + (-t153 / 0.2e1 + t154 / 0.2e1) * t136 + (t120 * t165 + t168 * t98 + t171 * t99) * qJD(6) / 0.2e1 - (t120 * t139 + t140 * t98 + t141 * t99) * qJD(2) / 0.2e1 + (t289 * t66 + t290 * t65) * g(3) + (t289 * t28 - t290 * (-t125 * t233 - t134 * t49)) * g(1) + (-t289 * t30 - t290 * (t128 * t233 - t134 * t46)) * g(2) + (m(7) * (-qJD(6) * t177 + t180) - t71 * t214 - t72 * t213 + t284) * pkin(9) + t6 * t174 + t280 * qJD(6) + Ifges(6,3) * qJDD(5) + (t288 / 0.2e1 + t280) * t210 + t54 * t213 / 0.2e1 + t164 * t270 + t167 * t271 + t170 * t272 + t130 * t273; -t35 * (mrSges(7,1) * t99 + mrSges(7,2) * t98) + (Ifges(7,1) * t98 - t258) * t300 + t53 * t268 + (Ifges(7,5) * t98 - Ifges(7,6) * t99) * t299 - t9 * t71 + t10 * t72 - g(1) * ((-t130 * t28 + t133 * t50) * mrSges(7,1) + (-t130 * t50 - t133 * t28) * mrSges(7,2)) - g(2) * ((t130 * t30 - t133 * t279) * mrSges(7,1) + (t130 * t279 + t133 * t30) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t22 - mrSges(7,2) * t23) + (t10 * t99 + t9 * t98) * mrSges(7,3) + t207 + (-Ifges(7,2) * t99 + t54 + t93) * t301 + t278;];
tau  = t12;
