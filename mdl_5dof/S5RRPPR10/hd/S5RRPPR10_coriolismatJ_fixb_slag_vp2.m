% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:05
% EndTime: 2019-12-31 19:43:12
% DurationCPUTime: 2.88s
% Computational Cost: add. (4778->374), mult. (10385->533), div. (0->0), fcn. (9702->6), ass. (0->181)
t197 = cos(pkin(8));
t196 = sin(pkin(8));
t220 = qJ(4) * t196 + pkin(2);
t164 = -pkin(3) * t197 - t220;
t168 = -mrSges(5,1) * t197 - mrSges(5,3) * t196;
t292 = m(5) * t164 + t168;
t198 = sin(qJ(5));
t267 = cos(qJ(5));
t208 = t198 * t196 + t197 * t267;
t291 = t208 / 0.2e1;
t290 = Ifges(5,4) + Ifges(4,5);
t289 = Ifges(4,6) - Ifges(5,6);
t278 = pkin(3) + pkin(4);
t288 = t196 * t278;
t154 = t267 * t196 - t198 * t197;
t199 = sin(qJ(2));
t200 = cos(qJ(2));
t246 = qJ(3) * t200;
t170 = t199 * pkin(2) - t246;
t238 = t197 * t170;
t240 = t196 * t199;
t118 = pkin(6) * t240 + t238;
t237 = t197 * t199;
t119 = -pkin(6) * t237 + t196 * t170;
t287 = -t118 * t196 + t119 * t197;
t226 = -pkin(6) * t196 - pkin(3);
t102 = t199 * t226 - t238;
t97 = t199 * qJ(4) + t119;
t286 = t102 * t196 + t197 * t97;
t285 = mrSges(5,3) / 0.2e1 - mrSges(4,2) / 0.2e1;
t284 = m(5) / 0.2e1;
t283 = m(6) / 0.2e1;
t282 = m(4) * pkin(6);
t281 = -mrSges(6,3) / 0.2e1;
t264 = -pkin(7) + qJ(3);
t167 = t264 * t196;
t169 = t264 * t197;
t91 = t167 * t267 - t198 * t169;
t280 = t91 / 0.2e1;
t92 = t198 * t167 + t169 * t267;
t279 = -t92 / 0.2e1;
t133 = t208 * t199;
t277 = -t133 / 0.2e1;
t134 = t154 * t200;
t276 = t134 / 0.2e1;
t135 = t208 * t200;
t275 = t135 / 0.2e1;
t274 = -t208 / 0.2e1;
t273 = t154 / 0.2e1;
t272 = -t196 / 0.2e1;
t271 = t196 / 0.2e1;
t270 = t197 / 0.2e1;
t269 = t199 / 0.2e1;
t268 = t200 / 0.2e1;
t265 = pkin(3) * t196;
t263 = Ifges(4,4) * t196;
t262 = Ifges(4,4) * t197;
t261 = Ifges(6,4) * t133;
t260 = Ifges(6,4) * t154;
t259 = Ifges(5,5) * t196;
t258 = Ifges(5,5) * t197;
t166 = -pkin(2) * t200 - t199 * qJ(3) - pkin(1);
t239 = t196 * t200;
t181 = pkin(6) * t239;
t106 = t197 * t166 - t181;
t236 = t197 * t200;
t107 = pkin(6) * t236 + t196 * t166;
t132 = t154 * t199;
t156 = t200 * mrSges(4,2) - mrSges(4,3) * t240;
t158 = -t200 * mrSges(4,1) - mrSges(4,3) * t237;
t159 = t200 * mrSges(5,1) + mrSges(5,2) * t237;
t163 = -mrSges(5,2) * t240 - t200 * mrSges(5,3);
t193 = t200 * pkin(3);
t72 = pkin(4) * t200 + t181 + t193 + (-pkin(7) * t199 - t166) * t197;
t95 = -qJ(4) * t200 + t107;
t77 = pkin(7) * t240 + t95;
t29 = -t198 * t77 + t267 * t72;
t30 = t198 * t72 + t267 * t77;
t96 = -t106 + t193;
t98 = -mrSges(6,2) * t200 + t132 * mrSges(6,3);
t248 = t200 * mrSges(6,1);
t255 = t133 * mrSges(6,3);
t99 = t248 - t255;
t9 = m(6) * (-t132 * t30 + t133 * t29) - t132 * t98 + t133 * t99 + ((-t158 + t159) * t197 + (-t156 - t163) * t196 + m(5) * (-t196 * t95 + t197 * t96) + m(4) * (-t106 * t197 - t107 * t196)) * t199;
t257 = qJD(1) * t9;
t100 = mrSges(6,2) * t199 + mrSges(6,3) * t134;
t101 = -mrSges(6,1) * t199 - mrSges(6,3) * t135;
t124 = Ifges(5,6) * t199 + (t196 * Ifges(5,3) + t258) * t200;
t125 = Ifges(4,6) * t199 + (-t196 * Ifges(4,2) + t262) * t200;
t126 = Ifges(5,4) * t199 + (t197 * Ifges(5,1) + t259) * t200;
t127 = Ifges(4,5) * t199 + (t197 * Ifges(4,1) - t263) * t200;
t219 = qJ(4) * t197 - pkin(6);
t128 = (-t219 + t265) * t199;
t179 = qJ(4) * t236;
t129 = -t179 + (pkin(6) + t265) * t200;
t216 = t196 * mrSges(5,1) - t197 * mrSges(5,3);
t146 = t216 * t199;
t147 = t216 * t200;
t148 = (t196 * mrSges(4,1) + t197 * mrSges(4,2)) * t200;
t157 = -t199 * mrSges(4,2) - mrSges(4,3) * t239;
t160 = t199 * mrSges(4,1) - mrSges(4,3) * t236;
t180 = mrSges(5,2) * t236;
t249 = t199 * mrSges(5,1);
t161 = t180 - t249;
t162 = -mrSges(5,2) * t239 + t199 * mrSges(5,3);
t209 = Ifges(6,5) * t275 + Ifges(6,6) * t276;
t73 = (-pkin(7) * t200 - t170) * t197 + (-pkin(4) + t226) * t199;
t78 = pkin(7) * t239 + t97;
t31 = -t198 * t78 + t267 * t73;
t32 = t198 * t73 + t267 * t78;
t58 = Ifges(6,2) * t132 + t200 * Ifges(6,6) + t261;
t59 = Ifges(6,4) * t135 + Ifges(6,2) * t134 - Ifges(6,6) * t199;
t123 = Ifges(6,4) * t132;
t60 = Ifges(6,1) * t133 + t200 * Ifges(6,5) + t123;
t61 = Ifges(6,1) * t135 + Ifges(6,4) * t134 - Ifges(6,5) * t199;
t64 = -mrSges(6,1) * t132 + mrSges(6,2) * t133;
t253 = t135 * mrSges(6,2);
t254 = t134 * mrSges(6,1);
t65 = t253 - t254;
t93 = (t219 - t288) * t199;
t94 = -t179 + (pkin(6) + t288) * t200;
t1 = (-pkin(1) * mrSges(3,2) + (Ifges(3,1) - Ifges(3,2) - Ifges(5,2) - Ifges(4,3) - Ifges(6,3) + m(4) * pkin(6) ^ 2 + (pkin(6) * mrSges(4,2) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t197) * t197 + (pkin(6) * mrSges(4,1) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t196 + (-Ifges(4,4) + Ifges(5,5)) * t197) * t196) * t199 + t209 + (t289 * t196 - t290 * t197 + Ifges(3,4)) * t200) * t200 + (Ifges(6,5) * t277 - Ifges(6,6) * t132 / 0.2e1 - pkin(1) * mrSges(3,1) + pkin(6) * t148 - Ifges(3,4) * t199 + (t126 / 0.2e1 + t127 / 0.2e1 + (Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t199) * t197 + (t124 / 0.2e1 - t125 / 0.2e1 + (Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * t199) * t196) * t199 + m(4) * (t106 * t118 + t107 * t119) + t102 * t159 + t106 * t160 + t96 * t161 + t95 * t162 + t97 * t163 + t119 * t156 + t107 * t157 + t118 * t158 + t129 * t146 + t128 * t147 + t132 * t59 / 0.2e1 + t133 * t61 / 0.2e1 + t29 * t101 + t93 * t65 - t94 * t64 + t32 * t98 + t31 * t99 + t30 * t100 + m(6) * (t29 * t31 + t30 * t32 - t93 * t94) + m(5) * (t102 * t96 + t128 * t129 + t95 * t97) + t60 * t275 + t58 * t276;
t256 = t1 * qJD(1);
t251 = t198 * mrSges(6,1);
t250 = t198 * t99;
t215 = t133 * mrSges(6,1) + t132 * mrSges(6,2);
t232 = Ifges(6,5) * t132 - Ifges(6,6) * t133;
t66 = -Ifges(6,2) * t133 + t123;
t67 = Ifges(6,1) * t132 - t261;
t4 = t29 * t98 - t30 * t99 + t93 * t215 + t232 * t268 + (-t30 * mrSges(6,3) + t67 / 0.2e1 - t58 / 0.2e1) * t133 - (t29 * mrSges(6,3) - t60 / 0.2e1 - t66 / 0.2e1) * t132;
t247 = t4 * qJD(1);
t228 = t267 * t98;
t12 = (-m(5) * t128 + m(6) * t93 - t146 + t64) * t237 + (m(6) * (t198 * t29 - t267 * t30) - t228 + t250 - t163 - m(5) * t95) * t200;
t245 = qJD(1) * t12;
t218 = t267 * t281;
t207 = -t228 / 0.2e1 - t132 * t218;
t229 = t267 * mrSges(6,2);
t217 = t229 / 0.2e1;
t13 = t200 * t217 + (t248 / 0.2e1 + t255 / 0.2e1 + t99 / 0.2e1) * t198 + t207;
t241 = t13 * qJD(1);
t235 = t198 * t132;
t234 = t198 * t208;
t231 = -Ifges(6,5) * t208 - Ifges(6,6) * t154;
t227 = t198 * t281;
t225 = t267 * t133;
t224 = t267 * t154;
t221 = -t237 / 0.2e1;
t80 = t154 * mrSges(6,1) - mrSges(6,2) * t208;
t145 = t197 * t278 + t220;
t151 = Ifges(6,4) * t208;
t82 = -Ifges(6,2) * t154 - t151;
t83 = -Ifges(6,2) * t208 + t260;
t84 = -Ifges(6,1) * t208 - t260;
t85 = Ifges(6,1) * t154 - t151;
t11 = t145 * t80 + (-t83 / 0.2e1 + t84 / 0.2e1) * t154 - (t82 / 0.2e1 + t85 / 0.2e1) * t208;
t203 = -(t60 / 0.4e1 + t66 / 0.4e1) * t208 + (t67 / 0.4e1 - t58 / 0.4e1) * t154 - (mrSges(6,3) * t280 - t85 / 0.4e1 - t82 / 0.4e1) * t132 + (mrSges(6,3) * t279 + t84 / 0.4e1 - t83 / 0.4e1) * t133 + t145 * t215 / 0.2e1 + t200 * t231 / 0.4e1 + t98 * t280 + t99 * t279 + t93 * t80 / 0.2e1;
t205 = Ifges(6,3) * t269 - t31 * mrSges(6,1) / 0.2e1 + t32 * mrSges(6,2) / 0.2e1 - t209;
t3 = t203 + t205;
t214 = t3 * qJD(1) + t11 * qJD(2);
t16 = (-t154 ^ 2 - t208 ^ 2) * mrSges(6,3) + m(6) * (t154 * t91 + t208 * t92) + (mrSges(4,3) + mrSges(5,2) + 0.4e1 * (m(4) / 0.4e1 + m(5) / 0.4e1) * qJ(3)) * (t196 ^ 2 + t197 ^ 2);
t202 = (-t158 / 0.2e1 + t159 / 0.2e1) * t196 + (t163 / 0.2e1 + t156 / 0.2e1) * t197 + (-t132 * t274 + t154 * t277) * mrSges(6,3) + m(4) * (-t106 * t196 + t107 * t197) / 0.2e1 + (t196 * t96 + t197 * t95) * t284 + (-t132 * t92 + t133 * t91 + t154 * t29 + t208 * t30) * t283 + t98 * t291 + t99 * t273;
t206 = -m(5) * t129 / 0.2e1 - m(6) * t94 / 0.2e1 - t254 / 0.2e1 + t253 / 0.2e1;
t6 = (-t282 / 0.2e1 + t285 * t197 + (-mrSges(5,1) / 0.2e1 - mrSges(4,1) / 0.2e1) * t196) * t200 + t202 + t206;
t213 = qJD(1) * t6 + qJD(2) * t16;
t81 = mrSges(6,1) * t208 + mrSges(6,2) * t154;
t24 = (m(6) * t145 - t292 + t81) * t196;
t201 = (-t196 * t128 + (-t164 * t199 - t246) * t197) * t284 + (t145 * t237 + t196 * t93) * t283 + t146 * t272 + t64 * t271 + t168 * t221 + t81 * t237 / 0.2e1 + ((t198 * t91 - t267 * t92) * t283 + t154 * t227 - t208 * t218) * t200;
t204 = t102 * t284 + (t198 * t32 + t267 * t31) * t283 + t198 * t100 / 0.2e1 - t249 / 0.2e1 + t267 * t101 / 0.2e1;
t7 = -t180 + t201 - t204;
t212 = -qJD(1) * t7 - qJD(2) * t24;
t211 = qJD(1) * t215 + qJD(2) * t80;
t35 = -m(5) * t237 + (t221 + t235 / 0.2e1 - t225 / 0.2e1) * m(6);
t56 = m(5) * t196 + (t234 / 0.2e1 + t224 / 0.2e1 + t271) * m(6);
t210 = qJD(1) * t35 - qJD(2) * t56;
t62 = m(6) * t272 + (t224 + t234) * t283;
t42 = m(6) * t221 + (t225 - t235) * t283;
t14 = t133 * t227 - t250 / 0.2e1 + (t251 / 0.2e1 + t217) * t200 - t207;
t8 = t201 + t204;
t5 = t268 * t282 + t202 - t206 - t285 * t236 + (mrSges(4,1) + mrSges(5,1)) * t239 / 0.2e1;
t2 = t203 - t205;
t10 = [qJD(2) * t1 + qJD(3) * t9 + qJD(4) * t12 + qJD(5) * t4, t5 * qJD(3) + t8 * qJD(4) + t2 * qJD(5) + t256 + (-t197 * t124 / 0.2e1 + t164 * t147 + t145 * t65 - pkin(2) * t148 + t91 * t101 - t94 * t81 + t92 * t100 + t125 * t270 + t61 * t273 + t59 * t274 + t85 * t275 + t83 * t276 + (t127 + t126) * t271 + (t290 * t196 + t289 * t197) * t269 + (-Ifges(3,6) - Ifges(6,5) * t154 / 0.2e1 + Ifges(6,6) * t291 + pkin(6) * mrSges(3,2)) * t199 + (-t31 * t154 - t208 * t32) * mrSges(6,3) + t287 * mrSges(4,3) + t286 * mrSges(5,2) + 0.2e1 * (-t145 * t94 + t31 * t91 + t32 * t92) * t283 + ((t157 + t162) * t197 + (-t160 + t161) * t196 + m(5) * t286 + m(4) * t287) * qJ(3) + (Ifges(3,5) + (Ifges(4,2) * t197 + t263) * t272 + (-Ifges(5,3) * t197 + t259) * t271 + (-m(4) * pkin(2) - mrSges(4,1) * t197 + mrSges(4,2) * t196 - mrSges(3,1)) * pkin(6) + (-t258 + t262 + (Ifges(4,1) + Ifges(5,1)) * t196) * t270) * t200 + t292 * t129) * qJD(2), qJD(2) * t5 + qJD(4) * t42 + t257, qJD(2) * t8 + qJD(3) * t42 + qJD(5) * t14 + t245, t247 + t2 * qJD(2) + t14 * qJD(4) + (-mrSges(6,1) * t30 - mrSges(6,2) * t29 + t232) * qJD(5); qJD(3) * t6 + qJD(4) * t7 + qJD(5) * t3 - t256, qJD(3) * t16 + qJD(4) * t24 + qJD(5) * t11, qJD(4) * t62 + t213, qJD(3) * t62 - t212, (-mrSges(6,1) * t92 - mrSges(6,2) * t91 + t231) * qJD(5) + t214; -qJD(2) * t6 + qJD(4) * t35 - qJD(5) * t215 - t257, -qJD(4) * t56 - qJD(5) * t80 - t213, 0, t210, -t211; -qJD(2) * t7 - qJD(3) * t35 - qJD(5) * t13 - t245, qJD(3) * t56 + t212, -t210, 0, -t241 + (-t229 - t251) * qJD(5); -qJD(2) * t3 + qJD(3) * t215 + qJD(4) * t13 - t247, qJD(3) * t80 - t214, t211, t241, 0;];
Cq = t10;
