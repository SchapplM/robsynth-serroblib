% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:17
% EndTime: 2019-03-08 20:09:29
% DurationCPUTime: 5.49s
% Computational Cost: add. (5041->427), mult. (13262->564), div. (0->0), fcn. (9996->10), ass. (0->206)
t146 = sin(pkin(11));
t148 = cos(pkin(11));
t151 = sin(qJ(4));
t154 = cos(qJ(4));
t127 = t146 * t151 - t154 * t148;
t147 = sin(pkin(6));
t155 = cos(qJ(2));
t212 = t147 * t155;
t161 = t127 * t212;
t237 = pkin(8) + qJ(3);
t136 = t237 * t146;
t137 = t237 * t148;
t267 = -t154 * t136 - t137 * t151;
t278 = qJD(1) * t161 - t127 * qJD(3) + qJD(4) * t267;
t122 = t127 * qJD(4);
t128 = t146 * t154 + t148 * t151;
t123 = t128 * qJD(4);
t152 = sin(qJ(2));
t210 = qJD(1) * t147;
t199 = t152 * t210;
t288 = pkin(4) * t123 + pkin(9) * t122 - t199;
t286 = Ifges(6,1) + Ifges(7,1);
t279 = Ifges(7,4) + Ifges(6,5);
t287 = Ifges(6,6) - Ifges(7,6);
t150 = sin(qJ(5));
t153 = cos(qJ(5));
t206 = qJD(5) * t153;
t207 = qJD(5) * t150;
t141 = -pkin(3) * t148 - pkin(2);
t91 = pkin(4) * t127 - pkin(9) * t128 + t141;
t97 = -t136 * t151 + t137 * t154;
t281 = t288 * t150 + t278 * t153 + t91 * t206 - t207 * t97;
t232 = t150 * t91 + t153 * t97;
t280 = -qJD(5) * t232 - t278 * t150 + t288 * t153;
t120 = t127 * qJD(2);
t273 = qJD(5) + t120;
t121 = t128 * qJD(2);
t284 = t121 * Ifges(5,1) / 0.2e1;
t283 = qJ(6) * t123 + qJD(6) * t127 + t281;
t282 = -pkin(5) * t123 - t280;
t76 = qJD(3) * t128 + qJD(4) * t97;
t162 = t128 * t212;
t98 = qJD(1) * t162;
t277 = t76 - t98;
t208 = qJD(1) * t155;
t198 = t147 * t208;
t174 = qJD(3) - t198;
t114 = qJD(2) * t141 + t174;
t131 = qJD(2) * qJ(3) + t199;
t149 = cos(pkin(6));
t209 = qJD(1) * t149;
t139 = t148 * t209;
t222 = pkin(8) * qJD(2);
t94 = t139 + (-t131 - t222) * t146;
t105 = t148 * t131 + t146 * t209;
t95 = t148 * t222 + t105;
t52 = t151 * t94 + t154 * t95;
t48 = qJD(4) * pkin(9) + t52;
t55 = pkin(4) * t120 - pkin(9) * t121 + t114;
t15 = -t150 * t48 + t153 * t55;
t16 = t150 * t55 + t153 * t48;
t171 = t15 * t153 + t16 * t150;
t268 = qJD(6) - t15;
t10 = -pkin(5) * t273 + t268;
t11 = qJ(6) * t273 + t16;
t173 = t10 * t153 - t11 * t150;
t223 = Ifges(7,5) * t153;
t178 = Ifges(7,3) * t150 + t223;
t226 = Ifges(6,4) * t153;
t184 = -Ifges(6,2) * t150 + t226;
t189 = mrSges(7,1) * t150 - mrSges(7,3) * t153;
t191 = mrSges(6,1) * t150 + mrSges(6,2) * t153;
t103 = qJD(4) * t150 + t121 * t153;
t164 = t153 * qJD(4) - t121 * t150;
t51 = -t151 * t95 + t154 * t94;
t47 = -qJD(4) * pkin(4) - t51;
t23 = -pkin(5) * t164 - qJ(6) * t103 + t47;
t245 = t153 / 0.2e1;
t247 = t150 / 0.2e1;
t248 = -t150 / 0.2e1;
t252 = t103 / 0.2e1;
t254 = -t164 / 0.2e1;
t255 = t164 / 0.2e1;
t224 = Ifges(7,5) * t150;
t227 = Ifges(6,4) * t150;
t265 = t153 * t286 + t224 - t227;
t266 = -t150 * t287 + t153 * t279;
t101 = Ifges(6,4) * t164;
t225 = Ifges(7,5) * t164;
t271 = t103 * t286 + t273 * t279 + t101 - t225;
t100 = Ifges(7,5) * t103;
t34 = Ifges(7,6) * t273 - Ifges(7,3) * t164 + t100;
t228 = Ifges(6,4) * t103;
t37 = Ifges(6,2) * t164 + Ifges(6,6) * t273 + t228;
t157 = t173 * mrSges(7,2) - t171 * mrSges(6,3) + t178 * t254 + t184 * t255 + t189 * t23 + t191 * t47 + t247 * t34 + t248 * t37 + t265 * t252 + t266 * t273 / 0.2e1 + t271 * t245;
t276 = t114 * mrSges(5,2) - t51 * mrSges(5,3) - Ifges(5,4) * t120 + Ifges(5,5) * qJD(4) + t157 + t284;
t275 = t279 * t150 + t287 * t153;
t274 = t286 * t150 - t223 + t226;
t253 = -t103 / 0.2e1;
t250 = -t273 / 0.2e1;
t113 = qJD(2) * t123;
t112 = qJD(2) * t122;
t68 = qJD(5) * t164 - t112 * t153;
t69 = qJD(5) * t103 - t112 * t150;
t272 = (-Ifges(6,4) + Ifges(7,5)) * t69 + t286 * t68 + t279 * t113;
t211 = t146 ^ 2 + t148 ^ 2;
t270 = mrSges(4,3) * t211;
t175 = pkin(5) * t150 - qJ(6) * t153;
t269 = -qJD(6) * t150 + t175 * t273 - t52;
t126 = (qJD(3) + t198) * qJD(2);
t27 = qJD(4) * t51 - t127 * t126;
t213 = t147 * t152;
t197 = qJD(2) * t213;
t195 = qJD(1) * t197;
t67 = pkin(4) * t113 + pkin(9) * t112 + t195;
t3 = t150 * t67 + t153 * t27 + t55 * t206 - t207 * t48;
t4 = -qJD(5) * t16 - t150 * t27 + t153 * t67;
t264 = -t150 * t4 + t153 * t3;
t1 = qJ(6) * t113 + qJD(6) * t273 + t3;
t2 = -pkin(5) * t113 - t4;
t263 = t1 * t153 + t150 * t2;
t202 = Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1;
t203 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t204 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t220 = t120 * Ifges(5,2);
t260 = t203 * t164 + t204 * t103 - t202 * t273 + t10 * mrSges(7,1) + t16 * mrSges(6,2) + t52 * mrSges(5,3) + Ifges(6,6) * t254 + Ifges(7,6) * t255 + Ifges(5,6) * qJD(4) - t220 / 0.2e1 + Ifges(5,4) * t121 - t11 * mrSges(7,3) - t114 * mrSges(5,1) - t15 * mrSges(6,1) + t279 * t253 + (Ifges(6,3) + Ifges(7,2)) * t250;
t259 = m(4) / 0.2e1;
t258 = t68 / 0.2e1;
t257 = -t69 / 0.2e1;
t256 = t69 / 0.2e1;
t251 = t113 / 0.2e1;
t246 = -t153 / 0.2e1;
t118 = -t146 * t213 + t148 * t149;
t119 = t146 * t149 + t148 * t213;
t165 = t154 * t118 - t119 * t151;
t28 = qJD(4) * t52 + t126 * t128;
t240 = t28 * t165;
t239 = t28 * t267;
t40 = mrSges(6,1) * t113 - mrSges(6,3) * t68;
t41 = -t113 * mrSges(7,1) + t68 * mrSges(7,2);
t236 = t41 - t40;
t42 = -mrSges(6,2) * t113 - mrSges(6,3) * t69;
t43 = -mrSges(7,2) * t69 + mrSges(7,3) * t113;
t235 = t42 + t43;
t89 = pkin(4) * t121 + pkin(9) * t120;
t26 = t150 * t89 + t153 * t51;
t71 = mrSges(7,2) * t164 + mrSges(7,3) * t273;
t231 = mrSges(6,3) * t164;
t72 = -mrSges(6,2) * t273 + t231;
t234 = t71 + t72;
t230 = mrSges(6,3) * t103;
t73 = mrSges(6,1) * t273 - t230;
t74 = -mrSges(7,1) * t273 + mrSges(7,2) * t103;
t233 = t73 - t74;
t217 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t164 + mrSges(6,2) * t103 + t121 * mrSges(5,3);
t58 = -mrSges(7,1) * t164 - mrSges(7,3) * t103;
t205 = t58 + t217;
t201 = t150 * t212;
t200 = t147 ^ 2 * t208;
t77 = t113 * mrSges(5,1) - t112 * mrSges(5,2);
t196 = t211 * t126;
t194 = -t1 * t150 + t2 * t153;
t193 = -t3 * t150 - t4 * t153;
t192 = mrSges(6,1) * t153 - mrSges(6,2) * t150;
t190 = mrSges(7,1) * t153 + mrSges(7,3) * t150;
t183 = Ifges(6,2) * t153 + t227;
t177 = -Ifges(7,3) * t153 + t224;
t176 = pkin(5) * t153 + qJ(6) * t150;
t172 = t10 * t150 + t11 * t153;
t170 = t15 * t150 - t16 * t153;
t25 = -t150 * t51 + t153 * t89;
t44 = -t150 * t97 + t153 * t91;
t166 = -(-t131 * t146 + t139) * t146 + t105 * t148;
t81 = t118 * t151 + t119 * t154;
t60 = t150 * t81 + t153 * t212;
t163 = t166 * t155;
t160 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t156 = qJD(2) ^ 2;
t132 = -pkin(4) - t176;
t130 = -qJD(2) * pkin(2) + t174;
t110 = Ifges(7,2) * t113;
t109 = Ifges(6,3) * t113;
t106 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t120;
t88 = mrSges(5,1) * t120 + mrSges(5,2) * t121;
t66 = Ifges(7,4) * t68;
t65 = Ifges(6,5) * t68;
t64 = Ifges(6,6) * t69;
t63 = Ifges(7,6) * t69;
t61 = t153 * t81 - t201;
t57 = pkin(5) * t103 - qJ(6) * t164;
t53 = t128 * t175 - t267;
t50 = qJD(2) * t162 + qJD(4) * t81;
t49 = -qJD(2) * t161 + qJD(4) * t165;
t33 = -pkin(5) * t127 - t44;
t32 = qJ(6) * t127 + t232;
t31 = mrSges(6,1) * t69 + mrSges(6,2) * t68;
t30 = mrSges(7,1) * t69 - mrSges(7,3) * t68;
t20 = t68 * Ifges(6,4) - t69 * Ifges(6,2) + t113 * Ifges(6,6);
t19 = t68 * Ifges(7,5) + t113 * Ifges(7,6) + t69 * Ifges(7,3);
t18 = -pkin(5) * t121 - t25;
t17 = qJ(6) * t121 + t26;
t14 = -qJD(5) * t201 + t150 * t49 - t153 * t197 + t206 * t81;
t13 = -qJD(5) * t60 + t150 * t197 + t153 * t49;
t12 = -t175 * t122 + (qJD(5) * t176 - qJD(6) * t153) * t128 + t76;
t5 = pkin(5) * t69 - qJ(6) * t68 - qJD(6) * t103 + t28;
t6 = [-t81 * t113 * mrSges(5,3) + t49 * t106 + t235 * t61 + t236 * t60 - t233 * t14 + t234 * t13 - (-t112 * mrSges(5,3) + t30 + t31) * t165 + t205 * t50 + ((-mrSges(3,1) * t156 + (qJD(2) * (-mrSges(4,1) * t148 + mrSges(4,2) * t146) + t88) * qJD(2)) * t152 + (-t77 + (-mrSges(3,2) + t270) * t156) * t155) * t147 + m(7) * (t1 * t61 + t10 * t14 + t11 * t13 - t165 * t5 + t2 * t60 + t23 * t50) + m(6) * (t13 * t16 - t14 * t15 + t3 * t61 - t4 * t60 + t47 * t50 - t240) + m(5) * (t27 * t81 + t49 * t52 - t50 * t51 - t240) + m(4) * (-t118 * t146 + t119 * t148) * t126 + 0.2e1 * (t147 * t163 * t259 + (m(5) * (t114 * t147 - t200) / 0.2e1 + (t130 * t147 - t200) * t259) * t152) * qJD(2); (t112 * t267 - t113 * t97) * mrSges(5,3) - t267 * t31 + (t1 * t32 + t2 * t33 + t5 * t53 + (t12 - t98) * t23 + t283 * t11 + t282 * t10) * m(7) + t283 * t71 - (t284 + t276) * t122 + t280 * t73 + (t15 * t280 + t16 * t281 + t232 * t3 + t277 * t47 + t4 * t44 - t239) * m(6) + t281 * t72 + t282 * t74 + (-t114 * t199 + t141 * t195 + t27 * t97 - t277 * t51 + t278 * t52 - t239) * m(5) + t278 * t106 + t232 * t42 + t217 * t76 - t205 * t98 - t88 * t199 + (t19 * t247 + t5 * t189 + t178 * t256 + t184 * t257 - Ifges(5,1) * t112 - Ifges(5,4) * t113 + mrSges(5,2) * t195 + (mrSges(5,3) + t191) * t28 + t193 * mrSges(6,3) + t194 * mrSges(7,2) + (-mrSges(7,2) * t172 + mrSges(6,3) * t170 + t177 * t255 + t183 * t254 + t190 * t23 + t192 * t47 + t246 * t37 + t250 * t275 + t253 * t274) * qJD(5) + t265 * t258 + t266 * t251 + (qJD(5) * t271 + t20) * t248 + (qJD(5) * t34 + t272) * t245) * t128 + (qJD(2) * t174 * t211 + t196) * mrSges(4,3) + (t220 / 0.2e1 - t260) * t123 + (-pkin(2) * t195 + qJ(3) * t196 + qJD(3) * t166 - (t130 * t152 + t163) * t210) * m(4) + (mrSges(5,1) * t195 - t27 * mrSges(5,3) + t65 / 0.2e1 - t64 / 0.2e1 + t109 / 0.2e1 + t66 / 0.2e1 + t110 / 0.2e1 + t63 / 0.2e1 + Ifges(5,4) * t112 + t203 * t69 - t204 * t68 + (Ifges(5,2) + t202) * t113 + t160) * t127 + t33 * t41 + t32 * t43 + t44 * t40 + t53 * t30 + t12 * t58 + t141 * t77; t120 * t106 + (m(4) + m(5)) * t195 - t156 * t270 - t205 * t121 + (t234 * t273 - t236) * t153 + (-t233 * t273 + t235) * t150 - m(5) * (-t120 * t52 - t121 * t51) - m(4) * t166 * qJD(2) + t77 + (-t23 * t121 + t172 * t273 - t194) * m(7) + (-t121 * t47 - t170 * t273 - t193) * m(6); t272 * t247 + (-t10 * t18 - t11 * t17 + t132 * t5 + t269 * t23) * m(7) + t269 * t58 - t217 * t52 - t5 * t190 + (-mrSges(5,1) - t192) * t28 + t274 * t258 + t275 * t251 - ((Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t121 - t276) * t120 + t263 * mrSges(7,2) + t264 * mrSges(6,3) + (m(7) * t263 + m(6) * t264 + (-m(6) * t171 + m(7) * t173 - t150 * t234 - t153 * t233) * qJD(5) + t150 * t236 + t153 * t235) * pkin(9) + t260 * t121 + (-pkin(4) * t28 - t15 * t25 - t16 * t26 - t47 * t52) * m(6) + t177 * t256 + t183 * t257 + t20 * t245 + t19 * t246 + t157 * qJD(5) - t27 * mrSges(5,2) - pkin(4) * t31 - t17 * t71 - t26 * t72 - t25 * t73 - t18 * t74 - t51 * t106 - Ifges(5,5) * t112 - Ifges(5,6) * t113 + t132 * t30; (t230 + t233) * t16 + (t231 - t234) * t15 + t160 + (-t10 * t164 + t103 * t11) * mrSges(7,2) + t109 + t110 + t66 + t65 - t64 + t63 + (Ifges(7,3) * t103 + t225) * t255 + t37 * t252 - pkin(5) * t41 + qJ(6) * t43 - t57 * t58 + qJD(6) * t71 - t23 * (mrSges(7,1) * t103 - mrSges(7,3) * t164) - t47 * (mrSges(6,1) * t103 + mrSges(6,2) * t164) + (-t103 * t287 + t279 * t164) * t250 + (-pkin(5) * t2 + qJ(6) * t1 - t10 * t16 + t11 * t268 - t23 * t57) * m(7) + (-Ifges(6,2) * t103 + t101 + t271) * t254 + (t164 * t286 + t100 - t228 + t34) * t253; t103 * t58 - t273 * t71 + 0.2e1 * (t2 / 0.2e1 + t23 * t252 + t11 * t250) * m(7) + t41;];
tauc  = t6(:);
