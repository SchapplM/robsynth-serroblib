% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:50
% EndTime: 2020-01-03 12:08:57
% DurationCPUTime: 3.16s
% Computational Cost: add. (4938->342), mult. (8596->485), div. (0->0), fcn. (5745->8), ass. (0->183)
t208 = cos(qJ(3));
t196 = t208 * qJD(4);
t205 = sin(qJ(3));
t257 = -qJ(4) - pkin(7);
t229 = qJD(3) * t257;
t160 = t205 * t229 + t196;
t161 = -t205 * qJD(4) + t208 * t229;
t202 = sin(pkin(9));
t203 = cos(pkin(9));
t171 = -t202 * t205 + t203 * t208;
t209 = cos(qJ(2));
t252 = qJD(1) * pkin(1);
t235 = t209 * t252;
t274 = t203 * t160 + t202 * t161 - t171 * t235;
t172 = t202 * t208 + t203 * t205;
t271 = -t160 * t202 + t203 * t161 + t172 * t235;
t166 = t171 * qJD(3);
t259 = pkin(8) * t166;
t289 = -t259 + t271;
t165 = t172 * qJD(3);
t157 = t165 * pkin(8);
t288 = t157 - t274;
t198 = qJD(3) + qJD(5);
t199 = qJD(1) + qJD(2);
t145 = t171 * t199;
t146 = t172 * t199;
t204 = sin(qJ(5));
t207 = cos(qJ(5));
t226 = t207 * t145 - t146 * t204;
t83 = t145 * t204 + t146 * t207;
t264 = Ifges(6,4) * t83;
t260 = pkin(8) * t146;
t206 = sin(qJ(2));
t236 = t206 * t252;
t179 = pkin(7) * t199 + t236;
t227 = qJ(4) * t199 + t179;
t132 = t227 * t208;
t111 = t202 * t132;
t131 = t227 * t205;
t117 = qJD(3) * pkin(3) - t131;
t62 = t203 * t117 - t111;
t45 = qJD(3) * pkin(4) - t260 + t62;
t261 = pkin(8) * t145;
t243 = t203 * t132;
t63 = t202 * t117 + t243;
t46 = t63 + t261;
t20 = -t204 * t46 + t207 * t45;
t136 = t199 * t166;
t233 = qJD(2) * t252;
t224 = t209 * t233;
t184 = t208 * t224;
t215 = qJD(3) * t227;
t88 = t196 * t199 - t205 * t215 + t184;
t89 = (-qJD(4) * t199 - t224) * t205 - t208 * t215;
t39 = -t202 * t88 + t203 * t89;
t25 = -pkin(8) * t136 + t39;
t135 = t199 * t165;
t40 = t202 * t89 + t203 * t88;
t26 = -pkin(8) * t135 + t40;
t3 = qJD(5) * t20 + t204 * t25 + t207 * t26;
t76 = Ifges(6,4) * t226;
t32 = Ifges(6,1) * t83 + Ifges(6,5) * t198 + t76;
t36 = qJD(5) * t226 - t135 * t204 + t136 * t207;
t37 = -qJD(5) * t83 - t135 * t207 - t136 * t204;
t21 = t204 * t45 + t207 * t46;
t4 = -qJD(5) * t21 - t204 * t26 + t207 * t25;
t193 = -t208 * pkin(3) - pkin(2);
t142 = t193 * t199 + qJD(4) - t235;
t90 = -t145 * pkin(4) + t142;
t287 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t36 + Ifges(6,6) * t37 - (Ifges(6,5) * t226 - Ifges(6,6) * t83) * t198 / 0.2e1 - (-Ifges(6,2) * t83 + t32 + t76) * t226 / 0.2e1 - t90 * (mrSges(6,1) * t83 + mrSges(6,2) * t226) - (Ifges(6,1) * t226 - t264) * t83 / 0.2e1;
t286 = -t165 / 0.2e1;
t285 = t166 / 0.2e1;
t284 = t226 / 0.2e1;
t240 = qJD(3) * t205;
t118 = -t179 * t240 + t184;
t239 = qJD(3) * t208;
t119 = -t179 * t239 - t205 * t224;
t218 = t118 * t208 - t119 * t205;
t283 = t20 * t226 + t21 * t83;
t31 = Ifges(6,2) * t226 + Ifges(6,6) * t198 + t264;
t281 = t31 / 0.2e1;
t185 = t257 * t205;
t197 = t208 * qJ(4);
t186 = pkin(7) * t208 + t197;
t115 = t203 * t185 - t186 * t202;
t258 = pkin(8) * t172;
t93 = t115 - t258;
t116 = t202 * t185 + t203 * t186;
t168 = t171 * pkin(8);
t94 = t168 + t116;
t41 = -t204 * t94 + t207 * t93;
t280 = qJD(5) * t41 + t204 * t289 - t288 * t207;
t42 = t204 * t93 + t207 * t94;
t279 = -qJD(5) * t42 + t288 * t204 + t207 * t289;
t190 = pkin(3) * t203 + pkin(4);
t262 = pkin(3) * t202;
t158 = t190 * t207 - t204 * t262;
t69 = t131 * t202 - t243;
t47 = t69 - t261;
t70 = -t203 * t131 - t111;
t48 = t70 - t260;
t273 = t158 * qJD(5) - t204 * t47 - t207 * t48;
t159 = t190 * t204 + t207 * t262;
t272 = -t159 * qJD(5) + t204 * t48 - t207 * t47;
t268 = t83 / 0.2e1;
t266 = t146 / 0.2e1;
t263 = pkin(1) * t209;
t256 = Ifges(4,4) * t205;
t254 = Ifges(5,4) * t146;
t253 = pkin(1) * qJD(2);
t251 = Ifges(4,5) * qJD(3);
t250 = Ifges(4,6) * qJD(3);
t246 = t199 * t205;
t245 = t199 * t208;
t191 = pkin(1) * t206 + pkin(7);
t241 = -qJ(4) - t191;
t225 = qJD(3) * t241;
t234 = t209 * t253;
t106 = t205 * t225 + t208 * t234 + t196;
t107 = (-qJD(4) - t234) * t205 + t208 * t225;
t53 = t203 * t106 + t202 * t107;
t169 = t241 * t205;
t170 = t191 * t208 + t197;
t98 = t202 * t169 + t203 * t170;
t189 = t206 * t233;
t232 = t199 * t240;
t162 = pkin(3) * t232 + t189;
t238 = -qJD(1) - t199;
t237 = -qJD(2) + t199;
t194 = pkin(3) * t240;
t12 = -t37 * mrSges(6,1) + t36 * mrSges(6,2);
t127 = pkin(4) * t165 + t194;
t75 = t135 * mrSges(5,1) + t136 * mrSges(5,2);
t228 = t179 * (t205 ^ 2 + t208 ^ 2);
t52 = -t106 * t202 + t203 * t107;
t97 = t203 * t169 - t170 * t202;
t100 = t171 * t204 + t172 * t207;
t99 = t171 * t207 - t172 * t204;
t50 = qJD(5) * t99 - t165 * t204 + t166 * t207;
t222 = -t4 * t100 - t20 * t50;
t220 = -mrSges(4,1) * t208 + mrSges(4,2) * t205;
t219 = -t166 * t62 - t172 * t39;
t73 = t97 - t258;
t74 = t168 + t98;
t29 = -t204 * t74 + t207 * t73;
t30 = t204 * t73 + t207 * t74;
t177 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t246;
t178 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t245;
t217 = t208 * t177 + t205 * t178;
t216 = t205 * t177 - t208 * t178;
t138 = -t171 * pkin(4) + t193;
t214 = (Ifges(4,2) * t208 + t256) * t199;
t213 = (mrSges(4,1) * t205 + mrSges(4,2) * t208) * qJD(3);
t149 = t214 + t250;
t188 = Ifges(4,4) * t245;
t150 = Ifges(4,1) * t246 + t188 + t251;
t180 = -t199 * pkin(2) - t235;
t51 = -qJD(5) * t100 - t165 * t207 - t166 * t204;
t78 = Ifges(5,2) * t145 + Ifges(5,6) * qJD(3) + t254;
t137 = Ifges(5,4) * t145;
t79 = Ifges(5,1) * t146 + Ifges(5,5) * qJD(3) + t137;
t95 = pkin(4) * t135 + t162;
t210 = (t51 * t284 + t99 * t37) * Ifges(6,2) + (t37 * t100 + t51 * t268 + t50 * t284 + t99 * t36) * Ifges(6,4) + (-t135 * t172 + t171 * t136 + t145 * t285 - t165 * t266) * Ifges(5,4) + (-t171 * t135 + t145 * t286) * Ifges(5,2) + t218 * mrSges(4,3) - (t214 + t149) * t240 / 0.2e1 + (t36 * t100 + t50 * t268) * Ifges(6,1) + (t136 * t172 + t166 * t266) * Ifges(5,1) + (t21 * t51 + t3 * t99) * mrSges(6,3) + (-t63 * t165 + t40 * t171) * mrSges(5,3) + ((0.3e1 * Ifges(4,4) * t208 + (Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t205) * t199 + t150) * t239 / 0.2e1 + (Ifges(4,1) * t208 - t256) * t232 + t50 * t32 / 0.2e1 + qJD(3) ^ 2 * (Ifges(4,5) * t208 - Ifges(4,6) * t205) / 0.2e1 + t79 * t285 + t78 * t286 + t220 * t189 + t51 * t281 + t95 * (-mrSges(6,1) * t99 + mrSges(6,2) * t100) + t90 * (-mrSges(6,1) * t51 + mrSges(6,2) * t50) + t162 * (-mrSges(5,1) * t171 + mrSges(5,2) * t172) + t198 * (Ifges(6,5) * t50 + Ifges(6,6) * t51) / 0.2e1 + t180 * t213 + qJD(3) * (Ifges(5,5) * t166 - Ifges(5,6) * t165) / 0.2e1 + t142 * (mrSges(5,1) * t165 + mrSges(5,2) * t166);
t195 = t206 * t253;
t192 = -pkin(2) - t263;
t181 = t193 - t263;
t176 = t195 + t194;
t163 = t220 * t199;
t151 = t199 * t213;
t124 = t138 - t263;
t121 = qJD(3) * mrSges(5,1) - t146 * mrSges(5,3);
t120 = -qJD(3) * mrSges(5,2) + t145 * mrSges(5,3);
t110 = t127 + t195;
t105 = pkin(3) * t246 + pkin(4) * t146;
t87 = -mrSges(5,1) * t145 + mrSges(5,2) * t146;
t68 = mrSges(6,1) * t198 - mrSges(6,3) * t83;
t67 = -mrSges(6,2) * t198 + mrSges(6,3) * t226;
t44 = -t157 + t53;
t43 = t52 - t259;
t38 = -mrSges(6,1) * t226 + mrSges(6,2) * t83;
t6 = -qJD(5) * t30 - t204 * t44 + t207 * t43;
t5 = qJD(5) * t29 + t204 * t43 + t207 * t44;
t1 = [(-t29 * t36 + t30 * t37 + t222) * mrSges(6,3) + (-t135 * t98 - t136 * t97 + t219) * mrSges(5,3) + m(6) * (t110 * t90 + t124 * t95 + t20 * t6 + t21 * t5 + t29 * t4 + t3 * t30) + m(5) * (t142 * t176 + t162 * t181 + t39 * t97 + t40 * t98 + t52 * t62 + t53 * t63) + (m(4) * t218 - qJD(3) * t217) * t191 + ((t163 + m(4) * (qJD(1) * t192 + t180) + t238 * mrSges(3,1)) * t206 + (m(4) * t228 + mrSges(3,2) * t238 - t216) * t209) * t253 + t210 + t53 * t120 + t52 * t121 + t124 * t12 + t110 * t38 + t5 * t67 + t6 * t68 + t176 * t87 + t181 * t75 + t192 * t151; t274 * t120 + (t205 * pkin(3) * t87 - pkin(7) * t217) * qJD(3) + t271 * t121 + t279 * t68 + t280 * t67 + (-t41 * t36 + t42 * t37 + t222) * mrSges(6,3) + (-t115 * t136 - t116 * t135 + t219) * mrSges(5,3) + t210 - pkin(2) * t151 + t138 * t12 + t127 * t38 + t193 * t75 + ((mrSges(3,2) * t237 + t216) * t209 + (mrSges(3,1) * t237 - t163 - t38 - t87) * t206) * t252 + (t138 * t95 + t3 * t42 + t4 * t41 + (t127 - t236) * t90 + t280 * t21 + t279 * t20) * m(6) + (t115 * t39 + t116 * t40 + t162 * t193 + t274 * t63 + t271 * t62 + (t194 - t236) * t142) * m(5) + (-pkin(2) * t189 + t218 * pkin(7) - (t180 * t206 + t209 * t228) * t252) * m(4); (-t158 * t36 + t159 * t37 + t283) * mrSges(6,3) + t83 * t281 + (t62 * t145 + t63 * t146 + (-t135 * t202 - t136 * t203) * pkin(3)) * mrSges(5,3) - t142 * (t146 * mrSges(5,1) + t145 * mrSges(5,2)) - qJD(3) * (Ifges(5,5) * t145 - Ifges(5,6) * t146) / 0.2e1 - (-Ifges(5,2) * t146 + t137 + t79) * t145 / 0.2e1 - t146 * (Ifges(5,1) * t145 - t254) / 0.2e1 + t287 + ((t251 / 0.2e1 - t180 * mrSges(4,2) - t150 / 0.2e1 - t188 / 0.2e1) * t208 + (-t250 / 0.2e1 - t180 * mrSges(4,1) + t149 / 0.2e1 + (t256 / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t208) * t199 + (-m(5) * t142 - t87) * pkin(3)) * t205) * t199 + m(5) * (t202 * t40 + t203 * t39) * pkin(3) + t39 * mrSges(5,1) - t40 * mrSges(5,2) + t272 * t68 + (-t90 * t105 + t4 * t158 + t3 * t159 + t20 * t272 + t21 * t273) * m(6) + t273 * t67 - m(5) * (t62 * t69 + t63 * t70) + t217 * t179 - Ifges(5,6) * t135 + Ifges(5,5) * t136 - t69 * t121 - t105 * t38 - t118 * mrSges(4,2) + t119 * mrSges(4,1) - t70 * t120 + t78 * t266; -t145 * t120 + t146 * t121 - t226 * t67 + t83 * t68 + t12 + t75 + (t20 * t83 - t21 * t226 + t95) * m(6) + (-t63 * t145 + t62 * t146 + t162) * m(5); t283 * mrSges(6,3) - t20 * t67 + t21 * t68 + t31 * t268 + t287;];
tauc = t1(:);
