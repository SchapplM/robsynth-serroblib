% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:42
% EndTime: 2019-12-31 20:05:57
% DurationCPUTime: 7.06s
% Computational Cost: add. (3632->472), mult. (9504->641), div. (0->0), fcn. (6310->6), ass. (0->211)
t287 = Ifges(5,1) + Ifges(6,1);
t286 = -Ifges(5,4) + Ifges(6,5);
t285 = Ifges(6,4) + Ifges(5,5);
t284 = -Ifges(5,6) + Ifges(6,6);
t291 = Ifges(5,3) + Ifges(6,2);
t177 = sin(pkin(8));
t178 = cos(pkin(8));
t179 = sin(qJ(4));
t181 = cos(qJ(4));
t144 = t177 * t181 + t178 * t179;
t182 = cos(qJ(2));
t189 = t144 * t182;
t113 = qJD(1) * t189;
t130 = t144 * qJD(4);
t227 = t113 - t130;
t232 = t177 * t179;
t143 = -t181 * t178 + t232;
t188 = t143 * t182;
t114 = qJD(1) * t188;
t129 = t143 * qJD(4);
t226 = -t114 + t129;
t180 = sin(qJ(2));
t225 = qJD(1) * t180;
t208 = t178 * t225;
t219 = t177 * qJD(2);
t141 = t208 + t219;
t207 = t177 * t225;
t191 = -t178 * qJD(2) + t207;
t184 = t181 * t141 - t179 * t191;
t186 = qJD(2) * t189;
t53 = qJD(1) * t186 + qJD(4) * t184;
t271 = t53 / 0.2e1;
t290 = Ifges(6,3) * t271;
t276 = -Ifges(3,6) / 0.2e1;
t289 = -Ifges(5,6) / 0.2e1;
t288 = Ifges(6,6) / 0.2e1;
t185 = qJD(2) * t188;
t88 = t179 * t141 + t181 * t191;
t52 = -qJD(1) * t185 - qJD(4) * t88;
t273 = t52 / 0.2e1;
t269 = -t88 / 0.2e1;
t268 = t88 / 0.2e1;
t224 = qJD(1) * t182;
t169 = qJD(4) - t224;
t256 = t169 / 0.2e1;
t266 = t184 / 0.2e1;
t218 = qJD(1) * qJD(2);
t204 = t180 * t218;
t283 = t285 * t204 + t286 * t53 + t287 * t52;
t248 = Ifges(6,5) * t88;
t86 = Ifges(5,4) * t88;
t282 = t285 * t169 + t287 * t184 + t248 - t86;
t174 = pkin(6) * t224;
t215 = pkin(3) * t224;
t136 = t177 * t215 + t174;
t281 = -t227 * pkin(4) + t226 * qJ(5) - qJD(5) * t144 - t136;
t150 = -pkin(2) * t182 - t180 * qJ(3) - pkin(1);
t135 = t150 * qJD(1);
t158 = qJD(2) * qJ(3) + t174;
t93 = t178 * t135 - t177 * t158;
t57 = -t141 * pkin(7) - t215 + t93;
t94 = t177 * t135 + t178 * t158;
t61 = -pkin(7) * t191 + t94;
t19 = -t179 * t61 + t181 * t57;
t280 = qJD(5) - t19;
t279 = t291 * t204 + t284 * t53 + t285 * t52;
t20 = t179 * t57 + t181 * t61;
t228 = t178 * t182;
t192 = pkin(3) * t180 - pkin(7) * t228;
t187 = t192 * qJD(2);
t197 = pkin(2) * t180 - qJ(3) * t182;
t126 = qJD(2) * t197 - t180 * qJD(3);
t115 = t126 * qJD(1);
t173 = pkin(6) * t225;
t148 = (qJD(3) - t173) * qJD(2);
t74 = t178 * t115 - t177 * t148;
t55 = qJD(1) * t187 + t74;
t203 = t182 * t218;
t201 = t177 * t203;
t75 = t177 * t115 + t178 * t148;
t60 = -pkin(7) * t201 + t75;
t4 = -qJD(4) * t20 - t179 * t60 + t181 * t55;
t140 = t178 * t150;
t229 = t178 * t180;
t92 = -pkin(7) * t229 + t140 + (-pkin(6) * t177 - pkin(3)) * t182;
t110 = pkin(6) * t228 + t177 * t150;
t231 = t177 * t180;
t99 = -pkin(7) * t231 + t110;
t244 = t179 * t92 + t181 * t99;
t223 = qJD(2) * t180;
t214 = pkin(6) * t223;
t97 = t178 * t126 + t177 * t214;
t69 = t187 + t97;
t117 = t177 * t126;
t230 = t177 * t182;
t190 = -pkin(6) * t229 - pkin(7) * t230;
t77 = qJD(2) * t190 + t117;
t9 = -qJD(4) * t244 - t179 * t77 + t181 * t69;
t220 = qJD(4) * t181;
t222 = qJD(4) * t179;
t3 = t179 * t55 + t181 * t60 + t57 * t220 - t222 * t61;
t1 = qJ(5) * t204 + qJD(5) * t169 + t3;
t2 = -pkin(4) * t204 - t4;
t278 = -t4 * mrSges(5,1) + t2 * mrSges(6,1) + t3 * mrSges(5,2) - t1 * mrSges(6,3);
t17 = -pkin(4) * t169 + t280;
t18 = qJ(5) * t169 + t20;
t211 = Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t212 = t288 + t289;
t213 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t277 = t211 * t169 + t213 * t184 + t212 * t88 + t18 * mrSges(6,3) + t19 * mrSges(5,1) + t93 * mrSges(4,1) + qJD(2) * t276 - (Ifges(3,4) * t180 + Ifges(3,2) * t182) * qJD(1) / 0.2e1 + Ifges(5,6) * t269 + Ifges(6,6) * t268 - Ifges(4,6) * t191 / 0.2e1 - Ifges(4,3) * t224 / 0.2e1 + Ifges(4,5) * t141 - pkin(6) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t224) - t17 * mrSges(6,1) - t20 * mrSges(5,2) - t94 * mrSges(4,2) + t285 * t266 + t291 * t256;
t275 = Ifges(6,5) * t273 + t204 * t288 + t290;
t274 = -Ifges(5,4) * t52 / 0.2e1 + Ifges(5,2) * t271 + t204 * t289;
t272 = -t53 / 0.2e1;
t267 = -t184 / 0.2e1;
t265 = pkin(1) * mrSges(3,1);
t264 = pkin(1) * mrSges(3,2);
t257 = -t169 / 0.2e1;
t255 = -t177 / 0.2e1;
t254 = t177 / 0.2e1;
t253 = t178 / 0.2e1;
t251 = mrSges(5,3) * t88;
t250 = mrSges(5,3) * t184;
t249 = Ifges(5,4) * t184;
t247 = pkin(7) + qJ(3);
t65 = -mrSges(5,2) * t169 - t251;
t68 = -mrSges(6,2) * t88 + mrSges(6,3) * t169;
t246 = t65 + t68;
t66 = mrSges(5,1) * t169 - t250;
t67 = -mrSges(6,1) * t169 + mrSges(6,2) * t184;
t245 = t66 - t67;
t146 = t197 * qJD(1);
t105 = pkin(6) * t207 + t178 * t146;
t76 = qJD(1) * t192 + t105;
t131 = t177 * t146;
t91 = qJD(1) * t190 + t131;
t31 = t179 * t76 + t181 * t91;
t243 = mrSges(4,2) * t178;
t242 = Ifges(4,1) * t141;
t241 = Ifges(4,1) * t178;
t240 = Ifges(4,4) * t141;
t239 = Ifges(4,4) * t177;
t238 = Ifges(4,4) * t178;
t236 = Ifges(4,5) * t178;
t235 = Ifges(4,2) * t177;
t234 = Ifges(4,6) * t177;
t233 = Ifges(3,5) * qJD(2);
t116 = mrSges(4,1) * t201 + t203 * t243;
t168 = pkin(6) * t203;
t125 = pkin(3) * t201 + t168;
t137 = (pkin(3) * t219 + pkin(6) * qJD(2)) * t182;
t147 = pkin(3) * t231 + t180 * pkin(6);
t221 = qJD(4) * t180;
t210 = Ifges(4,5) * t224;
t209 = Ifges(4,6) * t224;
t171 = -pkin(3) * t178 - pkin(2);
t15 = t53 * mrSges(5,1) + t52 * mrSges(5,2);
t14 = t53 * mrSges(6,1) - t52 * mrSges(6,3);
t149 = -qJD(2) * pkin(2) + qJD(3) + t173;
t202 = m(4) * t149 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t191 + t141 * mrSges(4,2) + mrSges(3,3) * t225;
t200 = mrSges(4,1) * t177 + t243;
t199 = -t239 + t241;
t198 = -t235 + t238;
t30 = -t179 * t91 + t181 * t76;
t41 = -t179 * t99 + t181 * t92;
t152 = t247 * t177;
t153 = t247 * t178;
t193 = -t181 * t152 - t153 * t179;
t103 = -t152 * t179 + t153 * t181;
t39 = -mrSges(6,1) * t204 + t52 * mrSges(6,2);
t8 = t179 * t69 + t181 * t77 + t92 * t220 - t222 * t99;
t104 = pkin(3) * t191 + t149;
t172 = Ifges(3,4) * t224;
t134 = Ifges(3,1) * t225 + t172 + t233;
t124 = (mrSges(4,1) * t180 - mrSges(4,3) * t228) * t218;
t123 = (-mrSges(4,2) * t180 - mrSges(4,3) * t230) * t218;
t122 = t143 * t180;
t121 = t144 * t180;
t112 = -mrSges(4,1) * t224 - t141 * mrSges(4,3);
t111 = mrSges(4,2) * t224 - mrSges(4,3) * t191;
t109 = -pkin(6) * t230 + t140;
t106 = -pkin(6) * t208 + t131;
t101 = (Ifges(4,5) * t180 + t182 * t199) * t218;
t100 = (Ifges(4,6) * t180 + t182 * t198) * t218;
t98 = -t178 * t214 + t117;
t85 = Ifges(6,5) * t184;
t84 = pkin(4) * t143 - qJ(5) * t144 + t171;
t82 = -Ifges(4,4) * t191 - t210 + t242;
t81 = -Ifges(4,2) * t191 - t209 + t240;
t72 = t220 * t229 - t221 * t232 + t186;
t71 = -t144 * t221 - t185;
t64 = qJD(3) * t144 + qJD(4) * t103;
t63 = -qJD(3) * t143 + qJD(4) * t193;
t58 = pkin(4) * t121 + qJ(5) * t122 + t147;
t40 = -mrSges(5,2) * t204 - mrSges(5,3) * t53;
t38 = mrSges(5,1) * t204 - mrSges(5,3) * t52;
t37 = -mrSges(6,2) * t53 + mrSges(6,3) * t204;
t36 = mrSges(5,1) * t88 + mrSges(5,2) * t184;
t35 = mrSges(6,1) * t88 - mrSges(6,3) * t184;
t34 = pkin(4) * t184 + qJ(5) * t88;
t33 = pkin(4) * t182 - t41;
t32 = -qJ(5) * t182 + t244;
t27 = -Ifges(5,2) * t88 + Ifges(5,6) * t169 + t249;
t24 = Ifges(6,6) * t169 + Ifges(6,3) * t88 + t85;
t23 = -pkin(4) * t225 - t30;
t22 = qJ(5) * t225 + t31;
t21 = t88 * pkin(4) - qJ(5) * t184 + t104;
t16 = pkin(4) * t72 - qJ(5) * t71 + qJD(5) * t122 + t137;
t7 = -pkin(4) * t223 - t9;
t6 = qJ(5) * t223 - qJD(5) * t182 + t8;
t5 = pkin(4) * t53 - qJ(5) * t52 - qJD(5) * t184 + t125;
t10 = [-t279 * t182 / 0.2e1 + (t101 * t253 + t100 * t255 + pkin(6) * t116 + (-t177 * t75 - t178 * t74) * mrSges(4,3) + ((Ifges(4,6) * t253 + t276) * qJD(2) + t277) * qJD(2) + (-0.2e1 * t265 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t236 / 0.2e1 - t234) * t180 - t213 * t122 + t212 * t121) * t218) * t180 + t110 * t123 + t109 * t124 + t98 * t111 + t97 * t112 + t8 * t65 + t9 * t66 + t7 * t67 + t6 * t68 + t58 * t14 + t16 * t35 + t32 * t37 + t33 * t39 + t41 * t38 + (t125 * mrSges(5,1) + t5 * mrSges(6,1) - mrSges(6,2) * t1 - mrSges(5,3) * t3 - Ifges(5,2) * t272 + t273 * t286 + t274 + t275 + t290) * t121 + t244 * t40 + m(5) * (t104 * t137 + t125 * t147 + t19 * t9 + t20 * t8 + t244 * t3 + t4 * t41) + (t104 * mrSges(5,1) + t21 * mrSges(6,1) + t24 / 0.2e1 - t27 / 0.2e1 - t18 * mrSges(6,2) - t20 * mrSges(5,3) + Ifges(6,3) * t268 - Ifges(5,2) * t269 + t286 * t266 + t284 * t256) * t72 + m(6) * (t1 * t32 + t16 * t21 + t17 * t7 + t18 * t6 + t2 * t33 + t5 * t58) + m(4) * (t74 * t109 + t75 * t110 + t93 * t97 + t94 * t98) - (t125 * mrSges(5,2) + mrSges(6,2) * t2 - mrSges(5,3) * t4 - t5 * mrSges(6,3) + Ifges(5,4) * t272 + Ifges(6,5) * t271 + t287 * t273) * t122 - t283 * t122 / 0.2e1 + ((t134 / 0.2e1 + t149 * t200 + t141 * t199 / 0.2e1 + (t198 * t253 + Ifges(3,5) / 0.2e1) * qJD(2) + t81 * t255 + t82 * t253 + (-t177 * t94 - t178 * t93) * mrSges(4,3) + t202 * pkin(6)) * qJD(2) + t75 * mrSges(4,2) - t74 * mrSges(4,1) - Ifges(6,6) * t271 - Ifges(5,6) * t272 - t285 * t273 + (-0.2e1 * t264 + (0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * t236 + 0.3e1 / 0.2e1 * t234) * t182 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(4,3) + Ifges(4,1) * t178 ^ 2 / 0.2e1 + (-0.3e1 / 0.2e1 * t238 + t235) * t177 + (m(4) * pkin(6) + t200) * pkin(6) - t211) * t180) * t218 + t278) * t182 + (t287 * t266 + mrSges(5,2) * t104 + t17 * mrSges(6,2) - t19 * mrSges(5,3) - mrSges(6,3) * t21 + Ifges(5,4) * t269 + Ifges(6,5) * t268 + t282 / 0.2e1 + t285 * t256) * t71 + t137 * t36 + t147 * t15; ((-t134 / 0.2e1 - t172 / 0.2e1 + t233 / 0.2e1 + qJD(1) * t264 + (-t149 * mrSges(4,2) + t93 * mrSges(4,3) - t242 / 0.2e1 - t82 / 0.2e1 + t210 / 0.2e1) * t178 + ((-m(4) * pkin(2) - mrSges(4,1) * t178 - mrSges(3,1)) * qJD(2) - t202) * pkin(6) + (-t149 * mrSges(4,1) + t94 * mrSges(4,3) + t240 / 0.2e1 + t81 / 0.2e1 - t209 / 0.2e1 + (pkin(6) * mrSges(4,2) + t241 / 0.2e1 - t239 / 0.2e1) * qJD(2)) * t177) * t182 + ((-Ifges(3,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + t198 * t254) * t224 + (t265 + (Ifges(3,4) / 0.2e1 + t234 / 0.2e1) * t180) * qJD(1) + (pkin(6) * mrSges(3,2) + Ifges(4,5) * t254 + t212 * t143 + t213 * t144 + t276) * qJD(2) - t277) * t180) * qJD(1) - pkin(2) * t116 - t106 * t111 - t105 * t112 + t84 * t14 - t31 * t65 - t30 * t66 - t23 * t67 - t22 * t68 - m(4) * (t93 * t105 + t94 * t106) + (-t143 * t3 - t144 * t4 + t19 * t226 + t20 * t227) * mrSges(5,3) + (-t1 * t143 + t144 * t2 - t17 * t226 + t18 * t227) * mrSges(6,2) + (-mrSges(6,1) * t227 + mrSges(6,3) * t226) * t21 + (-mrSges(5,1) * t227 - mrSges(5,2) * t226) * t104 + (t27 - t24) * (t113 / 0.2e1 - t130 / 0.2e1) + m(4) * ((-t177 * t93 + t178 * t94) * qJD(3) + (-t74 * t177 + t75 * t178) * qJ(3)) + (qJ(3) * t123 + qJD(3) * t111 + t75 * mrSges(4,3) + t100 / 0.2e1) * t178 + (-Ifges(5,4) * t129 - Ifges(6,5) * t114 - Ifges(5,2) * t130 + Ifges(6,3) * t113) * t269 + (-Ifges(5,4) * t114 - Ifges(6,5) * t129 - Ifges(5,2) * t113 + Ifges(6,3) * t130) * t268 + (t113 * t286 - t114 * t287) * t267 + t282 * (t114 / 0.2e1 - t129 / 0.2e1) + (t113 * t284 - t114 * t285) * t257 + (t193 * t4 + t103 * t3 - t104 * t136 + t125 * t171 + (-t31 + t63) * t20 + (-t30 - t64) * t19) * m(5) + (t1 * t103 - t193 * t2 + t5 * t84 + t281 * t21 + (-t22 + t63) * t18 + (t64 - t23) * t17) * m(6) - (t39 - t38) * t193 + (-qJ(3) * t124 - qJD(3) * t112 - t74 * mrSges(4,3) + t101 / 0.2e1) * t177 + (t143 * t286 + t144 * t287) * t273 + (-t129 * t287 + t130 * t286) * t266 + t281 * t35 + t283 * t144 / 0.2e1 + (-t129 * t285 + t130 * t284) * t256 + (t37 + t40) * t103 - t136 * t36 + t5 * (mrSges(6,1) * t143 - mrSges(6,3) * t144) + t125 * (mrSges(5,1) * t143 + mrSges(5,2) * t144) + t171 * t15 + (Ifges(6,5) * t144 + Ifges(6,3) * t143) * t271 + (Ifges(5,4) * t144 - Ifges(5,2) * t143) * t272 + t143 * t274 + t143 * t275 - t245 * t64 + t246 * t63; t191 * t111 + t141 * t112 + t245 * t184 + t246 * t88 + t116 + t14 + t15 + (-t17 * t184 + t18 * t88 + t5) * m(6) + (t184 * t19 + t20 * t88 + t125) * m(5) + (t93 * t141 + t191 * t94 + t168) * m(4); qJD(5) * t68 - t34 * t35 + qJ(5) * t37 - pkin(4) * t39 - t104 * (mrSges(5,1) * t184 - mrSges(5,2) * t88) - t21 * (mrSges(6,1) * t184 + mrSges(6,3) * t88) + (t17 * t88 + t18 * t184) * mrSges(6,2) + (Ifges(6,3) * t184 - t248) * t269 + (t184 * t284 - t285 * t88) * t257 + (-Ifges(5,2) * t184 + t282 - t86) * t268 + (t245 + t250) * t20 + (-t246 - t251) * t19 - t278 + (-t287 * t88 + t24 - t249 + t85) * t267 + (-pkin(4) * t2 + qJ(5) * t1 - t17 * t20 + t18 * t280 - t21 * t34) * m(6) + t27 * t266 + t279; -t169 * t68 + t184 * t35 + 0.2e1 * (t2 / 0.2e1 + t18 * t257 + t21 * t266) * m(6) + t39;];
tauc = t10(:);
