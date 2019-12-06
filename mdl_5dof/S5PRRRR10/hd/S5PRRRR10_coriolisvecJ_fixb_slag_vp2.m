% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:16
% EndTime: 2019-12-05 17:23:40
% DurationCPUTime: 8.61s
% Computational Cost: add. (5610->513), mult. (15700->764), div. (0->0), fcn. (12307->12), ass. (0->253)
t167 = cos(pkin(6));
t175 = cos(qJ(3));
t176 = cos(qJ(2));
t239 = t175 * t176;
t171 = sin(qJ(3));
t172 = sin(qJ(2));
t243 = t171 * t172;
t190 = -t167 * t243 + t239;
t166 = sin(pkin(5));
t238 = qJD(1) * t166;
t113 = t190 * t238;
t165 = sin(pkin(6));
t248 = t165 * t171;
t161 = pkin(8) * t248;
t244 = t167 * t175;
t145 = pkin(2) * t244 - t161;
t139 = t145 * qJD(3);
t308 = t113 - t139;
t209 = pkin(3) * t171 - pkin(9) * t175;
t189 = t209 * qJD(3);
t225 = t172 * t238;
t218 = t165 * t225;
t309 = t165 * t189 - t218;
t245 = t167 * t171;
t247 = t165 * t175;
t146 = pkin(2) * t245 + pkin(8) * t247;
t132 = pkin(9) * t167 + t146;
t210 = -pkin(3) * t175 - pkin(9) * t171;
t133 = (-pkin(2) + t210) * t165;
t170 = sin(qJ(4));
t174 = cos(qJ(4));
t231 = qJD(4) * t174;
t232 = qJD(4) * t170;
t298 = -t132 * t232 + t133 * t231 + t309 * t170 - t308 * t174;
t241 = t172 * t175;
t242 = t171 * t176;
t192 = t167 * t241 + t242;
t112 = t192 * t238;
t140 = t146 * qJD(3);
t295 = t140 - t112;
t234 = qJD(3) * t165;
t220 = t171 * t234;
t307 = pkin(10) * t220 + t298;
t143 = -t174 * t167 + t170 * t248;
t233 = qJD(3) * t175;
t221 = t174 * t233;
t106 = -qJD(4) * t143 + t165 * t221;
t144 = t167 * t170 + t174 * t248;
t222 = t170 * t233;
t107 = qJD(4) * t144 + t165 * t222;
t306 = pkin(4) * t107 - pkin(10) * t106 + t295;
t169 = sin(qJ(5));
t173 = cos(qJ(5));
t235 = qJD(2) * t175;
t223 = t165 * t235;
t156 = qJD(4) - t223;
t160 = qJD(2) * t167 + qJD(3);
t152 = qJD(2) * pkin(2) + t176 * t238;
t168 = cos(pkin(5));
t237 = qJD(1) * t168;
t187 = t152 * t167 + t165 * t237;
t236 = qJD(2) * t165;
t147 = pkin(8) * t236 + t225;
t249 = t147 * t175;
t85 = t171 * t187 + t249;
t71 = pkin(9) * t160 + t85;
t159 = t167 * t237;
t91 = t159 + (qJD(2) * t210 - t152) * t165;
t34 = t170 * t91 + t174 * t71;
t29 = pkin(10) * t156 + t34;
t224 = t171 * t236;
t122 = t160 * t174 - t170 * t224;
t123 = t160 * t170 + t174 * t224;
t141 = t171 * t147;
t84 = t175 * t187 - t141;
t70 = -pkin(3) * t160 - t84;
t35 = -pkin(4) * t122 - pkin(10) * t123 + t70;
t12 = -t169 * t29 + t173 * t35;
t13 = t169 * t35 + t173 * t29;
t198 = t12 * t173 + t13 * t169;
t200 = Ifges(6,5) * t173 - Ifges(6,6) * t169;
t262 = Ifges(6,4) * t173;
t202 = -Ifges(6,2) * t169 + t262;
t263 = Ifges(6,4) * t169;
t204 = Ifges(6,1) * t173 - t263;
t205 = mrSges(6,1) * t169 + mrSges(6,2) * t173;
t276 = t173 / 0.2e1;
t277 = -t169 / 0.2e1;
t33 = -t170 * t71 + t174 * t91;
t28 = -pkin(4) * t156 - t33;
t119 = qJD(5) - t122;
t280 = t119 / 0.2e1;
t94 = t123 * t173 + t156 * t169;
t284 = t94 / 0.2e1;
t93 = -t123 * t169 + t156 * t173;
t286 = t93 / 0.2e1;
t275 = Ifges(6,4) * t94;
t31 = Ifges(6,2) * t93 + Ifges(6,6) * t119 + t275;
t92 = Ifges(6,4) * t93;
t32 = Ifges(6,1) * t94 + Ifges(6,5) * t119 + t92;
t305 = -t198 * mrSges(6,3) + t200 * t280 + t202 * t286 + t204 * t284 + t205 * t28 + t276 * t32 + t277 * t31;
t296 = t174 * t132 + t170 * t133;
t297 = -qJD(4) * t296 + t308 * t170 + t309 * t174;
t261 = Ifges(5,5) * t156;
t118 = Ifges(5,4) * t122;
t265 = Ifges(5,1) * t123;
t67 = t118 + t261 + t265;
t184 = t33 * mrSges(5,3) - t67 / 0.2e1 - t261 / 0.2e1 - t70 * mrSges(5,2);
t304 = t184 - t305;
t303 = t224 / 0.2e1;
t302 = -t160 * Ifges(4,6) / 0.2e1;
t131 = t161 + (-pkin(2) * t175 - pkin(3)) * t167;
t74 = pkin(4) * t143 - pkin(10) * t144 + t131;
t76 = -pkin(10) * t247 + t296;
t26 = -t169 * t76 + t173 * t74;
t301 = qJD(5) * t26 + t306 * t169 + t307 * t173;
t27 = t169 * t74 + t173 * t76;
t300 = -qJD(5) * t27 - t307 * t169 + t306 * t173;
t299 = -pkin(4) * t220 - t297;
t294 = t167 * t239 - t243;
t186 = t192 * qJD(2);
t193 = t152 * t245 + t249;
t216 = t168 * t220;
t54 = t193 * qJD(3) + (t166 * t186 + t216) * qJD(1);
t98 = t160 * t231 + (-t171 * t232 + t221) * t236;
t99 = t160 * t232 + (t171 * t231 + t222) * t236;
t21 = pkin(4) * t99 - pkin(10) * t98 + t54;
t103 = (t189 + t225) * t236;
t185 = t190 * qJD(2);
t226 = t168 * t247;
t215 = qJD(3) * t226;
t53 = (t152 * t244 - t141) * qJD(3) + (t166 * t185 + t215) * qJD(1);
t14 = t170 * t103 + t174 * t53 + t91 * t231 - t232 * t71;
t219 = qJD(2) * t234;
t214 = t171 * t219;
t7 = pkin(10) * t214 + t14;
t1 = qJD(5) * t12 + t169 * t21 + t173 * t7;
t2 = -qJD(5) * t13 - t169 * t7 + t173 * t21;
t293 = t1 * t173 - t169 * t2;
t15 = -qJD(4) * t34 + t103 * t174 - t170 * t53;
t292 = -t15 * mrSges(5,1) + t14 * mrSges(5,2) - Ifges(5,5) * t98 + Ifges(5,6) * t99;
t42 = qJD(5) * t93 + t169 * t214 + t173 * t98;
t43 = -qJD(5) * t94 - t169 * t98 + t173 * t214;
t11 = Ifges(6,1) * t42 + Ifges(6,4) * t43 + Ifges(6,5) * t99;
t291 = t11 / 0.2e1;
t290 = -t31 / 0.2e1;
t289 = t42 / 0.2e1;
t288 = t43 / 0.2e1;
t287 = -t93 / 0.2e1;
t285 = -t94 / 0.2e1;
t283 = t99 / 0.2e1;
t282 = -t118 / 0.2e1;
t281 = -t119 / 0.2e1;
t279 = -t143 / 0.2e1;
t278 = t144 / 0.2e1;
t39 = Ifges(6,5) * t42;
t274 = Ifges(6,5) * t94;
t38 = Ifges(6,6) * t43;
t273 = Ifges(6,6) * t93;
t272 = pkin(9) * t174;
t269 = t98 * Ifges(5,1);
t268 = t98 * Ifges(5,4);
t267 = t99 * Ifges(5,4);
t16 = -mrSges(6,1) * t43 + mrSges(6,2) * t42;
t78 = mrSges(5,1) * t214 - mrSges(5,3) * t98;
t266 = t16 - t78;
t264 = Ifges(5,4) * t123;
t260 = Ifges(5,2) * t122;
t259 = Ifges(5,6) * t156;
t258 = Ifges(6,3) * t119;
t104 = -t166 * t294 - t226;
t257 = t104 * t54;
t102 = mrSges(5,1) * t156 - mrSges(5,3) * t123;
t48 = -mrSges(6,1) * t93 + mrSges(6,2) * t94;
t251 = t48 - t102;
t137 = t209 * t236;
t61 = t170 * t137 + t174 * t84;
t250 = -mrSges(4,1) * t160 - mrSges(5,1) * t122 + mrSges(5,2) * t123 + mrSges(4,3) * t224;
t246 = t166 * t172;
t240 = t174 * t175;
t155 = -pkin(4) * t174 - pkin(10) * t170 - pkin(3);
t230 = qJD(5) * t155;
t229 = qJD(5) * t174;
t9 = Ifges(6,3) * t99 + t38 + t39;
t217 = t236 * t246;
t213 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t212 = -t54 * mrSges(4,1) - t53 * mrSges(4,2);
t208 = pkin(4) * t170 - pkin(10) * t174;
t207 = -mrSges(4,1) * t175 + mrSges(4,2) * t171;
t206 = mrSges(6,1) * t173 - mrSges(6,2) * t169;
t203 = Ifges(6,1) * t169 + t262;
t201 = Ifges(6,2) * t173 + t263;
t199 = Ifges(6,5) * t169 + Ifges(6,6) * t173;
t191 = t167 * t242 + t241;
t105 = t166 * t191 + t168 * t248;
t142 = -t165 * t166 * t176 + t167 * t168;
t73 = t105 * t174 + t142 * t170;
t40 = t104 * t173 - t169 * t73;
t41 = t104 * t169 + t173 * t73;
t60 = t137 * t174 - t170 * t84;
t72 = t105 * t170 - t142 * t174;
t86 = -t132 * t170 + t133 * t174;
t108 = -t144 * t169 - t173 * t247;
t194 = -t144 * t173 + t169 * t247;
t117 = -t152 * t165 + t159;
t158 = Ifges(4,4) * t223;
t183 = -t84 * mrSges(4,3) + Ifges(4,1) * t303 + t158 / 0.2e1 + t160 * Ifges(4,5) + t117 * mrSges(4,2);
t182 = t117 * mrSges(4,1) + t33 * mrSges(5,1) + t302 - (Ifges(4,4) * t171 + Ifges(4,2) * t175) * t236 / 0.2e1 + t156 * Ifges(5,3) + t123 * Ifges(5,5) + t122 * Ifges(5,6) - t34 * mrSges(5,2) - t85 * mrSges(4,3);
t30 = t258 + t273 + t274;
t66 = t259 + t260 + t264;
t181 = t259 / 0.2e1 - t273 / 0.2e1 - t274 / 0.2e1 - t258 / 0.2e1 - t70 * mrSges(5,1) + t13 * mrSges(6,2) - t12 * mrSges(6,1) + t34 * mrSges(5,3) - t30 / 0.2e1 + t66 / 0.2e1 + t264 / 0.2e1;
t180 = -t260 / 0.2e1 - t181;
t154 = Ifges(4,5) * t175 * t219;
t153 = Ifges(5,3) * t214;
t150 = t208 * qJD(4);
t136 = t207 * t236;
t135 = -mrSges(4,2) * t160 + mrSges(4,3) * t223;
t129 = t155 * t169 + t173 * t272;
t128 = t155 * t173 - t169 * t272;
t127 = (mrSges(4,1) * t171 + mrSges(4,2) * t175) * t219;
t116 = (t169 * t171 + t173 * t240) * t236;
t115 = (-t169 * t240 + t171 * t173) * t236;
t101 = -mrSges(5,2) * t156 + mrSges(5,3) * t122;
t83 = pkin(4) * t123 - pkin(10) * t122;
t82 = -t169 * t230 + t150 * t173 + (t169 * t232 - t173 * t229) * pkin(9);
t81 = t173 * t230 + t150 * t169 + (-t169 * t229 - t173 * t232) * pkin(9);
t79 = -mrSges(5,2) * t214 - mrSges(5,3) * t99;
t75 = pkin(4) * t247 - t86;
t69 = t215 + (qJD(3) * t294 + t185) * t166;
t68 = t216 + (qJD(3) * t191 + t186) * t166;
t64 = mrSges(6,1) * t119 - mrSges(6,3) * t94;
t63 = -mrSges(6,2) * t119 + mrSges(6,3) * t93;
t62 = (t171 * t237 + t208 * t235) * t165 + t193;
t57 = qJD(5) * t194 - t106 * t169 + t173 * t220;
t56 = qJD(5) * t108 + t106 * t173 + t169 * t220;
t55 = mrSges(5,1) * t99 + mrSges(5,2) * t98;
t51 = pkin(10) * t224 + t61;
t50 = -pkin(4) * t224 - t60;
t47 = Ifges(5,5) * t214 - t267 + t269;
t46 = -t99 * Ifges(5,2) + Ifges(5,6) * t214 + t268;
t25 = -mrSges(6,2) * t99 + mrSges(6,3) * t43;
t24 = mrSges(6,1) * t99 - mrSges(6,3) * t42;
t23 = -qJD(4) * t72 + t170 * t217 + t174 * t69;
t22 = qJD(4) * t73 + t170 * t69 - t174 * t217;
t20 = t169 * t62 + t173 * t51;
t19 = -t169 * t51 + t173 * t62;
t18 = t169 * t83 + t173 * t33;
t17 = -t169 * t33 + t173 * t83;
t10 = Ifges(6,4) * t42 + Ifges(6,2) * t43 + Ifges(6,6) * t99;
t8 = -pkin(4) * t214 - t15;
t6 = qJD(5) * t40 + t169 * t68 + t173 * t23;
t5 = -qJD(5) * t41 - t169 * t23 + t173 * t68;
t3 = [t23 * t101 + t104 * t55 + t142 * t127 + t69 * t135 + t40 * t24 + t41 * t25 + t5 * t64 + t6 * t63 + t73 * t79 + t266 * t72 + t250 * t68 + t251 * t22 + (-mrSges(3,1) * t172 - mrSges(3,2) * t176) * qJD(2) ^ 2 * t166 + (t136 * t246 + (t104 * t175 - t105 * t171) * qJD(3) * mrSges(4,3)) * t236 + m(4) * (t257 + t105 * t53 - t68 * t84 + t69 * t85 + (qJD(1) * t142 + t117) * t217) + m(5) * (t14 * t73 - t15 * t72 - t22 * t33 + t23 * t34 + t68 * t70 + t257) + m(6) * (t1 * t41 + t12 * t5 + t13 * t6 + t2 * t40 + t22 * t28 + t72 * t8); ((-m(4) * pkin(2) + t207) * t218 + ((Ifges(4,5) * t167 / 0.2e1 - t145 * mrSges(4,3) + 0.3e1 / 0.2e1 * Ifges(4,4) * t247) * t175 + (-Ifges(4,6) * t167 + Ifges(5,5) * t278 + Ifges(5,6) * t279 - t146 * mrSges(4,3) - 0.3e1 / 0.2e1 * Ifges(4,4) * t248 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - Ifges(5,3) / 0.2e1) * t247) * t171) * qJD(3)) * t236 + (t154 / 0.2e1 + t212) * t167 + (t131 * t54 + t14 * t296 + t15 * t86 + t295 * t70 + t297 * t33 + t298 * t34) * m(5) + t296 * t79 - t194 * t291 + (-Ifges(6,5) * t194 + Ifges(6,6) * t108 + Ifges(6,3) * t143) * t283 + (-Ifges(6,4) * t194 + Ifges(6,2) * t108 + Ifges(6,6) * t143) * t288 + (-Ifges(6,1) * t194 + Ifges(6,4) * t108 + Ifges(6,5) * t143) * t289 + t8 * (-mrSges(6,1) * t108 - mrSges(6,2) * t194) + t2 * (mrSges(6,1) * t143 + mrSges(6,3) * t194) + (Ifges(6,1) * t56 + Ifges(6,4) * t57 + Ifges(6,5) * t107) * t284 + (Ifges(6,4) * t56 + Ifges(6,2) * t57 + Ifges(6,6) * t107) * t286 + (-t136 * t225 + t54 * mrSges(4,3) * t171 - pkin(2) * t127 + (-t153 / 0.2e1 + t53 * mrSges(4,3) + t292) * t175 + (t183 * t175 + (t302 + t182) * t171) * qJD(3)) * t165 + t47 * t278 + t46 * t279 + (Ifges(6,5) * t56 + Ifges(6,6) * t57 + Ifges(6,3) * t107) * t280 - m(4) * (-t112 * t84 + t113 * t85 + t117 * t218) + t156 * (Ifges(5,5) * t106 - Ifges(5,6) * t107) / 0.2e1 + t300 * t64 + (t1 * t27 + t12 * t300 + t13 * t301 + t2 * t26 + t28 * t299 + t75 * t8) * m(6) + t301 * t63 + t297 * t102 + t298 * t101 + t299 * t48 + t295 * t250 + t122 * (Ifges(5,4) * t106 - Ifges(5,2) * t107) / 0.2e1 - t99 * (Ifges(5,4) * t144 - Ifges(5,2) * t143) / 0.2e1 + m(4) * (t139 * t85 - t140 * t84 - t145 * t54 + t146 * t53) + t123 * (Ifges(5,1) * t106 - Ifges(5,4) * t107) / 0.2e1 + t98 * (Ifges(5,1) * t144 - Ifges(5,4) * t143) / 0.2e1 + (-t106 * t33 - t107 * t34 - t14 * t143 - t144 * t15) * mrSges(5,3) + t26 * t24 + t27 * t25 + t56 * t32 / 0.2e1 + t57 * t31 / 0.2e1 + t28 * (-mrSges(6,1) * t57 + mrSges(6,2) * t56) + t75 * t16 + t86 * t78 + t106 * t67 / 0.2e1 + t107 * t30 / 0.2e1 + t12 * (mrSges(6,1) * t107 - mrSges(6,3) * t56) + t13 * (-mrSges(6,2) * t107 + mrSges(6,3) * t57) - t107 * t66 / 0.2e1 + t70 * (mrSges(5,1) * t107 + mrSges(5,2) * t106) + t108 * t10 / 0.2e1 + t131 * t55 + t143 * t9 / 0.2e1 + t1 * (-mrSges(6,2) * t143 + mrSges(6,3) * t108) + t54 * (mrSges(5,1) * t143 + mrSges(5,2) * t144) - t308 * t135; (t47 / 0.2e1 + t8 * t205 + t204 * t289 + t202 * t288 + t200 * t283 + t10 * t277 + t11 * t276 + t54 * mrSges(5,2) - t15 * mrSges(5,3) + t269 / 0.2e1 - t267 / 0.2e1 + (-t1 * t169 - t2 * t173) * mrSges(6,3) + (-m(5) * t15 + m(6) * t8 + t266) * pkin(9) + (t173 * t290 + t32 * t277 + t28 * t206 + t199 * t281 + t201 * t287 + t203 * t285 + (t12 * t169 - t13 * t173) * mrSges(6,3)) * qJD(5)) * t170 + (t14 * mrSges(5,3) - t9 / 0.2e1 + t46 / 0.2e1 - t54 * mrSges(5,1) + t268 / 0.2e1 - t39 / 0.2e1 - t38 / 0.2e1 + (-Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t99 + (m(5) * t14 + t79) * pkin(9) + t213) * t174 - t250 * t85 + (Ifges(6,1) * t116 + Ifges(6,4) * t115) * t285 + (Ifges(6,4) * t116 + Ifges(6,2) * t115) * t287 + t115 * t290 + ((qJD(3) * (Ifges(5,5) * t170 + Ifges(5,6) * t174) / 0.2e1 + Ifges(4,4) * t303 + (-qJD(3) + t160 / 0.2e1) * Ifges(4,6) - t182) * t171 + (-t158 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t224 + (-t265 / 0.2e1 + t282 + t184) * t174 - t180 * t170 - t183) * t175) * t236 + (Ifges(6,5) * t116 + Ifges(6,6) * t115) * t281 + (-t19 + t82) * t64 + t154 + (-t20 + t81) * t63 - m(6) * (t12 * t19 + t13 * t20 + t28 * t50) - m(5) * (t33 * t60 + t34 * t61 + t70 * t85) + ((t180 + (-m(5) * t34 - t101) * pkin(9)) * t170 + (t118 / 0.2e1 + t265 / 0.2e1 + (-m(5) * t33 + m(6) * t28 + t251) * pkin(9) - t304) * t174) * qJD(4) + (-t13 * t115 + t12 * t116) * mrSges(6,3) + (-t54 * m(5) - t55) * pkin(3) + m(6) * (t1 * t129 + t12 * t82 + t128 * t2 + t13 * t81) + t212 - t50 * t48 - t61 * t101 - t60 * t102 - t116 * t32 / 0.2e1 - t28 * (-mrSges(6,1) * t115 + mrSges(6,2) * t116) + t128 * t24 + t129 * t25 - t84 * t135; (t282 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t123 + t304) * t122 + t293 * mrSges(6,3) - t251 * t34 - t8 * t206 + t199 * t283 + t201 * t288 + t203 * t289 + t169 * t291 + t10 * t276 + t153 + (-pkin(4) * t8 - t12 * t17 - t13 * t18 - t28 * t34) * m(6) + t181 * t123 - t292 + t305 * qJD(5) + ((-m(6) * t198 - t169 * t63 - t173 * t64) * qJD(5) + m(6) * t293 - t169 * t24 + t173 * t25) * pkin(10) - pkin(4) * t16 - t18 * t63 - t17 * t64 - t33 * t101; -t28 * (mrSges(6,1) * t94 + mrSges(6,2) * t93) + (Ifges(6,1) * t93 - t275) * t285 + t31 * t284 + (Ifges(6,5) * t93 - Ifges(6,6) * t94) * t281 - t12 * t63 + t13 * t64 + (t12 * t93 + t13 * t94) * mrSges(6,3) - t213 + t9 + (-Ifges(6,2) * t94 + t32 + t92) * t287;];
tauc = t3(:);
