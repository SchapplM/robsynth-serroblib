% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:03
% EndTime: 2019-12-05 17:19:24
% DurationCPUTime: 6.34s
% Computational Cost: add. (4078->449), mult. (10473->649), div. (0->0), fcn. (7288->10), ass. (0->224)
t158 = sin(qJ(4));
t160 = sin(qJ(2));
t162 = cos(qJ(4));
t155 = sin(pkin(5));
t211 = qJD(1) * t155;
t163 = cos(qJ(3));
t164 = cos(qJ(2));
t213 = t163 * t164;
t100 = (t158 * t160 + t162 * t213) * t211;
t159 = sin(qJ(3));
t186 = pkin(3) * t159 - pkin(8) * t163;
t137 = t186 * qJD(3);
t141 = -pkin(3) * t163 - pkin(8) * t159 - pkin(2);
t202 = qJD(4) * t162;
t203 = qJD(4) * t158;
t205 = qJD(3) * t162;
t56 = t158 * t137 + t141 * t202 + (-t159 * t205 - t163 * t203) * pkin(7);
t294 = t56 - t100;
t214 = t162 * t163;
t151 = pkin(7) * t214;
t176 = pkin(4) * t159 - pkin(9) * t214;
t206 = qJD(3) * t159;
t245 = pkin(7) * t158;
t212 = t162 * t137 + t206 * t245;
t244 = pkin(9) * t159;
t99 = (-t158 * t213 + t160 * t162) * t211;
t293 = t176 * qJD(3) + (-t151 + (-t141 + t244) * t158) * qJD(4) + t212 - t99;
t204 = qJD(3) * t163;
t170 = t158 * t204 + t159 * t202;
t292 = pkin(9) * t170 - t294;
t260 = -pkin(9) - pkin(8);
t197 = qJD(4) * t260;
t207 = qJD(2) * t163;
t196 = t160 * t211;
t139 = qJD(2) * pkin(7) + t196;
t156 = cos(pkin(5));
t210 = qJD(1) * t156;
t104 = -t159 * t139 + t163 * t210;
t134 = t186 * qJD(2);
t61 = t162 * t104 + t158 * t134;
t291 = -t61 + (pkin(9) * t207 + t197) * t158;
t60 = -t104 * t158 + t162 * t134;
t290 = -qJD(2) * t176 + t162 * t197 - t60;
t289 = -Ifges(4,1) / 0.2e1;
t154 = Ifges(4,4) * t207;
t288 = -t154 / 0.2e1;
t287 = qJD(3) / 0.2e1;
t161 = cos(qJ(5));
t157 = sin(qJ(5));
t209 = qJD(2) * t159;
t129 = -t158 * t209 + t205;
t195 = t164 * t211;
t107 = qJD(2) * t141 - t195;
t194 = t159 * t210;
t105 = t139 * t163 + t194;
t91 = qJD(3) * pkin(8) + t105;
t51 = t107 * t158 + t162 * t91;
t34 = pkin(9) * t129 + t51;
t226 = t157 * t34;
t150 = qJD(4) - t207;
t130 = qJD(3) * t158 + t162 * t209;
t50 = t162 * t107 - t158 * t91;
t33 = -pkin(9) * t130 + t50;
t30 = pkin(4) * t150 + t33;
t12 = t161 * t30 - t226;
t224 = t161 * t34;
t13 = t157 * t30 + t224;
t201 = qJD(2) * qJD(3);
t189 = t159 * t201;
t148 = Ifges(6,3) * t189;
t187 = t161 * t129 - t130 * t157;
t75 = t129 * t157 + t130 * t161;
t247 = Ifges(6,4) * t75;
t146 = qJD(5) + t150;
t252 = -t146 / 0.2e1;
t264 = -t75 / 0.2e1;
t266 = -t187 / 0.2e1;
t69 = Ifges(6,4) * t187;
t29 = Ifges(6,1) * t75 + Ifges(6,5) * t146 + t69;
t90 = -qJD(3) * pkin(3) - t104;
t63 = -pkin(4) * t129 + t90;
t286 = t148 + (Ifges(6,5) * t187 - Ifges(6,6) * t75) * t252 + (t12 * t187 + t13 * t75) * mrSges(6,3) + (-Ifges(6,2) * t75 + t29 + t69) * t266 - t63 * (mrSges(6,1) * t75 + mrSges(6,2) * t187) + (Ifges(6,1) * t187 - t247) * t264;
t191 = Ifges(4,5) * t287;
t128 = t162 * t141;
t76 = -t162 * t244 + t128 + (-pkin(4) - t245) * t163;
t109 = t158 * t141 + t151;
t85 = -t158 * t244 + t109;
t37 = -t157 * t85 + t161 * t76;
t285 = qJD(5) * t37 + t293 * t157 - t292 * t161;
t38 = t157 * t76 + t161 * t85;
t284 = -qJD(5) * t38 + t292 * t157 + t293 * t161;
t144 = t260 * t158;
t145 = t260 * t162;
t96 = t144 * t157 - t145 * t161;
t283 = -qJD(5) * t96 - t291 * t157 + t290 * t161;
t95 = t144 * t161 + t145 * t157;
t282 = qJD(5) * t95 + t290 * t157 + t291 * t161;
t234 = qJD(2) * pkin(2);
t140 = -t195 - t234;
t178 = t158 * t51 + t162 * t50;
t237 = Ifges(5,4) * t162;
t181 = -Ifges(5,2) * t158 + t237;
t238 = Ifges(5,4) * t158;
t183 = Ifges(5,1) * t162 - t238;
t184 = mrSges(5,1) * t158 + mrSges(5,2) * t162;
t235 = Ifges(5,6) * t158;
t236 = Ifges(5,5) * t162;
t248 = t162 / 0.2e1;
t249 = -t158 / 0.2e1;
t253 = t130 / 0.2e1;
t239 = Ifges(5,4) * t130;
t65 = Ifges(5,2) * t129 + Ifges(5,6) * t150 + t239;
t126 = Ifges(5,4) * t129;
t66 = Ifges(5,1) * t130 + Ifges(5,5) * t150 + t126;
t166 = -t178 * mrSges(5,3) + t65 * t249 + t66 * t248 + t90 * t184 + t129 * t181 / 0.2e1 + t183 * t253 + t150 * (-t235 + t236) / 0.2e1;
t281 = -t140 * mrSges(4,2) + t104 * mrSges(4,3) + t209 * t289 - t166 - t191 + t288;
t101 = (t137 + t196) * qJD(2);
t215 = t155 * t164;
t192 = qJD(2) * t215;
t169 = qJD(1) * (qJD(3) * t156 + t192);
t67 = -t139 * t206 + t163 * t169;
t18 = -qJD(4) * t51 + t162 * t101 - t158 * t67;
t200 = qJD(3) * qJD(4);
t97 = t162 * t200 + (-t159 * t203 + t162 * t204) * qJD(2);
t11 = pkin(4) * t189 - pkin(9) * t97 + t18;
t17 = t158 * t101 + t107 * t202 + t162 * t67 - t203 * t91;
t98 = -qJD(2) * t170 - t158 * t200;
t14 = pkin(9) * t98 + t17;
t2 = qJD(5) * t12 + t11 * t157 + t14 * t161;
t25 = qJD(5) * t187 + t157 * t98 + t161 * t97;
t26 = -qJD(5) * t75 - t157 * t97 + t161 * t98;
t3 = -qJD(5) * t13 + t11 * t161 - t14 * t157;
t280 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t25 + Ifges(6,6) * t26;
t28 = Ifges(6,2) * t187 + Ifges(6,6) * t146 + t247;
t278 = t28 / 0.2e1;
t190 = -Ifges(4,6) * qJD(3) / 0.2e1;
t35 = -mrSges(6,1) * t187 + mrSges(6,2) * t75;
t198 = m(6) * t63 + t35;
t177 = t157 * t158 - t161 * t162;
t113 = t177 * t159;
t273 = -t158 * t18 + t162 * t17;
t272 = qJD(4) + qJD(5);
t269 = -t18 * mrSges(5,1) + t17 * mrSges(5,2) - Ifges(5,5) * t97 - Ifges(5,6) * t98 - t280;
t268 = t25 / 0.2e1;
t267 = t26 / 0.2e1;
t265 = t187 / 0.2e1;
t263 = t75 / 0.2e1;
t262 = t97 / 0.2e1;
t261 = t98 / 0.2e1;
t132 = t157 * t162 + t158 * t161;
t112 = t132 * t159;
t257 = -t112 / 0.2e1;
t256 = -t113 / 0.2e1;
t255 = -t129 / 0.2e1;
t254 = -t130 / 0.2e1;
t251 = t146 / 0.2e1;
t250 = -t150 / 0.2e1;
t246 = pkin(4) * t158;
t240 = Ifges(4,4) * t159;
t216 = t155 * t160;
t115 = -t156 * t163 + t159 * t216;
t68 = t139 * t204 + t159 * t169;
t232 = t115 * t68;
t173 = t132 * t163;
t110 = qJD(2) * t173;
t79 = t272 * t132;
t222 = -t110 + t79;
t172 = t177 * t163;
t111 = qJD(2) * t172;
t78 = t272 * t177;
t221 = -t111 + t78;
t220 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t129 + mrSges(5,2) * t130 + mrSges(4,3) * t209;
t208 = qJD(2) * t160;
t199 = mrSges(4,3) * t207;
t193 = t155 * t208;
t143 = -qJD(3) * mrSges(4,2) + t199;
t188 = -m(4) * t105 - t143;
t185 = mrSges(5,1) * t162 - mrSges(5,2) * t158;
t182 = Ifges(5,1) * t158 + t237;
t180 = Ifges(5,2) * t162 + t238;
t179 = Ifges(5,5) * t158 + Ifges(5,6) * t162;
t116 = t156 * t159 + t163 * t216;
t175 = -t116 * t162 + t158 * t215;
t83 = -t116 * t158 - t162 * t215;
t39 = t157 * t175 + t161 * t83;
t40 = t157 * t83 - t161 * t175;
t174 = -m(4) * t104 + m(5) * t90 + t220;
t167 = t12 * mrSges(6,1) + t140 * mrSges(4,1) + t50 * mrSges(5,1) + t190 - (t163 * Ifges(4,2) + t240) * qJD(2) / 0.2e1 + t146 * Ifges(6,3) + t75 * Ifges(6,5) + t187 * Ifges(6,6) + t150 * Ifges(5,3) + t130 * Ifges(5,5) + t129 * Ifges(5,6) - t105 * mrSges(4,3) - t13 * mrSges(6,2) - t51 * mrSges(5,2);
t165 = qJD(2) ^ 2;
t153 = -pkin(4) * t162 - pkin(3);
t149 = Ifges(5,3) * t189;
t138 = (pkin(7) + t246) * t159;
t133 = (-mrSges(4,1) * t163 + mrSges(4,2) * t159) * qJD(2);
t123 = (mrSges(4,1) * t159 + mrSges(4,2) * t163) * t201;
t108 = -t163 * t245 + t128;
t106 = pkin(4) * t170 + pkin(7) * t204;
t103 = mrSges(5,1) * t150 - mrSges(5,3) * t130;
t102 = -mrSges(5,2) * t150 + mrSges(5,3) * t129;
t82 = -qJD(3) * t115 + t163 * t192;
t81 = qJD(3) * t116 + t159 * t192;
t77 = t194 + (qJD(2) * t246 + t139) * t163;
t71 = -mrSges(5,2) * t189 + mrSges(5,3) * t98;
t70 = mrSges(5,1) * t189 - mrSges(5,3) * t97;
t59 = mrSges(6,1) * t146 - mrSges(6,3) * t75;
t58 = -mrSges(6,2) * t146 + mrSges(6,3) * t187;
t57 = -qJD(4) * t109 + t212;
t54 = -mrSges(5,1) * t98 + mrSges(5,2) * t97;
t46 = t97 * Ifges(5,1) + t98 * Ifges(5,4) + Ifges(5,5) * t189;
t45 = t97 * Ifges(5,4) + t98 * Ifges(5,2) + Ifges(5,6) * t189;
t44 = -qJD(3) * t173 + t113 * t272;
t43 = -qJD(3) * t172 - t159 * t79;
t42 = -pkin(4) * t98 + t68;
t32 = qJD(4) * t83 + t158 * t193 + t162 * t82;
t31 = qJD(4) * t175 - t158 * t82 + t162 * t193;
t22 = -mrSges(6,2) * t189 + mrSges(6,3) * t26;
t21 = mrSges(6,1) * t189 - mrSges(6,3) * t25;
t16 = t161 * t33 - t226;
t15 = -t157 * t33 - t224;
t10 = -mrSges(6,1) * t26 + mrSges(6,2) * t25;
t9 = t25 * Ifges(6,1) + t26 * Ifges(6,4) + Ifges(6,5) * t189;
t8 = t25 * Ifges(6,4) + t26 * Ifges(6,2) + Ifges(6,6) * t189;
t5 = -qJD(5) * t40 - t157 * t32 + t161 * t31;
t4 = qJD(5) * t39 + t157 * t31 + t161 * t32;
t1 = [-t116 * mrSges(4,3) * t189 + t32 * t102 + t31 * t103 + t82 * t143 + t39 * t21 + t40 * t22 + t4 * t58 + t5 * t59 + t83 * t70 - t175 * t71 + (t35 + t220) * t81 + (qJD(3) * t199 + t10 + t54) * t115 + ((-mrSges(3,2) * t165 - t123) * t164 + (-mrSges(3,1) * t165 + qJD(2) * t133) * t160) * t155 + m(4) * (-t104 * t81 + t105 * t82 + t232 + t116 * t67 + (t140 - t195) * t193) + m(5) * (-t17 * t175 + t18 * t83 + t31 * t50 + t32 * t51 + t81 * t90 + t232) + m(6) * (t115 * t42 + t12 * t5 + t13 * t4 + t2 * t40 + t3 * t39 + t63 * t81); (t57 - t99) * t103 + t44 * t278 + t294 * t102 + (0.2e1 * (-t234 / 0.2e1 - t140 / 0.2e1) * m(4) - t133) * t196 + (-t112 * t2 + t113 * t3 - t12 * t43 + t13 * t44) * mrSges(6,3) + t42 * (mrSges(6,1) * t112 - mrSges(6,2) * t113) + (-Ifges(6,4) * t113 - Ifges(6,2) * t112) * t267 + (-Ifges(6,1) * t113 - Ifges(6,4) * t112) * t268 + t284 * t59 + t285 * t58 + (t106 * t63 + t284 * t12 + t285 * t13 + t138 * t42 + t2 * t38 + t3 * t37) * m(6) - m(5) * (t100 * t51 + t50 * t99) + m(5) * (t108 * t18 + t109 * t17 + t50 * t57 + t51 * t56) + (t183 * t262 + t181 * t261 + t45 * t249 + t46 * t248 + (-t158 * t17 - t162 * t18) * mrSges(5,3) + (mrSges(4,3) + t184) * t68 + (mrSges(4,2) * t208 + (-t174 - t198) * t164) * t211 + (t90 * t185 + t180 * t255 + t182 * t254 + t179 * t250 - t162 * t65 / 0.2e1 + t66 * t249 + (t158 * t50 - t162 * t51) * mrSges(5,3)) * qJD(4) + (((-0.3e1 / 0.2e1 * Ifges(4,4) + t236 / 0.2e1 - t235 / 0.2e1) * t159 + Ifges(6,5) * t256 + Ifges(6,6) * t257 + (0.3e1 / 0.2e1 * Ifges(4,1) - Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,2) - Ifges(5,3) / 0.2e1) * t163) * qJD(2) + t167 + t190) * qJD(3) + (t54 + (m(5) + m(4)) * t68 + t188 * qJD(3)) * pkin(7)) * t159 + (-t148 / 0.2e1 - t149 / 0.2e1 + (m(4) * pkin(7) + mrSges(4,3)) * t67 + (-mrSges(4,1) * t208 + t164 * t188) * t211 + (0.3e1 / 0.2e1 * t154 + t191 + t174 * pkin(7) - t281) * qJD(3) + t269) * t163 + t37 * t21 + t38 * t22 + t43 * t29 / 0.2e1 + t63 * (-mrSges(6,1) * t44 + mrSges(6,2) * t43) + t106 * t35 + t108 * t70 + t109 * t71 - pkin(2) * t123 + t138 * t10 + (Ifges(6,5) * t43 + Ifges(6,6) * t44) * t251 + t9 * t256 + t8 * t257 + (Ifges(6,1) * t43 + Ifges(6,4) * t44) * t263 + (Ifges(6,4) * t43 + Ifges(6,2) * t44) * t265; m(5) * (-pkin(3) * t68 + pkin(8) * t273) + t273 * mrSges(5,3) + t42 * (mrSges(6,1) * t177 + mrSges(6,2) * t132) + (t12 * t221 - t13 * t222 - t132 * t3 - t177 * t2) * mrSges(6,3) + (Ifges(6,4) * t132 - Ifges(6,2) * t177) * t267 + (Ifges(6,1) * t132 - Ifges(6,4) * t177) * t268 - t177 * t8 / 0.2e1 + (-Ifges(6,5) * t78 - Ifges(6,6) * t79) * t251 + (-Ifges(6,1) * t78 - Ifges(6,4) * t79) * t263 + (-Ifges(6,4) * t78 - Ifges(6,2) * t79) * t265 + (-t79 / 0.2e1 + t110 / 0.2e1) * t28 + (-t78 / 0.2e1 + t111 / 0.2e1) * t29 + (-Ifges(6,5) * t111 - Ifges(6,6) * t110) * t252 + (-Ifges(6,1) * t111 - Ifges(6,4) * t110) * t264 + (-Ifges(6,4) * t111 - Ifges(6,2) * t110) * t266 + (-t158 * t70 + t162 * t71) * pkin(8) - m(5) * (t105 * t90 + t50 * t60 + t51 * t61) + t282 * t58 + t283 * t59 + (t283 * t12 + t282 * t13 + t153 * t42 + t2 * t96 + t3 * t95 - t63 * t77) * m(6) + ((t191 + t288 + t281) * t163 + ((t240 / 0.2e1 + (Ifges(4,2) / 0.2e1 + t289) * t163) * qJD(2) - t167 + t190 + (Ifges(6,5) * t132 - Ifges(6,6) * t177 + t179) * t287) * t159) * qJD(2) + (-mrSges(4,1) - t185) * t68 - pkin(3) * t54 - t67 * mrSges(4,2) - t77 * t35 + t95 * t21 + t96 * t22 - t61 * t102 - t60 * t103 + t132 * t9 / 0.2e1 - t104 * t143 + t153 * t10 + t158 * t46 / 0.2e1 - t220 * t105 + (mrSges(6,1) * t222 - mrSges(6,2) * t221) * t63 + (t198 * t246 + (-m(5) * t178 - t158 * t102 - t162 * t103) * pkin(8) + t166) * qJD(4) + t45 * t248 + t180 * t261 + t182 * t262; (-Ifges(5,2) * t130 + t126 + t66) * t255 - t269 - m(6) * (t12 * t15 + t13 * t16) + t75 * t278 + (t129 * t50 + t130 * t51) * mrSges(5,3) + (t157 * t22 + t161 * t21 + m(6) * (t157 * t2 + t161 * t3) - t198 * t130 + (-t157 * t59 + t161 * t58 + m(6) * (-t12 * t157 + t13 * t161)) * qJD(5)) * pkin(4) + t149 - t16 * t58 - t15 * t59 - t50 * t102 + t51 * t103 - t90 * (mrSges(5,1) * t130 + mrSges(5,2) * t129) + (Ifges(5,5) * t129 - Ifges(5,6) * t130) * t250 + t65 * t253 + (Ifges(5,1) * t129 - t239) * t254 + t286; -t12 * t58 + t13 * t59 + t28 * t263 + t280 + t286;];
tauc = t1(:);
