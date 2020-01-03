% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:54
% EndTime: 2019-12-31 20:17:05
% DurationCPUTime: 4.47s
% Computational Cost: add. (6888->407), mult. (18190->559), div. (0->0), fcn. (13439->8), ass. (0->198)
t177 = sin(qJ(2));
t180 = cos(qJ(2));
t243 = -qJ(3) - pkin(6);
t212 = qJD(2) * t243;
t136 = qJD(3) * t180 + t177 * t212;
t137 = -t177 * qJD(3) + t180 * t212;
t173 = sin(pkin(9));
t174 = cos(pkin(9));
t88 = t174 * t136 + t173 * t137;
t225 = t174 * t180;
t150 = -t173 * t177 + t225;
t139 = t150 * qJD(1);
t168 = -pkin(2) * t180 - pkin(1);
t223 = qJD(1) * t168;
t157 = qJD(3) + t223;
t106 = -t139 * pkin(3) + t157;
t172 = qJD(2) + qJD(4);
t175 = sin(qJ(5));
t178 = cos(qJ(5));
t176 = sin(qJ(4));
t179 = cos(qJ(4));
t161 = t243 * t180;
t156 = qJD(1) * t161;
t143 = t173 * t156;
t160 = t243 * t177;
t155 = qJD(1) * t160;
t149 = qJD(2) * pkin(2) + t155;
t100 = t174 * t149 + t143;
t221 = qJD(1) * t180;
t222 = qJD(1) * t177;
t140 = -t173 * t221 - t174 * t222;
t252 = pkin(7) * t140;
t74 = qJD(2) * pkin(3) + t100 + t252;
t226 = t174 * t156;
t101 = t173 * t149 - t226;
t253 = pkin(7) * t139;
t77 = t101 + t253;
t44 = t176 * t74 + t179 * t77;
t42 = pkin(8) * t172 + t44;
t196 = t139 * t176 - t179 * t140;
t210 = t179 * t139 + t140 * t176;
t45 = -pkin(4) * t210 - pkin(8) * t196 + t106;
t15 = -t175 * t42 + t178 * t45;
t16 = t175 * t45 + t178 * t42;
t199 = t15 * t178 + t16 * t175;
t191 = t199 * mrSges(6,3);
t244 = t196 * Ifges(5,1);
t91 = Ifges(5,4) * t210;
t273 = t244 / 0.2e1 + t91 / 0.2e1;
t201 = Ifges(6,5) * t178 - Ifges(6,6) * t175;
t239 = Ifges(6,4) * t178;
t203 = -Ifges(6,2) * t175 + t239;
t240 = Ifges(6,4) * t175;
t205 = Ifges(6,1) * t178 - t240;
t206 = mrSges(6,1) * t175 + mrSges(6,2) * t178;
t257 = t178 / 0.2e1;
t258 = -t175 / 0.2e1;
t82 = t172 * t175 + t178 * t196;
t267 = t82 / 0.2e1;
t275 = qJD(5) - t210;
t255 = Ifges(6,4) * t82;
t80 = t172 * t178 - t175 * t196;
t34 = t80 * Ifges(6,2) + Ifges(6,6) * t275 + t255;
t78 = Ifges(6,4) * t80;
t35 = t82 * Ifges(6,1) + Ifges(6,5) * t275 + t78;
t43 = -t176 * t77 + t179 * t74;
t41 = -pkin(4) * t172 - t43;
t279 = t34 * t258 + t35 * t257 + t41 * t206 + t80 * t203 / 0.2e1 + t205 * t267 + t275 * t201 / 0.2e1;
t280 = -t106 * mrSges(5,2) - t172 * Ifges(5,5) + t191 - t273 - t279;
t167 = pkin(2) * t174 + pkin(3);
t254 = pkin(2) * t173;
t135 = t176 * t167 + t179 * t254;
t104 = -t155 * t173 + t226;
t193 = t104 - t253;
t105 = t174 * t155 + t143;
t83 = t105 + t252;
t278 = t135 * qJD(4) - t176 * t83 + t179 * t193;
t277 = -t91 / 0.2e1 + t280;
t276 = -t15 * t175 + t16 * t178;
t142 = t150 * qJD(2);
t128 = qJD(1) * t142;
t151 = t173 * t180 + t174 * t177;
t152 = t173 * t160;
t79 = (-t151 * qJD(3) + (t225 * t243 - t152) * qJD(2)) * qJD(1);
t185 = -t128 * pkin(7) + t79;
t141 = t151 * qJD(2);
t127 = qJD(1) * t141;
t81 = t88 * qJD(1);
t69 = -pkin(7) * t127 + t81;
t13 = qJD(4) * t43 + t176 * t185 + t179 * t69;
t170 = pkin(2) * t222;
t166 = qJD(2) * t170;
t107 = pkin(3) * t127 + t166;
t62 = qJD(4) * t210 - t127 * t176 + t128 * t179;
t63 = qJD(4) * t196 + t179 * t127 + t128 * t176;
t25 = pkin(4) * t63 - pkin(8) * t62 + t107;
t2 = qJD(5) * t15 + t13 * t178 + t175 * t25;
t3 = -qJD(5) * t16 - t13 * t175 + t178 * t25;
t38 = qJD(5) * t80 + t178 * t62;
t39 = -qJD(5) * t82 - t175 * t62;
t274 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t38 + Ifges(6,6) * t39;
t65 = pkin(4) * t196 - pkin(8) * t210;
t272 = t38 / 0.2e1;
t271 = t39 / 0.2e1;
t270 = t63 / 0.2e1;
t269 = -t80 / 0.2e1;
t268 = -t82 / 0.2e1;
t265 = -t275 / 0.2e1;
t264 = pkin(1) * mrSges(3,1);
t263 = pkin(1) * mrSges(3,2);
t261 = -t140 / 0.2e1;
t260 = -t141 / 0.2e1;
t259 = t142 / 0.2e1;
t14 = qJD(4) * t44 + t176 * t69 - t179 * t185;
t108 = t174 * t160 + t161 * t173;
t89 = -pkin(7) * t151 + t108;
t109 = -t174 * t161 + t152;
t90 = pkin(7) * t150 + t109;
t54 = t176 * t90 - t179 * t89;
t251 = t14 * t54;
t250 = t178 * t2;
t249 = t3 * t175;
t245 = t210 * Ifges(5,2);
t242 = -mrSges(5,1) * t172 - mrSges(6,1) * t80 + mrSges(6,2) * t82 + mrSges(5,3) * t196;
t241 = Ifges(3,4) * t177;
t238 = t140 * Ifges(4,4);
t231 = Ifges(3,5) * qJD(2);
t230 = Ifges(3,6) * qJD(2);
t229 = qJD(2) * mrSges(3,1);
t228 = qJD(2) * mrSges(3,2);
t220 = qJD(2) * t177;
t219 = qJD(5) * t175;
t218 = qJD(5) * t178;
t215 = t231 / 0.2e1;
t214 = -t230 / 0.2e1;
t213 = t63 * mrSges(5,1) + t62 * mrSges(5,2);
t112 = -pkin(3) * t140 + t170;
t113 = pkin(2) * t220 + pkin(3) * t141;
t211 = t127 * mrSges(4,1) + t128 * mrSges(4,2);
t87 = -t136 * t173 + t174 * t137;
t208 = -t175 * t2 - t178 * t3;
t207 = mrSges(6,1) * t178 - mrSges(6,2) * t175;
t204 = Ifges(6,1) * t175 + t239;
t202 = Ifges(6,2) * t178 + t240;
t200 = Ifges(6,5) * t175 + Ifges(6,6) * t178;
t17 = mrSges(6,1) * t63 - mrSges(6,3) * t38;
t18 = -mrSges(6,2) * t63 + mrSges(6,3) * t39;
t197 = -t175 * t17 + t178 * t18;
t103 = t150 * t176 + t151 * t179;
t117 = -t150 * pkin(3) + t168;
t195 = t179 * t150 - t151 * t176;
t53 = -pkin(4) * t195 - t103 * pkin(8) + t117;
t55 = t176 * t89 + t179 * t90;
t26 = -t175 * t55 + t178 * t53;
t27 = t175 * t53 + t178 * t55;
t134 = t167 * t179 - t176 * t254;
t194 = -pkin(7) * t142 + t87;
t187 = -qJD(5) * t199 - t249;
t8 = t38 * Ifges(6,4) + t39 * Ifges(6,2) + t63 * Ifges(6,6);
t9 = t38 * Ifges(6,1) + t39 * Ifges(6,4) + t63 * Ifges(6,5);
t186 = -t13 * mrSges(5,2) + mrSges(6,3) * t250 + t175 * t9 / 0.2e1 + t204 * t272 + t202 * t271 + t200 * t270 - Ifges(5,6) * t63 + Ifges(5,5) * t62 + t8 * t257 + (-t207 - mrSges(5,1)) * t14 + t279 * qJD(5);
t184 = t16 * mrSges(6,2) - t275 * Ifges(6,3) - t82 * Ifges(6,5) - t80 * Ifges(6,6) + t172 * Ifges(5,6) + t245 / 0.2e1 + Ifges(5,4) * t196 - t106 * mrSges(5,1) - t15 * mrSges(6,1);
t183 = -t245 / 0.2e1 - t184;
t169 = Ifges(3,4) * t221;
t159 = mrSges(3,3) * t221 - t228;
t158 = -mrSges(3,3) * t222 + t229;
t148 = Ifges(3,1) * t222 + t169 + t231;
t147 = t230 + (t180 * Ifges(3,2) + t241) * qJD(1);
t133 = Ifges(4,4) * t139;
t132 = pkin(8) + t135;
t131 = -pkin(4) - t134;
t122 = t134 * qJD(4);
t116 = qJD(2) * mrSges(4,1) + t140 * mrSges(4,3);
t115 = -qJD(2) * mrSges(4,2) + t139 * mrSges(4,3);
t99 = -mrSges(4,1) * t139 - mrSges(4,2) * t140;
t94 = -t140 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t133;
t93 = t139 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t238;
t84 = -mrSges(5,2) * t172 + mrSges(5,3) * t210;
t73 = -pkin(7) * t141 + t88;
t67 = qJD(4) * t103 + t179 * t141 + t142 * t176;
t66 = qJD(4) * t195 - t141 * t176 + t142 * t179;
t64 = -mrSges(5,1) * t210 + mrSges(5,2) * t196;
t59 = Ifges(6,3) * t63;
t52 = mrSges(6,1) * t275 - mrSges(6,3) * t82;
t51 = -mrSges(6,2) * t275 + mrSges(6,3) * t80;
t48 = t112 + t65;
t47 = t176 * t193 + t179 * t83;
t28 = pkin(4) * t67 - pkin(8) * t66 + t113;
t24 = t175 * t65 + t178 * t43;
t23 = -t175 * t43 + t178 * t65;
t22 = qJD(4) * t55 + t176 * t73 - t179 * t194;
t21 = -qJD(4) * t54 + t176 * t194 + t179 * t73;
t20 = t175 * t48 + t178 * t47;
t19 = -t175 * t47 + t178 * t48;
t10 = -mrSges(6,1) * t39 + mrSges(6,2) * t38;
t5 = -qJD(5) * t27 - t175 * t21 + t178 * t28;
t4 = qJD(5) * t26 + t175 * t28 + t178 * t21;
t1 = [(-t43 * t66 - t44 * t67 + t54 * t62 - t55 * t63) * mrSges(5,3) + m(4) * (t100 * t87 + t101 * t88 + t108 * t79 + t109 * t81) + (t273 - t280) * t66 + t183 * t67 + (-Ifges(5,4) * t63 + Ifges(5,1) * t62 + t107 * mrSges(5,2) + t8 * t258 + t9 * t257 + t201 * t270 + t203 * t271 + t205 * t272 + (mrSges(5,3) + t206) * t14 + t208 * mrSges(6,3) + (t35 * t258 - t178 * t34 / 0.2e1 + t200 * t265 + t202 * t269 + t204 * t268 + t41 * t207 - t276 * mrSges(6,3)) * qJD(5)) * t103 + t157 * (mrSges(4,1) * t141 + mrSges(4,2) * t142) + (t128 * t151 + t142 * t261) * Ifges(4,1) + (Ifges(4,5) * t259 + Ifges(4,6) * t260 + (t148 / 0.2e1 - pkin(6) * t158 + t215 + (-0.2e1 * t263 + 0.3e1 / 0.2e1 * Ifges(3,4) * t180) * qJD(1)) * t180) * qJD(2) + (-t147 / 0.2e1 - pkin(6) * t159 + t214 + (-0.2e1 * t264 - 0.3e1 / 0.2e1 * t241 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t180) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t150 + mrSges(4,2) * t151) + m(4) * (t157 + t223) + t99) * pkin(2)) * t220 + t117 * t213 + m(6) * (t15 * t5 + t16 * t4 + t2 * t27 + t22 * t41 + t26 * t3 + t251) + m(5) * (t106 * t113 + t107 * t117 + t13 * t55 + t21 * t44 - t22 * t43 + t251) + t168 * t211 + t94 * t259 + t93 * t260 + t242 * t22 + t113 * t64 + t88 * t115 + t87 * t116 - (-Ifges(5,4) * t62 + t107 * mrSges(5,1) + t59 / 0.2e1 - t13 * mrSges(5,3) + (Ifges(5,2) + Ifges(6,3) / 0.2e1) * t63 + t274) * t195 + (-t100 * t142 - t101 * t141 - t108 * t128 - t109 * t127 + t150 * t81 - t151 * t79) * mrSges(4,3) + (-t127 * t151 + t150 * t128 + t139 * t259 - t141 * t261) * Ifges(4,4) + (-t150 * t127 + t139 * t260) * Ifges(4,2) + t21 * t84 + t4 * t51 + t5 * t52 + t54 * t10 + t26 * t17 + t27 * t18; -t157 * (-mrSges(4,1) * t140 + mrSges(4,2) * t139) + (-t175 * t52 + t178 * t51 + t84) * t122 + t186 + m(4) * (t173 * t81 + t174 * t79) * pkin(2) - (Ifges(4,2) * t140 + t133 + t94) * t139 / 0.2e1 - t183 * t196 + (-t134 * t62 - t135 * t63 + t196 * t44 + t210 * t43) * mrSges(5,3) + ((-t175 * t51 - t178 * t52) * t132 - t191) * qJD(5) + ((-t148 / 0.2e1 - t169 / 0.2e1 + t215 + qJD(1) * t263 + (t158 - t229) * pkin(6)) * t180 + (t147 / 0.2e1 + t214 + (t264 + t241 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t180) * qJD(1) + (t159 + t228) * pkin(6) + (-m(4) * t157 - t99) * pkin(2)) * t177) * qJD(1) - qJD(2) * (Ifges(4,5) * t139 + Ifges(4,6) * t140) / 0.2e1 - mrSges(6,3) * t249 + t93 * t261 + t197 * t132 + t140 * (Ifges(4,1) * t139 + t238) / 0.2e1 - m(4) * (t100 * t104 + t101 * t105) - Ifges(4,6) * t127 + Ifges(4,5) * t128 + t131 * t10 - t112 * t64 - t105 * t115 - t104 * t116 + (t100 * t139 - t101 * t140 + (-t127 * t173 - t128 * t174) * pkin(2)) * mrSges(4,3) - t47 * t84 + t79 * mrSges(4,1) - t81 * mrSges(4,2) - t20 * t51 - t19 * t52 + (-t244 / 0.2e1 + t277) * t210 + t278 * t242 + (-t15 * t19 - t16 * t20 + (t187 + t250) * t132 + t131 * t14 + t278 * t41 + t276 * t122) * m(6) + (-t106 * t112 + t13 * t135 - t134 * t14 + (t122 - t47) * t44 - t278 * t43) * m(5); -t139 * t115 - t140 * t116 - t210 * t84 - t242 * t196 + (t275 * t51 + t17) * t178 + (-t275 * t52 + t18) * t175 + t211 + t213 + (-t41 * t196 + t275 * t276 - t208) * m(6) + (t196 * t43 - t210 * t44 + t107) * m(5) + (-t100 * t140 - t101 * t139 + t166) * m(4); (-t52 * t218 - t51 * t219 + t197) * pkin(8) + (t43 * mrSges(5,3) + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t196 + t277) * t210 + (t44 * mrSges(5,3) + t184) * t196 + t186 + t187 * mrSges(6,3) - t242 * t44 - t43 * t84 - t24 * t51 - t23 * t52 - pkin(4) * t10 + ((-t15 * t218 - t16 * t219 - t249 + t250) * pkin(8) - pkin(4) * t14 - t15 * t23 - t16 * t24 - t41 * t44) * m(6); t59 - t41 * (mrSges(6,1) * t82 + mrSges(6,2) * t80) + (Ifges(6,1) * t80 - t255) * t268 + t34 * t267 + (Ifges(6,5) * t80 - Ifges(6,6) * t82) * t265 - t15 * t51 + t16 * t52 + (t15 * t80 + t16 * t82) * mrSges(6,3) + (-Ifges(6,2) * t82 + t35 + t78) * t269 + t274;];
tauc = t1(:);
