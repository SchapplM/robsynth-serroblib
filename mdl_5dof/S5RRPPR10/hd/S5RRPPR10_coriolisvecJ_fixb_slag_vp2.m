% Calculate vector of centrifugal and Coriolis load on the joints for
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:02
% EndTime: 2019-12-31 19:43:15
% DurationCPUTime: 5.25s
% Computational Cost: add. (2601->445), mult. (6803->618), div. (0->0), fcn. (4191->6), ass. (0->216)
t279 = Ifges(4,1) + Ifges(5,1);
t275 = Ifges(5,4) + Ifges(4,5);
t278 = -qJD(2) / 0.2e1;
t277 = qJD(2) / 0.2e1;
t276 = -Ifges(4,4) + Ifges(5,5);
t154 = sin(qJ(5));
t156 = cos(qJ(5));
t152 = sin(pkin(8));
t153 = cos(pkin(8));
t173 = t152 * t154 + t153 * t156;
t155 = sin(qJ(2));
t198 = -pkin(6) * t152 - pkin(3);
t157 = cos(qJ(2));
t212 = t153 * t157;
t203 = pkin(7) * t212;
t163 = -t203 + (-pkin(4) + t198) * t155;
t179 = pkin(2) * t155 - qJ(3) * t157;
t125 = t179 * qJD(1);
t214 = t153 * t125;
t38 = qJD(1) * t163 - t214;
t110 = t152 * t125;
t210 = qJD(1) * t155;
t144 = qJ(4) * t210;
t213 = t153 * t155;
t216 = t152 * t157;
t170 = -pkin(6) * t213 + pkin(7) * t216;
t48 = qJD(1) * t170 + t110 + t144;
t237 = -pkin(7) + qJ(3);
t131 = t237 * t152;
t132 = t237 * t153;
t78 = t131 * t156 - t132 * t154;
t274 = qJD(3) * t173 + qJD(5) * t78 - t154 * t38 - t156 * t48;
t122 = t152 * t156 - t153 * t154;
t79 = t131 * t154 + t132 * t156;
t273 = qJD(3) * t122 - qJD(5) * t79 + t154 * t48 - t156 * t38;
t227 = Ifges(5,5) * t152;
t231 = Ifges(4,4) * t152;
t272 = t279 * t153 + t227 - t231;
t168 = t173 * t157;
t164 = qJD(2) * t168;
t197 = t152 * t210;
t119 = -t153 * qJD(2) + t197;
t196 = t153 * t210;
t205 = t152 * qJD(2);
t120 = t196 + t205;
t60 = t119 * t156 - t120 * t154;
t25 = qJD(1) * t164 + qJD(5) * t60;
t258 = t25 / 0.2e1;
t169 = t122 * t157;
t165 = qJD(2) * t169;
t174 = t119 * t154 + t120 * t156;
t26 = qJD(1) * t165 - qJD(5) * t174;
t257 = t26 / 0.2e1;
t270 = t120 / 0.2e1;
t204 = qJD(1) * qJD(2);
t193 = t155 * t204;
t269 = -t193 / 0.2e1;
t268 = t204 / 0.2e1;
t194 = Ifges(3,6) * t278;
t266 = Ifges(3,5) * t277;
t252 = pkin(3) + pkin(4);
t265 = t152 * t252;
t206 = qJD(4) * t157;
t208 = qJD(2) * t155;
t264 = qJ(4) * t208 - t206;
t209 = qJD(1) * t157;
t143 = qJD(5) + t209;
t147 = pkin(6) * t210;
t126 = (qJD(3) - t147) * qJD(2);
t107 = qJD(2) * t179 - t155 * qJD(3);
t95 = t107 * qJD(1);
t44 = -t152 * t126 + t153 * t95;
t20 = (-t155 * t252 - t203) * t204 - t44;
t45 = t153 * t126 + t152 * t95;
t199 = qJ(4) * t193 + t45;
t21 = (pkin(7) * t205 - qJD(4)) * t209 + t199;
t129 = -pkin(2) * t157 - t155 * qJ(3) - pkin(1);
t115 = t129 * qJD(1);
t148 = pkin(6) * t209;
t138 = qJD(2) * qJ(3) + t148;
t64 = t153 * t115 - t152 * t138;
t47 = pkin(3) * t209 + qJD(4) - t64;
t24 = pkin(4) * t209 - t120 * pkin(7) + t47;
t65 = t152 * t115 + t153 * t138;
t49 = -qJ(4) * t209 + t65;
t29 = t119 * pkin(7) + t49;
t8 = -t154 * t29 + t156 * t24;
t1 = qJD(5) * t8 + t154 * t20 + t156 * t21;
t9 = t154 * t24 + t156 * t29;
t2 = -qJD(5) * t9 - t154 * t21 + t156 * t20;
t263 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t25 + Ifges(6,6) * t26;
t232 = Ifges(3,4) * t155;
t262 = -0.2e1 * (Ifges(4,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t119 + (Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t120 + t49 * mrSges(5,3) + t64 * mrSges(4,1) + t9 * mrSges(6,2) + t194 - (t157 * Ifges(3,2) + t232) * qJD(1) / 0.2e1 - t143 * Ifges(6,3) - t174 * Ifges(6,5) - t60 * Ifges(6,6) - t47 * mrSges(5,1) - t65 * mrSges(4,2) - t8 * mrSges(6,1) + t275 * t270 - (Ifges(4,3) + Ifges(5,2)) * t209 / 0.2e1;
t146 = Ifges(3,4) * t209;
t226 = Ifges(5,5) * t153;
t180 = Ifges(5,3) * t152 + t226;
t230 = Ifges(4,4) * t153;
t181 = -Ifges(4,2) * t152 + t230;
t184 = qJD(2) * pkin(2) - qJD(3) - t147;
t233 = mrSges(5,3) * t153;
t244 = t153 / 0.2e1;
t245 = t152 / 0.2e1;
t246 = -t152 / 0.2e1;
t166 = qJ(4) * t120 + t184;
t42 = pkin(3) * t119 - t166;
t261 = -(t49 * t152 - t47 * t153) * mrSges(5,2) - (t65 * t152 + t64 * t153) * mrSges(4,3) + (-m(4) * t184 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t119 + mrSges(4,2) * t120 + mrSges(3,3) * t210) * pkin(6) - t184 * (mrSges(4,1) * t152 + mrSges(4,2) * t153) + t42 * (mrSges(5,1) * t152 - t233) + Ifges(3,1) * t210 / 0.2e1 + t146 / 0.2e1 + t266 + (Ifges(5,5) * t120 - Ifges(5,6) * t209) * t245 + (Ifges(4,4) * t120 - Ifges(4,6) * t209) * t246 + t272 * t270 + (t279 * t120 - t275 * t209) * t244 + (Ifges(5,3) * t245 - Ifges(4,2) * t246 + t276 * t244 - t181 / 0.2e1 + t180 / 0.2e1) * t119;
t260 = Ifges(6,4) * t258 + Ifges(6,2) * t257 + Ifges(6,6) * t269;
t259 = Ifges(6,1) * t258 + Ifges(6,4) * t257 + Ifges(6,5) * t269;
t256 = -t60 / 0.2e1;
t255 = t60 / 0.2e1;
t254 = -t174 / 0.2e1;
t253 = t174 / 0.2e1;
t251 = pkin(1) * mrSges(3,1);
t250 = pkin(1) * mrSges(3,2);
t248 = -t143 / 0.2e1;
t247 = t143 / 0.2e1;
t243 = Ifges(6,4) * t174;
t242 = pkin(3) * t152;
t241 = pkin(7) * t155;
t17 = -mrSges(6,1) * t60 + mrSges(6,2) * t174;
t69 = mrSges(5,1) * t119 - mrSges(5,3) * t120;
t236 = -t69 + t17;
t89 = -t119 * mrSges(5,2) - mrSges(5,3) * t209;
t90 = mrSges(4,2) * t209 - t119 * mrSges(4,3);
t235 = t89 + t90;
t91 = -mrSges(4,1) * t209 - t120 * mrSges(4,3);
t92 = mrSges(5,1) * t209 + t120 * mrSges(5,2);
t234 = t91 - t92;
t229 = Ifges(5,4) * t153;
t228 = Ifges(4,5) * t153;
t225 = Ifges(4,6) * t152;
t224 = Ifges(5,6) * t152;
t108 = t173 * qJD(5);
t94 = qJD(1) * t168;
t222 = -t108 - t94;
t109 = t122 * qJD(5);
t93 = qJD(1) * t169;
t221 = -t109 - t93;
t218 = qJ(4) * t153;
t217 = qJD(2) * mrSges(3,2);
t215 = t153 * t107;
t192 = t157 * t204;
t188 = t153 * t192;
t211 = -qJ(4) * t188 - t120 * qJD(4);
t189 = t152 * t192;
t97 = mrSges(4,1) * t189 + mrSges(4,2) * t188;
t88 = pkin(6) * t212 + t152 * t129;
t207 = qJD(2) * t157;
t202 = pkin(6) * t208;
t195 = qJD(4) * t213;
t191 = qJ(4) * t152 + pkin(2);
t141 = pkin(6) * t216;
t87 = t153 * t129 - t141;
t186 = t198 * t155;
t7 = -t26 * mrSges(6,1) + t25 * mrSges(6,2);
t81 = -qJ(4) * t157 + t88;
t178 = -t218 + t242;
t36 = -mrSges(6,2) * t143 + mrSges(6,3) * t60;
t37 = mrSges(6,1) * t143 - mrSges(6,3) * t174;
t175 = -t154 * t37 + t156 * t36;
t151 = t157 * pkin(3);
t50 = pkin(4) * t157 + t141 + t151 + (-t129 - t241) * t153;
t63 = t152 * t241 + t81;
t15 = -t154 * t63 + t156 * t50;
t16 = t154 * t50 + t156 * t63;
t98 = t152 * t107;
t73 = -t153 * t202 + t98;
t84 = -pkin(6) * t196 + t110;
t172 = pkin(6) + t178;
t171 = t218 - t265;
t101 = t122 * t155;
t167 = -pkin(6) + t171;
t28 = -qJD(1) * t206 + t199;
t39 = (pkin(6) + t242) * t192 + t211;
t162 = t39 * mrSges(5,1) + (Ifges(5,6) * t155 + t157 * t180) * t268 - (Ifges(4,6) * t155 + t157 * t181) * t204 / 0.2e1 - t45 * mrSges(4,3) - t28 * mrSges(5,2);
t33 = -pkin(3) * t193 - t44;
t161 = t33 * mrSges(5,2) - t44 * mrSges(4,3) - t39 * mrSges(5,3) + (t275 * t155 + t272 * t157) * t268;
t160 = t229 / 0.2e1 + t224 / 0.2e1 + t228 / 0.2e1 - t225 / 0.2e1;
t139 = mrSges(3,3) * t209 - t217;
t135 = mrSges(5,2) * t188;
t133 = mrSges(5,1) * t189;
t127 = -pkin(3) * t153 - t191;
t111 = t153 * t252 + t191;
t106 = (-mrSges(5,2) * t216 + mrSges(5,3) * t155) * t204;
t105 = -mrSges(5,1) * t193 + t135;
t104 = (mrSges(4,1) * t155 - mrSges(4,3) * t212) * t204;
t103 = (-mrSges(4,2) * t155 - mrSges(4,3) * t216) * t204;
t102 = t173 * t155;
t99 = t172 * t155;
t96 = -mrSges(5,3) * t188 + t133;
t86 = t178 * t209 + t148;
t83 = pkin(6) * t197 + t214;
t82 = t151 - t87;
t80 = t167 * t155;
t72 = t152 * t202 + t215;
t71 = qJD(1) * t186 - t214;
t68 = t144 + t84;
t67 = t172 * t207 - t195;
t66 = t171 * t209 - t148;
t58 = qJD(2) * t186 - t215;
t57 = Ifges(6,4) * t60;
t46 = t167 * t207 + t195;
t43 = t73 + t264;
t41 = qJD(5) * t101 + t164;
t40 = -t108 * t155 + t165;
t32 = qJD(2) * t170 + t264 + t98;
t31 = qJD(2) * t163 - t215;
t30 = (pkin(6) + t265) * t192 + t211;
t27 = -t119 * t252 + t166;
t19 = mrSges(6,2) * t193 + mrSges(6,3) * t26;
t18 = -mrSges(6,1) * t193 - mrSges(6,3) * t25;
t14 = Ifges(6,1) * t174 + Ifges(6,5) * t143 + t57;
t13 = Ifges(6,2) * t60 + Ifges(6,6) * t143 + t243;
t4 = -qJD(5) * t16 - t154 * t32 + t156 * t31;
t3 = qJD(5) * t15 + t154 * t31 + t156 * t32;
t5 = [(-t44 * mrSges(4,1) + t33 * mrSges(5,1) + t45 * mrSges(4,2) - t28 * mrSges(5,3) + ((-0.2e1 * t250 + (-0.3e1 / 0.2e1 * t229 - 0.3e1 / 0.2e1 * t224 - 0.3e1 / 0.2e1 * t228 + 0.3e1 / 0.2e1 * t225 + 0.3e1 / 0.2e1 * Ifges(3,4)) * t157 + (-0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(6,3) + m(4) * pkin(6) ^ 2 + (pkin(6) * mrSges(4,2) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t153) * t153 + (pkin(6) * mrSges(4,1) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t152 + t276 * t153) * t152) * t155) * qJD(1) + t266 + t261) * qJD(2) + t263) * t157 + (pkin(6) * t97 + t161 * t153 + t162 * t152 + (((-0.3e1 / 0.2e1 * Ifges(3,4) + t160) * t155 - 0.2e1 * t251 - Ifges(6,5) * t102 / 0.2e1 - Ifges(6,6) * t101 / 0.2e1) * qJD(1) - pkin(6) * t139 + t194 + t262) * qJD(2)) * t155 + (Ifges(6,5) * t41 + Ifges(6,6) * t40) * t247 + (Ifges(6,1) * t41 + Ifges(6,4) * t40) * t253 + (Ifges(6,4) * t41 + Ifges(6,2) * t40) * t255 + (Ifges(6,4) * t102 + Ifges(6,2) * t101) * t257 + (Ifges(6,1) * t102 + Ifges(6,4) * t101) * t258 - t30 * (-mrSges(6,1) * t101 + mrSges(6,2) * t102) + t88 * t103 + t87 * t104 + t82 * t105 + t81 * t106 + t72 * t91 + t58 * t92 + t99 * t96 + t80 * t7 + t43 * t89 + t73 * t90 + t67 * t69 + t40 * t13 / 0.2e1 + t27 * (-mrSges(6,1) * t40 + mrSges(6,2) * t41) + t41 * t14 / 0.2e1 + t46 * t17 + t3 * t36 + t4 * t37 + t15 * t18 + t16 * t19 + m(5) * (t28 * t81 + t33 * t82 + t39 * t99 + t42 * t67 + t43 * t49 + t47 * t58) + m(6) * (t1 * t16 + t15 * t2 + t27 * t46 + t3 * t9 - t30 * t80 + t4 * t8) + m(4) * (t44 * t87 + t45 * t88 + t64 * t72 + t65 * t73) + t102 * t259 + t101 * t260 + (t1 * t101 - t102 * t2 + t40 * t9 - t41 * t8) * mrSges(6,3); (-t108 / 0.2e1 - t94 / 0.2e1) * t14 + (-Ifges(6,5) * t108 - Ifges(6,6) * t109) * t247 + (-Ifges(6,1) * t108 - Ifges(6,4) * t109) * t253 + (-Ifges(6,4) * t108 - Ifges(6,2) * t109) * t255 + (-t109 / 0.2e1 - t93 / 0.2e1) * t13 + m(4) * ((-t152 * t64 + t153 * t65) * qJD(3) + (-t152 * t44 + t153 * t45) * qJ(3)) + (Ifges(6,4) * t122 - Ifges(6,2) * t173) * t257 + (Ifges(6,1) * t122 - Ifges(6,4) * t173) * t258 - t30 * (mrSges(6,1) * t173 + mrSges(6,2) * t122) + (-t1 * t173 - t2 * t122 + t221 * t9 - t222 * t8) * mrSges(6,3) - t173 * t260 + t273 * t37 + (t1 * t79 - t111 * t30 + t2 * t78 + t274 * t9 + t273 * t8 + (qJD(4) * t152 - t66) * t27) * m(6) + t274 * t36 + (((Ifges(6,5) * t122 - Ifges(6,6) * t173) * t278 + (t232 / 0.2e1 + t251) * qJD(1) + (t139 + t217) * pkin(6) + t194 + ((Ifges(4,6) - Ifges(5,6)) * t153 + t275 * t152) * t277 - t262) * t155 + ((t157 * t160 + t250) * qJD(1) - t146 / 0.2e1 + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t210 + ((Ifges(4,2) * t153 + t231) * t246 + (-Ifges(5,3) * t153 + t227) * t245 + Ifges(3,5) / 0.2e1 + (-m(4) * pkin(2) - mrSges(4,1) * t153 + mrSges(4,2) * t152 - mrSges(3,1)) * pkin(6) + (t279 * t152 - t226 + t230) * t244) * qJD(2) - t261) * t157) * qJD(1) + (-t42 * t86 - t47 * t71 - t49 * t68 + t127 * t39 + (qJ(3) * t28 + qJD(3) * t49) * t153 + (qJ(3) * t33 + qJD(3) * t47 - qJD(4) * t42) * t152) * m(5) + (Ifges(6,5) * t94 + Ifges(6,6) * t93) * t248 + (Ifges(6,1) * t94 + Ifges(6,4) * t93) * t254 + (Ifges(6,4) * t94 + Ifges(6,2) * t93) * t256 + t122 * t259 - m(4) * (t64 * t83 + t65 * t84) + t127 * t96 + t111 * t7 - t83 * t91 - t71 * t92 - pkin(2) * t97 + t78 * t18 + t79 * t19 - t86 * t69 - t68 * t89 - t84 * t90 - t66 * t17 + (-mrSges(6,1) * t221 + mrSges(6,2) * t222) * t27 + (t235 * qJD(3) + (t103 + t106) * qJ(3) - t162) * t153 + (t236 * qJD(4) - t234 * qJD(3) + (-t104 + t105) * qJ(3) + t161) * t152; t60 * t36 - t174 * t37 + t133 + t234 * t120 + t235 * t119 + (m(4) * pkin(6) - t233) * t192 - m(4) * (-t119 * t65 - t120 * t64) - t7 + t97 + (-t174 * t8 + t60 * t9 + t30) * m(6) + (t119 * t49 - t120 * t47 + t39) * m(5); t154 * t19 + t156 * t18 + t135 - t236 * t120 + t175 * qJD(5) + (-mrSges(5,1) * t208 + (t175 + t89) * t157) * qJD(1) + (t1 * t154 - t27 * t120 + t156 * t2 - t143 * (t154 * t8 - t156 * t9)) * m(6) + (t42 * t120 + t209 * t49 + t33) * m(5); -Ifges(6,3) * t193 - t27 * (mrSges(6,1) * t174 + mrSges(6,2) * t60) + (Ifges(6,1) * t60 - t243) * t254 + t13 * t253 + (Ifges(6,5) * t60 - Ifges(6,6) * t174) * t248 - t8 * t36 + t9 * t37 + (t174 * t9 + t60 * t8) * mrSges(6,3) + (-Ifges(6,2) * t174 + t14 + t57) * t256 + t263;];
tauc = t5(:);
