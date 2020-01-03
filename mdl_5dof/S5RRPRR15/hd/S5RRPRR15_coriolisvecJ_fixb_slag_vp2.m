% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR15_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR15_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:12
% EndTime: 2019-12-31 20:41:24
% DurationCPUTime: 5.29s
% Computational Cost: add. (4141->454), mult. (9564->633), div. (0->0), fcn. (5578->6), ass. (0->231)
t172 = cos(qJ(4));
t169 = sin(qJ(4));
t247 = Ifges(5,4) * t169;
t193 = Ifges(5,2) * t172 + t247;
t307 = -t193 / 0.2e1;
t293 = -qJD(1) / 0.2e1;
t170 = sin(qJ(2));
t173 = cos(qJ(2));
t184 = -pkin(8) * t169 * t170 + pkin(4) * t173;
t221 = qJD(4) * t169;
t174 = -pkin(2) - pkin(7);
t251 = pkin(8) - t174;
t226 = qJD(1) * t173;
t160 = pkin(6) * t226;
t129 = pkin(3) * t226 + t160;
t163 = t170 * qJD(1);
t159 = pkin(2) * t163;
t191 = pkin(7) * t170 - qJ(3) * t173;
t97 = qJD(1) * t191 + t159;
t57 = t172 * t129 - t169 * t97;
t306 = -qJD(1) * t184 + t251 * t221 - t57;
t135 = t251 * t172;
t209 = t172 * t163;
t58 = t169 * t129 + t172 * t97;
t305 = pkin(8) * t209 + qJD(4) * t135 + t58;
t167 = qJD(2) * qJ(3);
t108 = t167 + t129;
t205 = -qJ(3) * t170 - pkin(1);
t115 = t173 * t174 + t205;
t89 = t115 * qJD(1);
t158 = pkin(6) * t163;
t128 = -pkin(3) * t163 - t158;
t282 = qJD(3) - t128;
t93 = qJD(2) * t174 + t282;
t46 = -t169 * t89 + t172 * t93;
t47 = t169 * t93 + t172 * t89;
t189 = t46 * t169 - t47 * t172;
t198 = mrSges(5,1) * t172 - mrSges(5,2) * t169;
t255 = -t172 / 0.2e1;
t256 = -t169 / 0.2e1;
t120 = -qJD(2) * t169 - t172 * t226;
t153 = t163 + qJD(4);
t224 = qJD(2) * t172;
t121 = -t169 * t226 + t224;
t248 = Ifges(5,4) * t121;
t54 = t120 * Ifges(5,2) + t153 * Ifges(5,6) + t248;
t116 = Ifges(5,4) * t120;
t55 = t121 * Ifges(5,1) + t153 * Ifges(5,5) + t116;
t304 = t189 * mrSges(5,3) + t108 * t198 + t255 * t54 + t256 * t55;
t168 = sin(qJ(5));
t171 = cos(qJ(5));
t39 = pkin(8) * t120 + t47;
t234 = t171 * t39;
t38 = -pkin(8) * t121 + t46;
t34 = pkin(4) * t153 + t38;
t10 = t168 * t34 + t234;
t217 = qJD(1) * qJD(2);
t206 = t173 * t217;
t150 = Ifges(6,3) * t206;
t202 = t171 * t120 - t121 * t168;
t65 = t120 * t168 + t121 * t171;
t254 = Ifges(6,4) * t65;
t215 = qJD(4) + qJD(5);
t145 = t163 + t215;
t260 = -t145 / 0.2e1;
t273 = -t65 / 0.2e1;
t275 = -t202 / 0.2e1;
t59 = Ifges(6,4) * t202;
t31 = Ifges(6,1) * t65 + Ifges(6,5) * t145 + t59;
t72 = -pkin(4) * t120 + t108;
t236 = t168 * t39;
t9 = t171 * t34 - t236;
t303 = t150 + (Ifges(6,5) * t202 - Ifges(6,6) * t65) * t260 + (t10 * t65 + t202 * t9) * mrSges(6,3) + (-Ifges(6,2) * t65 + t31 + t59) * t275 - t72 * (mrSges(6,1) * t65 + mrSges(6,2) * t202) + (Ifges(6,1) * t202 - t254) * t273;
t269 = pkin(3) + pkin(6);
t302 = -mrSges(3,1) + mrSges(4,2);
t136 = -pkin(2) * t173 + t205;
t109 = t136 * qJD(1);
t139 = -t160 - t167;
t241 = Ifges(5,6) * t172;
t245 = Ifges(5,5) * t169;
t192 = t241 + t245;
t246 = Ifges(5,4) * t172;
t195 = Ifges(5,1) * t169 + t246;
t257 = t153 / 0.2e1;
t263 = t121 / 0.2e1;
t291 = qJD(2) / 0.2e1;
t292 = -qJD(2) / 0.2e1;
t301 = -t139 * mrSges(4,1) + t109 * mrSges(4,2) - Ifges(4,5) * t291 - Ifges(3,6) * t292 + t120 * t307 - t192 * t257 - t195 * t263 + t304 + ((-Ifges(3,2) - Ifges(4,3)) * t173 + (-Ifges(3,4) - Ifges(4,6)) * t170) * t293;
t223 = qJD(2) * t173;
t131 = t269 * t223;
t114 = qJD(1) * t131;
t207 = t170 * t217;
t152 = pkin(2) * t207;
t222 = qJD(3) * t170;
t177 = qJD(2) * t191 - t222;
t76 = qJD(1) * t177 + t152;
t20 = -qJD(4) * t47 + t172 * t114 - t169 * t76;
t216 = qJD(2) * qJD(4);
t219 = qJD(4) * t173;
t225 = qJD(2) * t170;
t80 = -t169 * t216 + (t169 * t225 - t172 * t219) * qJD(1);
t11 = pkin(4) * t206 - pkin(8) * t80 + t20;
t220 = qJD(4) * t172;
t19 = t169 * t114 + t172 * t76 + t93 * t220 - t221 * t89;
t208 = t169 * t219;
t178 = t170 * t224 + t208;
t81 = qJD(1) * t178 - t172 * t216;
t14 = pkin(8) * t81 + t19;
t2 = qJD(5) * t9 + t11 * t168 + t14 * t171;
t25 = qJD(5) * t202 + t168 * t81 + t171 * t80;
t26 = -qJD(5) * t65 - t168 * t80 + t171 * t81;
t284 = qJD(5) * t10;
t3 = t11 * t171 - t14 * t168 - t284;
t300 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t25 + Ifges(6,6) * t26;
t299 = qJD(3) + t158;
t228 = t171 * t172;
t229 = t168 * t169;
t185 = -t228 + t229;
t186 = t168 * t172 + t171 * t169;
t70 = t215 * t186;
t180 = t186 * t170;
t91 = qJD(1) * t180;
t249 = -t70 - t91;
t218 = qJD(5) * t168;
t69 = -t168 * t221 - t169 * t218 + t215 * t228;
t90 = -t163 * t229 + t171 * t209;
t250 = t69 + t90;
t298 = -t10 * t250 + t185 * t3 - t186 * t2 - t249 * t9;
t30 = Ifges(6,2) * t202 + Ifges(6,6) * t145 + t254;
t296 = t30 / 0.2e1;
t134 = t251 * t169;
t73 = t134 * t168 - t135 * t171;
t290 = qJD(5) * t73 + t306 * t168 - t305 * t171;
t74 = -t134 * t171 - t135 * t168;
t289 = -qJD(5) * t74 + t305 * t168 + t306 * t171;
t210 = -pkin(4) * t172 - pkin(3);
t285 = pkin(4) * t220 - t163 * t210 + t299;
t143 = t269 * t170;
t124 = t169 * t143;
t68 = t172 * t115 + t124;
t283 = t169 * t19 + t172 * t20;
t279 = t173 * t215;
t278 = t20 * mrSges(5,1) - t19 * mrSges(5,2) + Ifges(5,5) * t80 + Ifges(5,6) * t81 + t300;
t277 = t25 / 0.2e1;
t276 = t26 / 0.2e1;
t274 = t202 / 0.2e1;
t272 = t65 / 0.2e1;
t98 = t185 * t173;
t271 = t98 / 0.2e1;
t99 = t186 * t173;
t270 = -t99 / 0.2e1;
t268 = pkin(1) * mrSges(3,1);
t267 = pkin(1) * mrSges(3,2);
t265 = -t120 / 0.2e1;
t264 = -t121 / 0.2e1;
t262 = -t186 / 0.2e1;
t261 = -t185 / 0.2e1;
t259 = t145 / 0.2e1;
t258 = -t153 / 0.2e1;
t244 = Ifges(5,5) * t172;
t243 = Ifges(4,6) * t173;
t242 = Ifges(5,6) * t169;
t141 = -mrSges(4,1) * t226 - qJD(2) * mrSges(4,3);
t71 = -mrSges(5,1) * t120 + mrSges(5,2) * t121;
t232 = -t141 + t71;
t231 = qJD(2) * mrSges(3,2);
t227 = t172 * t173;
t144 = t269 * t173;
t214 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t213 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t212 = -0.3e1 / 0.2e1 * Ifges(4,6) - 0.3e1 / 0.2e1 * Ifges(3,4);
t211 = m(4) * pkin(6) + mrSges(4,1);
t204 = pkin(8) * t173 - t115;
t162 = pkin(2) * t225;
t85 = t162 + t177;
t203 = t172 * t131 - t169 * t85;
t130 = t269 * t225;
t133 = -qJD(2) * pkin(2) + t299;
t201 = m(4) * t133 + (mrSges(4,1) + mrSges(3,3)) * t163 + t302 * qJD(2);
t200 = m(4) * t139 - mrSges(3,3) * t226 + t141 + t231;
t197 = mrSges(5,1) * t169 + mrSges(5,2) * t172;
t196 = Ifges(5,1) * t172 - t247;
t194 = -Ifges(5,2) * t169 + t246;
t125 = t172 * t143;
t52 = pkin(4) * t170 + t169 * t204 + t125;
t56 = -pkin(8) * t227 + t68;
t27 = -t168 * t56 + t171 * t52;
t28 = t168 * t52 + t171 * t56;
t60 = mrSges(5,1) * t206 - mrSges(5,3) * t80;
t61 = -mrSges(5,2) * t206 + mrSges(5,3) * t81;
t188 = t169 * t61 + t172 * t60;
t83 = -mrSges(5,2) * t153 + mrSges(5,3) * t120;
t84 = mrSges(5,1) * t153 - mrSges(5,3) * t121;
t187 = -t169 * t84 + t172 * t83;
t181 = -qJ(3) * t223 - t222;
t32 = -t115 * t221 + t169 * t131 + t143 * t220 + t172 * t85;
t166 = qJD(2) * qJD(3);
t94 = -qJD(1) * t130 + t166;
t157 = Ifges(3,4) * t226;
t176 = t133 * mrSges(4,1) + t46 * mrSges(5,1) + t9 * mrSges(6,1) + Ifges(3,1) * t163 / 0.2e1 + Ifges(3,5) * t291 + t157 / 0.2e1 + Ifges(4,4) * t292 + (-Ifges(4,2) * t170 - t243) * t293 + t145 * Ifges(6,3) + t65 * Ifges(6,5) + t202 * Ifges(6,6) + t153 * Ifges(5,3) + t121 * Ifges(5,5) + t120 * Ifges(5,6) - t10 * mrSges(6,2) - t109 * mrSges(4,3) - t47 * mrSges(5,2);
t154 = pkin(4) * t169 + qJ(3);
t151 = Ifges(5,3) * t206;
t132 = pkin(6) * t207 - t166;
t126 = (mrSges(4,2) * t173 - mrSges(4,3) * t170) * qJD(1);
t103 = pkin(4) * t227 + t144;
t101 = t162 + t181;
t87 = qJD(1) * t181 + t152;
t75 = -pkin(4) * t208 + (-pkin(6) + t210) * t225;
t67 = -t115 * t169 + t125;
t51 = -pkin(4) * t81 + t94;
t50 = mrSges(6,1) * t145 - mrSges(6,3) * t65;
t49 = -mrSges(6,2) * t145 + mrSges(6,3) * t202;
t44 = -mrSges(5,1) * t81 + mrSges(5,2) * t80;
t43 = t80 * Ifges(5,1) + t81 * Ifges(5,4) + Ifges(5,5) * t206;
t42 = t80 * Ifges(5,4) + t81 * Ifges(5,2) + Ifges(5,6) * t206;
t41 = -t185 * t225 + t186 * t279;
t40 = qJD(2) * t180 + t185 * t279;
t35 = -mrSges(6,1) * t202 + mrSges(6,2) * t65;
t33 = -qJD(4) * t68 + t203;
t24 = pkin(8) * t178 + t32;
t21 = t184 * qJD(2) + (t172 * t204 - t124) * qJD(4) + t203;
t18 = -mrSges(6,2) * t206 + mrSges(6,3) * t26;
t17 = mrSges(6,1) * t206 - mrSges(6,3) * t25;
t13 = t171 * t38 - t236;
t12 = -t168 * t38 - t234;
t8 = -mrSges(6,1) * t26 + mrSges(6,2) * t25;
t7 = t25 * Ifges(6,1) + t26 * Ifges(6,4) + Ifges(6,5) * t206;
t6 = t25 * Ifges(6,4) + t26 * Ifges(6,2) + Ifges(6,6) * t206;
t5 = -qJD(5) * t28 - t168 * t24 + t171 * t21;
t4 = qJD(5) * t27 + t168 * t21 + t171 * t24;
t1 = [(t10 * t41 + t2 * t98 + t3 * t99 - t40 * t9) * mrSges(6,3) + t51 * (-mrSges(6,1) * t98 - mrSges(6,2) * t99) + (-Ifges(6,4) * t99 + Ifges(6,2) * t98) * t276 + (-Ifges(6,1) * t99 + Ifges(6,4) * t98) * t277 + m(5) * (-t108 * t130 + t144 * t94 + t19 * t68 + t20 * t67 + t32 * t47 + t33 * t46) + m(4) * (t101 * t109 + t136 * t87) + (-t87 * mrSges(4,3) + t150 / 0.2e1 + t151 / 0.2e1 + ((-t136 * mrSges(4,2) + t170 * t212 - 0.2e1 * t268) * qJD(1) + t213 * qJD(2) + t200 * pkin(6) - t301) * qJD(2) + t278) * t170 + (t81 * t307 - t80 * t195 / 0.2e1 + t43 * t256 + t42 * t255 + t87 * mrSges(4,2) + t94 * t198 - t211 * t132 + (t169 * t20 - t172 * t19) * mrSges(5,3) + (t196 * t264 + t194 * t265 - t108 * t197 + (t242 - t244) * t257 + t169 * t54 / 0.2e1 + t55 * t255 + (t169 * t47 + t172 * t46) * mrSges(5,3)) * qJD(4) + (t214 * qJD(2) + t201 * pkin(6) + ((-t245 / 0.2e1 - t241 / 0.2e1 - t212) * t173 + Ifges(6,5) * t270 + Ifges(6,6) * t271 - 0.2e1 * t267 - t136 * mrSges(4,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + t211 * pkin(6)) * t170) * qJD(1) + t176) * qJD(2)) * t173 + t144 * t44 - t130 * t71 + t101 * t126 + t103 * t8 + t32 * t83 + t33 * t84 + t72 * (-mrSges(6,1) * t41 + mrSges(6,2) * t40) + t75 * t35 + t67 * t60 + t68 * t61 + t4 * t49 + t5 * t50 + t40 * t31 / 0.2e1 + t27 * t17 + t28 * t18 + m(6) * (t10 * t4 + t103 * t51 + t2 * t28 + t27 * t3 + t5 * t9 + t72 * t75) + (Ifges(6,4) * t40 + Ifges(6,2) * t41) * t274 + t7 * t270 + t6 * t271 + (Ifges(6,1) * t40 + Ifges(6,4) * t41) * t272 + (Ifges(6,5) * t40 + Ifges(6,6) * t41) * t259 + t41 * t296; (-m(4) * t109 - t126) * (-qJ(3) * t226 + t159) + m(4) * (-qJ(3) * t132 - qJD(3) * t139) + (((t268 + (Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1) * t170) * qJD(1) + (-t200 + t231) * pkin(6) + (-qJ(3) * mrSges(4,1) + t213) * qJD(2) + t301) * t170 + (((-m(4) * pkin(2) + t302) * qJD(2) - t201) * pkin(6) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t163 - t176 - t157 / 0.2e1 + (-t243 / 0.2e1 + t267) * qJD(1) + (-pkin(2) * mrSges(4,1) + Ifges(6,5) * t261 + Ifges(6,6) * t262 + t244 / 0.2e1 - t242 / 0.2e1 + t214) * qJD(2)) * t173) * qJD(1) + t298 * mrSges(6,3) + (-t90 / 0.2e1 - t69 / 0.2e1) * t30 + t285 * t35 + (t94 * qJ(3) + t108 * t282 - t46 * t57 - t47 * t58) * m(5) + ((-m(5) * t189 + t187) * qJD(4) + m(5) * t283 + t188) * t174 - t283 * mrSges(5,3) + (t192 * t258 + t193 * t265 + t195 * t264 + t304) * qJD(4) - t128 * t71 - t132 * mrSges(4,3) - t58 * t83 - t57 * t84 + t73 * t17 + t74 * t18 + qJ(3) * t44 + t289 * t50 + t290 * t49 + (t10 * t290 + t154 * t51 + t2 * t74 + t285 * t72 + t289 * t9 + t3 * t73) * m(6) + (Ifges(6,4) * t91 + Ifges(6,2) * t90) * t275 + t7 * t261 + t6 * t262 + (Ifges(6,1) * t91 + Ifges(6,4) * t90) * t273 + t42 * t256 + (Ifges(6,5) * t91 + Ifges(6,6) * t90) * t260 + t51 * (mrSges(6,1) * t186 - mrSges(6,2) * t185) + (-Ifges(6,4) * t185 - Ifges(6,2) * t186) * t276 + (-Ifges(6,1) * t185 - Ifges(6,4) * t186) * t277 + (-t91 / 0.2e1 - t70 / 0.2e1) * t31 + (-Ifges(6,4) * t70 - Ifges(6,2) * t69) * t274 + (-Ifges(6,1) * t70 - Ifges(6,4) * t69) * t272 + (-Ifges(6,5) * t70 - Ifges(6,6) * t69) * t259 + (mrSges(6,1) * t250 + mrSges(6,2) * t249) * t72 + t154 * t8 + t172 * t43 / 0.2e1 + t81 * t194 / 0.2e1 + t80 * t196 / 0.2e1 + t94 * t197 + t232 * qJD(3); t186 * t18 - t185 * t17 + t249 * t50 + t250 * t49 + t187 * qJD(4) + (-t35 - t232) * qJD(2) + (t211 * t223 + (t126 + t187) * t170) * qJD(1) - m(4) * (-qJD(2) * t139 - t109 * t163) + t188 + (-qJD(2) * t72 - t298) * m(6) + (-qJD(2) * t108 - t153 * t189 + t283) * m(5); (-Ifges(5,2) * t121 + t116 + t55) * t265 + t278 + (-t121 * t35 + t168 * t18 + t171 * t17 + (-t168 * t50 + t171 * t49) * qJD(5) + (-t121 * t72 + t168 * t2 - t9 * t218 + (t3 + t284) * t171) * m(6)) * pkin(4) - m(6) * (t10 * t13 + t12 * t9) - t108 * (mrSges(5,1) * t121 + mrSges(5,2) * t120) - t46 * t83 + t47 * t84 - t13 * t49 - t12 * t50 + t151 + t54 * t263 + (Ifges(5,1) * t120 - t248) * t264 + (Ifges(5,5) * t120 - Ifges(5,6) * t121) * t258 + t65 * t296 + (t120 * t46 + t121 * t47) * mrSges(5,3) + t303; t10 * t50 + t30 * t272 - t9 * t49 + t300 + t303;];
tauc = t1(:);
