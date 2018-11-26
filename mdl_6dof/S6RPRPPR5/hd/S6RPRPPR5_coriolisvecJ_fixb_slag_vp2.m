% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:54:13
% EndTime: 2018-11-23 15:54:20
% DurationCPUTime: 6.60s
% Computational Cost: add. (6279->521), mult. (16434->700), div. (0->0), fcn. (12107->8), ass. (0->227)
t189 = cos(pkin(9));
t268 = cos(qJ(3));
t228 = t268 * t189;
t223 = qJD(1) * t228;
t187 = sin(pkin(9));
t192 = sin(qJ(3));
t241 = t187 * t192;
t227 = qJD(3) * t241;
t143 = qJD(1) * t227 - qJD(3) * t223;
t163 = t187 * t268 + t192 * t189;
t157 = t163 * qJD(3);
t144 = qJD(1) * t157;
t186 = sin(pkin(10));
t152 = qJD(1) * t241 - t223;
t153 = t163 * qJD(1);
t209 = qJ(4) * t143 - qJD(4) * t153;
t260 = pkin(3) + qJ(5);
t47 = qJD(5) * t152 + t144 * t260 + t209;
t188 = cos(pkin(10));
t199 = t163 * qJD(2);
t197 = qJD(1) * t199;
t259 = pkin(7) + qJ(2);
t169 = t259 * t187;
t164 = qJD(1) * t169;
t170 = t259 * t189;
t165 = qJD(1) * t170;
t291 = -t192 * t164 + t268 * t165;
t59 = -t143 * pkin(4) + t197 + (-qJD(5) + t291) * qJD(3);
t55 = t188 * t59;
t13 = -pkin(5) * t143 + t55 + (-pkin(8) * t144 - t47) * t186;
t18 = t186 * t59 + t188 * t47;
t245 = t144 * t188;
t15 = pkin(8) * t245 + t18;
t191 = sin(qJ(6));
t193 = cos(qJ(6));
t132 = -t188 * qJD(3) - t186 * t152;
t230 = -pkin(2) * t189 - pkin(1);
t168 = qJD(1) * t230 + qJD(2);
t198 = -qJ(4) * t153 + t168;
t69 = t152 * t260 + t198;
t125 = t268 * t164 + t165 * t192;
t205 = pkin(4) * t153 + t125;
t300 = qJD(4) + t205;
t77 = -qJD(3) * t260 + t300;
t32 = -t186 * t69 + t188 * t77;
t20 = pkin(5) * t153 + pkin(8) * t132 + t32;
t133 = -qJD(3) * t186 + t152 * t188;
t33 = t186 * t77 + t188 * t69;
t21 = pkin(8) * t133 + t33;
t5 = -t191 * t21 + t193 * t20;
t1 = qJD(6) * t5 + t13 * t191 + t15 * t193;
t6 = t191 * t20 + t193 * t21;
t2 = -qJD(6) * t6 + t13 * t193 - t15 * t191;
t211 = t193 * t186 + t191 * t188;
t155 = t211 * qJD(6);
t200 = t211 * t153;
t235 = -t155 - t200;
t289 = -t186 * t191 + t188 * t193;
t154 = t289 * qJD(6);
t292 = t289 * t153;
t290 = t154 + t292;
t303 = -t1 * t211 - t2 * t289 - t235 * t5 - t290 * t6;
t262 = mrSges(5,2) - mrSges(4,1);
t302 = mrSges(4,3) + mrSges(5,1);
t225 = t132 * t191 + t193 * t133;
t41 = qJD(6) * t225 + t211 * t144;
t283 = t41 / 0.2e1;
t78 = t132 * t193 - t133 * t191;
t42 = qJD(6) * t78 + t144 * t289;
t282 = t42 / 0.2e1;
t299 = -t143 / 0.2e1;
t258 = -pkin(8) - t260;
t166 = t258 * t186;
t167 = t258 * t188;
t128 = t166 * t193 + t167 * t191;
t95 = -pkin(4) * t152 + t291;
t88 = t188 * t95;
t248 = qJ(4) * t152;
t89 = t153 * t260 + t248;
t24 = -pkin(5) * t152 + t88 + (-pkin(8) * t153 - t89) * t186;
t266 = pkin(8) * t188;
t37 = t186 * t95 + t188 * t89;
t31 = t153 * t266 + t37;
t298 = -qJD(5) * t289 - qJD(6) * t128 + t191 * t31 - t193 * t24;
t127 = -t166 * t191 + t167 * t193;
t297 = -qJD(5) * t211 + qJD(6) * t127 - t191 * t24 - t193 * t31;
t295 = -Ifges(4,1) - Ifges(6,3);
t294 = Ifges(5,4) - Ifges(4,5);
t293 = Ifges(5,5) - Ifges(4,6);
t173 = qJD(2) * t228;
t226 = qJD(3) * t268;
t232 = qJD(1) * qJD(2);
t90 = -t192 * (qJD(3) * t165 + t187 * t232) + qJD(1) * t173 - t164 * t226;
t195 = qJD(3) * qJD(4) + t90;
t288 = -qJD(4) - t125;
t287 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t41 + Ifges(7,6) * t42;
t286 = (m(3) * qJ(2) + mrSges(3,3)) * (t187 ^ 2 + t189 ^ 2);
t102 = t192 * (qJD(2) * t187 + qJD(3) * t170) + t169 * t226 - t173;
t285 = Ifges(7,4) * t283 + Ifges(7,2) * t282 + Ifges(7,6) * t299;
t284 = Ifges(7,1) * t283 + Ifges(7,4) * t282 + Ifges(7,5) * t299;
t281 = -t225 / 0.2e1;
t280 = t225 / 0.2e1;
t279 = t78 / 0.2e1;
t278 = -t78 / 0.2e1;
t277 = t5 * mrSges(7,1);
t276 = t6 * mrSges(7,2);
t148 = qJD(6) + t153;
t275 = -t148 / 0.2e1;
t274 = t148 / 0.2e1;
t271 = t186 / 0.2e1;
t270 = -t188 / 0.2e1;
t269 = t188 / 0.2e1;
t267 = Ifges(7,4) * t78;
t265 = t225 * Ifges(7,6);
t264 = t78 * Ifges(7,5);
t261 = Ifges(4,4) + Ifges(5,6);
t161 = -t228 + t241;
t156 = -t189 * t226 + t227;
t208 = qJ(4) * t156 - qJD(4) * t163;
t51 = qJD(5) * t161 + t157 * t260 + t208;
t131 = -t192 * t169 + t170 * t268;
t103 = qJD(3) * t131 + t199;
t71 = -t156 * pkin(4) + t103;
t23 = t186 * t71 + t188 * t51;
t257 = Ifges(6,1) * t186;
t256 = Ifges(6,4) * t186;
t255 = Ifges(6,4) * t188;
t254 = Ifges(6,5) * t186;
t253 = Ifges(6,6) * t188;
t130 = t268 * t169 + t170 * t192;
t91 = qJD(3) * t291 + t197;
t252 = t130 * t91;
t251 = t132 * Ifges(6,5);
t250 = t133 * Ifges(6,6);
t249 = t148 * Ifges(7,3);
t112 = pkin(4) * t163 + t130;
t206 = -qJ(4) * t163 + t230;
t93 = t161 * t260 + t206;
t45 = t186 * t112 + t188 * t93;
t246 = t144 * t186;
t244 = t157 * t186;
t243 = t157 * t188;
t135 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t152;
t137 = t152 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t239 = t135 - t137;
t238 = -qJD(3) * t262 - t302 * t153;
t34 = -mrSges(7,1) * t225 - mrSges(7,2) * t78;
t85 = -mrSges(6,1) * t133 - mrSges(6,2) * t132;
t231 = t137 - t34 - t85;
t229 = -pkin(5) * t188 - pkin(4);
t14 = -t42 * mrSges(7,1) + t41 * mrSges(7,2);
t98 = -mrSges(6,1) * t245 + mrSges(6,2) * t246;
t220 = -mrSges(6,1) * t188 + mrSges(6,2) * t186;
t219 = t255 + t257;
t218 = Ifges(6,2) * t188 + t256;
t217 = t253 + t254;
t17 = -t186 * t47 + t55;
t216 = t17 * t188 + t18 * t186;
t215 = -t17 * t186 + t18 * t188;
t214 = t186 * t33 + t188 * t32;
t213 = t186 * t32 - t188 * t33;
t101 = t188 * t112;
t28 = pkin(5) * t163 + t101 + (-pkin(8) * t161 - t93) * t186;
t35 = t161 * t266 + t45;
t9 = -t191 * t35 + t193 * t28;
t10 = t191 * t28 + t193 * t35;
t212 = -t130 * t143 - t131 * t144;
t204 = -t254 / 0.2e1 - t253 / 0.2e1;
t185 = qJD(3) * qJ(4);
t84 = qJD(5) + t185 + t95;
t203 = t84 * t220;
t114 = t289 * t161;
t177 = pkin(5) * t186 + qJ(4);
t146 = Ifges(4,4) * t152;
t145 = Ifges(5,6) * t152;
t141 = Ifges(7,3) * t143;
t140 = t143 * mrSges(4,2);
t139 = t143 * mrSges(5,3);
t124 = pkin(3) * t161 + t206;
t123 = -mrSges(5,2) * t152 - mrSges(5,3) * t153;
t122 = pkin(3) * t153 + t248;
t121 = -t185 - t291;
t120 = -qJD(3) * pkin(3) - t288;
t119 = t153 * Ifges(4,1) + Ifges(4,5) * qJD(3) - t146;
t118 = Ifges(4,4) * t153 - t152 * Ifges(4,2) + Ifges(4,6) * qJD(3);
t117 = Ifges(5,4) * qJD(3) - t153 * Ifges(5,2) + t145;
t116 = Ifges(5,5) * qJD(3) - Ifges(5,6) * t153 + t152 * Ifges(5,3);
t115 = t211 * t161;
t113 = -t161 * pkin(4) + t131;
t111 = mrSges(6,2) * t143 + mrSges(6,3) * t245;
t110 = -mrSges(6,1) * t143 - mrSges(6,3) * t246;
t105 = mrSges(6,1) * t153 + mrSges(6,3) * t132;
t104 = -mrSges(6,2) * t153 + mrSges(6,3) * t133;
t99 = pkin(3) * t152 + t198;
t92 = pkin(3) * t157 + t208;
t76 = Ifges(7,4) * t225;
t75 = pkin(3) * t144 + t209;
t74 = t161 * t229 + t131;
t70 = -pkin(4) * t157 - t102;
t68 = t188 * t71;
t65 = -t143 * Ifges(6,5) + t144 * t219;
t64 = -t143 * Ifges(6,6) + t144 * t218;
t63 = t153 * t229 - t125;
t62 = -t132 * Ifges(6,1) + t133 * Ifges(6,4) + Ifges(6,5) * t153;
t61 = -t132 * Ifges(6,4) + t133 * Ifges(6,2) + Ifges(6,6) * t153;
t60 = t153 * Ifges(6,3) + t250 - t251;
t58 = -pkin(4) * t144 + t195;
t57 = mrSges(7,1) * t148 + mrSges(7,3) * t78;
t56 = -mrSges(7,2) * t148 + mrSges(7,3) * t225;
t53 = -pkin(5) * t133 + t84;
t52 = t157 * t229 - t102;
t50 = -t155 * t161 + t157 * t289;
t49 = qJD(6) * t114 + t157 * t211;
t46 = t144 * t229 + t195;
t44 = -t186 * t93 + t101;
t36 = -t186 * t89 + t88;
t30 = mrSges(7,2) * t143 + mrSges(7,3) * t42;
t29 = -mrSges(7,1) * t143 - mrSges(7,3) * t41;
t27 = -Ifges(7,1) * t78 + Ifges(7,5) * t148 + t76;
t26 = Ifges(7,2) * t225 + Ifges(7,6) * t148 - t267;
t25 = t249 - t264 + t265;
t22 = -t186 * t51 + t68;
t19 = pkin(8) * t243 + t23;
t16 = -pkin(5) * t156 + t68 + (-pkin(8) * t157 - t51) * t186;
t4 = -qJD(6) * t10 + t16 * t193 - t19 * t191;
t3 = qJD(6) * t9 + t16 * t191 + t19 * t193;
t7 = [m(5) * (t102 * t121 + t103 * t120 + t124 * t75 + t131 * t195 + t92 * t99 + t252) + (t58 * t220 - t75 * mrSges(5,2) - t195 * mrSges(5,1) - t90 * mrSges(4,3) + t65 * t271 + t64 * t269 + (t204 + t261) * t143 + t215 * mrSges(6,3) + (Ifges(5,3) + Ifges(4,2) + Ifges(6,2) * t188 ^ 2 / 0.2e1 + (t255 + t257 / 0.2e1) * t186) * t144) * t161 + (Ifges(7,4) * t49 + Ifges(7,2) * t50 - Ifges(7,6) * t156) * t280 + (Ifges(7,4) * t115 + Ifges(7,2) * t114) * t282 + (Ifges(7,1) * t115 + Ifges(7,4) * t114) * t283 + t115 * t284 + t114 * t285 + (Ifges(7,5) * t49 + Ifges(7,6) * t50 - Ifges(7,3) * t156) * t274 + t156 * t276 + (Ifges(7,1) * t49 + Ifges(7,4) * t50 - Ifges(7,5) * t156) * t278 + t157 * t203 + t133 * (-Ifges(6,6) * t156 + t157 * t218) / 0.2e1 - t132 * (-Ifges(6,5) * t156 + t157 * t219) / 0.2e1 + (-t120 * t156 + t121 * t157 + t212) * mrSges(5,1) - (t119 + t60 + t25) * t156 / 0.2e1 + m(7) * (t1 * t10 + t2 * t9 + t3 * t6 + t4 * t5 + t46 * t74 + t52 * t53) + m(6) * (t113 * t58 + t17 * t44 + t18 * t45 + t22 * t32 + t23 * t33 + t70 * t84) + (t1 * t114 - t115 * t2 - t49 * t5 + t50 * t6) * mrSges(7,3) + (-t18 * mrSges(6,2) + t17 * mrSges(6,1) - t75 * mrSges(5,3) - t141 / 0.2e1 + t302 * t91 + (t217 - t261) * t144 + (-Ifges(5,2) - Ifges(7,3) / 0.2e1 + t295) * t143 + t287) * t163 + t61 * t243 / 0.2e1 + t33 * (mrSges(6,2) * t156 + mrSges(6,3) * t243) + t62 * t244 / 0.2e1 + t32 * (-mrSges(6,1) * t156 - mrSges(6,3) * t244) + (Ifges(7,5) * t115 + Ifges(7,6) * t114) * t299 - t238 * t103 - t239 * t102 + t230 * (t144 * mrSges(4,1) - t140) + 0.2e1 * t286 * t232 + t168 * (mrSges(4,1) * t157 - mrSges(4,2) * t156) + t99 * (-mrSges(5,2) * t157 + mrSges(5,3) * t156) + t157 * t116 / 0.2e1 - t152 * (-Ifges(4,4) * t156 - Ifges(4,2) * t157) / 0.2e1 - t157 * t118 / 0.2e1 - t153 * (Ifges(5,2) * t156 + Ifges(5,6) * t157) / 0.2e1 + t152 * (Ifges(5,6) * t156 + Ifges(5,3) * t157) / 0.2e1 + t156 * t117 / 0.2e1 + t124 * (-t144 * mrSges(5,2) + t139) + t92 * t123 + t45 * t111 + t113 * t98 + t46 * (-mrSges(7,1) * t114 + mrSges(7,2) * t115) + (t156 * t294 + t157 * t293) * qJD(3) / 0.2e1 + ((-Ifges(4,4) + t217) * t157 + t295 * t156) * t153 / 0.2e1 - t156 * t277 + (-t125 * t156 - t157 * t291 + t212) * mrSges(4,3) + m(4) * (-t102 * t291 + t103 * t125 + t131 * t90 + t252) + t9 * t29 + t10 * t30 + t49 * t27 / 0.2e1 + t50 * t26 / 0.2e1 + t52 * t34 + t53 * (-mrSges(7,1) * t50 + mrSges(7,2) * t49) + t3 * t56 + t4 * t57 + t74 * t14 + t70 * t85 + t23 * t104 + t22 * t105 + t44 * t110; -t186 * t110 + t188 * t111 - t211 * t29 + t289 * t30 + t139 - t140 - t290 * t57 + t235 * t56 - t262 * t144 + (-t104 * t186 - t105 * t188 + t238) * t153 + (t135 - t231) * t152 - m(4) * (t125 * t153 - t152 * t291) - t286 * qJD(1) ^ 2 + (t1 * t289 + t152 * t53 - t2 * t211 + t235 * t6 - t290 * t5) * m(7) + (t152 * t84 - t153 * t214 + t215) * m(6) + (-t120 * t153 - t121 * t152 + t75) * m(5); (-pkin(3) * t91 + qJ(4) * t195 - t120 * t291 + t121 * t288 - t122 * t99) * m(5) + t195 * mrSges(5,3) + (mrSges(7,1) * t290 + mrSges(7,2) * t235) * t53 + (-Ifges(7,4) * t155 - Ifges(7,2) * t154) * t280 + (-Ifges(7,5) * t155 - Ifges(7,6) * t154) * t274 + (-Ifges(7,1) * t155 - Ifges(7,4) * t154) * t278 + (-t155 / 0.2e1 - t200 / 0.2e1) * t27 + (qJ(4) * t58 - qJD(5) * t214 - t260 * t216 + t300 * t84 - t32 * t36 - t33 * t37) * m(6) + (Ifges(7,1) * t289 - Ifges(7,4) * t211) * t283 + t46 * (mrSges(7,1) * t211 + mrSges(7,2) * t289) + (Ifges(6,5) * t270 + Ifges(6,6) * t271 - Ifges(7,5) * t289 / 0.2e1 + Ifges(7,6) * t211 / 0.2e1 + pkin(3) * mrSges(5,1) + t294) * t143 + (Ifges(7,4) * t289 - Ifges(7,2) * t211) * t282 + t289 * t284 + (t58 * mrSges(6,2) + t65 / 0.2e1 - t17 * mrSges(6,3) - qJD(5) * t105 - t260 * t110) * t188 + (t58 * mrSges(6,1) - t64 / 0.2e1 - t18 * mrSges(6,3) - qJD(5) * t104 - t260 * t111) * t186 + t239 * t125 - t231 * qJD(4) + (-t154 / 0.2e1 - t292 / 0.2e1) * t26 + (Ifges(7,1) * t200 + Ifges(7,4) * t292) * t279 + (Ifges(7,4) * t200 + Ifges(7,2) * t292) * t281 + (Ifges(7,5) * t200 + Ifges(7,6) * t292) * t275 + t177 * t14 + t127 * t29 + t128 * t30 - t122 * t123 + t262 * t91 + (-qJ(4) * mrSges(5,1) + (Ifges(6,1) * t188 - t256) * t271 + (-Ifges(6,2) * t186 + t255) * t269 + t293) * t144 + t205 * t85 + t238 * t291 + (t132 * t219 / 0.2e1 - t133 * t218 / 0.2e1 - t203 - t168 * mrSges(4,1) + t118 / 0.2e1 + t99 * mrSges(5,2) - t116 / 0.2e1 - t121 * mrSges(5,1) + t291 * mrSges(4,3) - t186 * t62 / 0.2e1 + t61 * t270 + (Ifges(4,4) / 0.2e1 + Ifges(5,6) / 0.2e1 + t204) * t153 + (Ifges(4,6) / 0.2e1 - Ifges(5,5) / 0.2e1) * qJD(3) + t213 * mrSges(6,3) + (Ifges(6,3) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(5,3) / 0.2e1) * t152) * t153 - t211 * t285 + t297 * t56 + t298 * t57 + (t1 * t128 + t127 * t2 + t177 * t46 + t297 * t6 + (qJD(4) - t63) * t53 + t298 * t5) * m(7) + (-t33 * mrSges(6,2) + t32 * mrSges(6,1) - t99 * mrSges(5,3) + t168 * mrSges(4,2) - t146 / 0.2e1 - t145 / 0.2e1 + t277 + t265 / 0.2e1 - t264 / 0.2e1 + t249 / 0.2e1 + t119 / 0.2e1 + t25 / 0.2e1 + t60 / 0.2e1 - t276 - t117 / 0.2e1 - t251 / 0.2e1 + t250 / 0.2e1 + t120 * mrSges(5,1) + t125 * mrSges(4,3) + (-Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * qJD(3)) * t152 + t303 * mrSges(7,3) - t63 * t34 - t90 * mrSges(4,2) + qJ(4) * t98 - t37 * t104 - t36 * t105; -t143 * mrSges(5,1) + t188 * t110 + t186 * t111 + t211 * t30 + t289 * t29 + t235 * t57 + t290 * t56 + (t104 * t188 - t105 * t186 + t123) * t153 + t231 * qJD(3) + (-qJD(3) * t53 - t303) * m(7) + (-qJD(3) * t84 - t153 * t213 + t216) * m(6) + (qJD(3) * t121 + t153 * t99 + t91) * m(5); -t133 * t104 - t132 * t105 - t225 * t56 - t78 * t57 + t14 + t98 + (-t225 * t6 - t5 * t78 + t46) * m(7) + (-t132 * t32 - t133 * t33 + t58) * m(6); -t141 - t53 * (-mrSges(7,1) * t78 + mrSges(7,2) * t225) + (Ifges(7,1) * t225 + t267) * t279 + t26 * t278 + (Ifges(7,5) * t225 + Ifges(7,6) * t78) * t275 - t5 * t56 + t6 * t57 + (t225 * t5 - t6 * t78) * mrSges(7,3) + (Ifges(7,2) * t78 + t27 + t76) * t281 + t287;];
tauc  = t7(:);
