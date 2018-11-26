% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 15:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:51:10
% EndTime: 2018-11-23 15:51:14
% DurationCPUTime: 3.97s
% Computational Cost: add. (9050->404), mult. (20594->543), div. (0->0), fcn. (15108->8), ass. (0->189)
t162 = sin(pkin(10));
t163 = cos(pkin(10));
t167 = sin(qJ(4));
t170 = cos(qJ(4));
t139 = -t170 * t162 - t167 * t163;
t166 = sin(qJ(5));
t169 = cos(qJ(5));
t225 = t163 * t170;
t190 = t162 * t167 - t225;
t107 = t139 * t166 - t169 * t190;
t220 = qJD(4) * t170;
t221 = qJD(4) * t167;
t134 = -t162 * t221 + t163 * t220;
t135 = t139 * qJD(4);
t180 = qJD(5) * t107 + t169 * t134 + t135 * t166;
t132 = t139 * qJD(1);
t164 = -pkin(1) - qJ(3);
t145 = qJD(1) * t164 + qJD(2);
t212 = -pkin(7) * qJD(1) + t145;
t126 = t212 * t162;
t127 = t212 * t163;
t94 = t126 * t170 + t127 * t167;
t81 = pkin(8) * t132 + t94;
t232 = t166 * t81;
t222 = qJD(1) * t162;
t133 = qJD(1) * t225 - t167 * t222;
t93 = -t126 * t167 + t170 * t127;
t80 = -pkin(8) * t133 + t93;
t77 = qJD(4) * pkin(4) + t80;
t43 = t169 * t77 - t232;
t231 = t169 * t81;
t44 = t166 * t77 + t231;
t191 = t169 * t139 + t166 * t190;
t69 = qJD(5) * t191 - t134 * t166 + t169 * t135;
t280 = t180 * t44 + t43 * t69;
t124 = qJD(1) * t135;
t187 = qJD(1) * t190;
t75 = qJD(3) * t187 - qJD(4) * t94;
t177 = -pkin(8) * t124 + t75;
t125 = qJD(4) * t187;
t186 = t139 * qJD(3);
t74 = qJD(1) * t186 - t126 * t221 + t127 * t220;
t65 = pkin(8) * t125 + t74;
t13 = qJD(5) * t44 + t166 * t65 - t169 * t177;
t168 = cos(qJ(6));
t211 = t169 * t132 - t133 * t166;
t63 = qJD(5) * t211 + t124 * t169 + t125 * t166;
t159 = qJD(4) + qJD(5);
t165 = sin(qJ(6));
t192 = t132 * t166 + t169 * t133;
t84 = t159 * t168 - t165 * t192;
t39 = qJD(6) * t84 + t168 * t63;
t85 = t159 * t165 + t168 * t192;
t40 = -qJD(6) * t85 - t165 * t63;
t14 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t279 = m(7) * t13 + t14;
t199 = Ifges(7,5) * t168 - Ifges(7,6) * t165;
t240 = Ifges(7,4) * t168;
t201 = -Ifges(7,2) * t165 + t240;
t241 = Ifges(7,4) * t165;
t203 = Ifges(7,1) * t168 - t241;
t204 = mrSges(7,1) * t165 + mrSges(7,2) * t168;
t258 = t168 / 0.2e1;
t259 = -t165 / 0.2e1;
t266 = t85 / 0.2e1;
t256 = Ifges(7,4) * t85;
t96 = qJD(6) - t211;
t35 = Ifges(7,2) * t84 + Ifges(7,6) * t96 + t256;
t83 = Ifges(7,4) * t84;
t36 = Ifges(7,1) * t85 + Ifges(7,5) * t96 + t83;
t41 = -pkin(5) * t159 - t43;
t278 = t35 * t259 + t36 * t258 + t41 * t204 + t84 * t201 / 0.2e1 + t203 * t266 + t96 * t199 / 0.2e1;
t277 = -Ifges(6,2) * t211 / 0.2e1;
t42 = pkin(9) * t159 + t44;
t161 = qJD(1) * qJ(2);
t155 = qJD(3) + t161;
t144 = pkin(3) * t222 + t155;
t110 = -pkin(4) * t132 + t144;
t47 = -pkin(5) * t211 - pkin(9) * t192 + t110;
t16 = t165 * t47 + t168 * t42;
t15 = -t165 * t42 + t168 * t47;
t234 = t15 * t168;
t197 = t16 * t165 + t234;
t275 = t197 * mrSges(7,3);
t247 = -pkin(7) + t164;
t142 = t247 * t162;
t143 = t247 * t163;
t109 = t170 * t142 + t167 * t143;
t274 = -t15 * t165 + t16 * t168;
t223 = t162 ^ 2 + t163 ^ 2;
t209 = qJD(1) * t223;
t12 = qJD(5) * t43 + t166 * t177 + t169 * t65;
t160 = qJD(1) * qJD(2);
t111 = -pkin(4) * t125 + t160;
t64 = qJD(5) * t192 + t124 * t166 - t169 * t125;
t25 = pkin(5) * t64 - pkin(9) * t63 + t111;
t2 = qJD(6) * t15 + t12 * t168 + t165 * t25;
t3 = -qJD(6) * t16 - t12 * t165 + t168 * t25;
t273 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t39 + Ifges(7,6) * t40;
t95 = Ifges(6,4) * t211;
t272 = Ifges(6,1) * t192 / 0.2e1 + t95 / 0.2e1;
t67 = pkin(5) * t192 - pkin(9) * t211;
t271 = t39 / 0.2e1;
t270 = t40 / 0.2e1;
t269 = t64 / 0.2e1;
t268 = -t84 / 0.2e1;
t267 = -t85 / 0.2e1;
t265 = -t96 / 0.2e1;
t262 = -t133 / 0.2e1;
t261 = -t134 / 0.2e1;
t260 = t135 / 0.2e1;
t257 = m(3) * qJ(2);
t108 = -t142 * t167 + t170 * t143;
t91 = pkin(8) * t190 + t108;
t92 = pkin(8) * t139 + t109;
t50 = t166 * t92 - t169 * t91;
t252 = t13 * t50;
t251 = t168 * t2;
t250 = t3 * t165;
t249 = t43 * mrSges(6,3);
t248 = t44 * mrSges(6,3);
t246 = -mrSges(6,1) * t159 - mrSges(7,1) * t84 + mrSges(7,2) * t85 + mrSges(6,3) * t192;
t245 = m(4) * qJD(3);
t243 = Ifges(5,4) * t133;
t236 = t13 * t107;
t229 = t139 * t125;
t227 = t190 * t124;
t149 = t162 * pkin(3) + qJ(2);
t206 = mrSges(4,1) * t162 + mrSges(4,2) * t163;
t224 = -mrSges(5,1) * t132 + mrSges(5,2) * t133 + qJD(1) * t206;
t219 = qJD(6) * t165;
t218 = qJD(6) * t168;
t214 = t64 * mrSges(6,1) + t63 * mrSges(6,2);
t213 = -t125 * mrSges(5,1) + t124 * mrSges(5,2);
t119 = pkin(4) * t134 + qJD(2);
t117 = -pkin(4) * t139 + t149;
t207 = -t2 * t165 - t3 * t168;
t205 = mrSges(7,1) * t168 - mrSges(7,2) * t165;
t202 = Ifges(7,1) * t165 + t240;
t200 = Ifges(7,2) * t168 + t241;
t198 = Ifges(7,5) * t165 + Ifges(7,6) * t168;
t21 = mrSges(7,1) * t64 - mrSges(7,3) * t39;
t22 = -mrSges(7,2) * t64 + mrSges(7,3) * t40;
t195 = -t165 * t21 + t168 * t22;
t51 = t166 * t91 + t169 * t92;
t55 = -pkin(5) * t191 - pkin(9) * t107 + t117;
t27 = t165 * t55 + t168 * t51;
t26 = -t165 * t51 + t168 * t55;
t53 = -mrSges(7,2) * t96 + mrSges(7,3) * t84;
t54 = mrSges(7,1) * t96 - mrSges(7,3) * t85;
t194 = -t165 * t54 + t168 * t53;
t193 = -t165 * t53 - t168 * t54;
t89 = -mrSges(6,2) * t159 + mrSges(6,3) * t211;
t189 = -t194 - t89;
t182 = -qJD(6) * t197 - t250;
t181 = t134 * t94 + t135 * t93 - t139 * t74 - t190 * t75;
t87 = qJD(3) * t190 - qJD(4) * t109;
t86 = -t142 * t221 + t143 * t220 + t186;
t179 = -pkin(8) * t135 + t87;
t8 = t39 * Ifges(7,4) + t40 * Ifges(7,2) + t64 * Ifges(7,6);
t9 = t39 * Ifges(7,1) + t40 * Ifges(7,4) + t64 * Ifges(7,5);
t178 = -t12 * mrSges(6,2) + mrSges(7,3) * t251 + t202 * t271 + t200 * t270 + t198 * t269 + t165 * t9 / 0.2e1 - Ifges(6,6) * t64 + Ifges(6,5) * t63 + t8 * t258 + (-t205 - mrSges(6,1)) * t13 + t278 * qJD(6);
t176 = t110 * mrSges(6,1) + t15 * mrSges(7,1) - t16 * mrSges(7,2) - Ifges(6,4) * t192 + Ifges(7,5) * t85 - Ifges(6,6) * t159 + Ifges(7,6) * t84 + Ifges(7,3) * t96 + t277;
t175 = m(7) * (-t15 * t218 - t16 * t219 - t250 + t251) - t53 * t219 - t54 * t218 + t195;
t174 = t277 + t176;
t173 = t110 * mrSges(6,2) + Ifges(6,5) * t159 + t272 + t278;
t172 = t173 + t272;
t171 = qJD(1) ^ 2;
t128 = Ifges(5,4) * t132;
t116 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t133;
t115 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t132;
t98 = t133 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t128;
t97 = t132 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t243;
t73 = -pkin(8) * t134 + t86;
t66 = -mrSges(6,1) * t211 + mrSges(6,2) * t192;
t60 = Ifges(7,3) * t64;
t49 = pkin(4) * t133 + t67;
t46 = t169 * t80 - t232;
t45 = t166 * t80 + t231;
t28 = pkin(5) * t180 - pkin(9) * t69 + t119;
t24 = t165 * t67 + t168 * t43;
t23 = -t165 * t43 + t168 * t67;
t20 = t165 * t49 + t168 * t46;
t19 = -t165 * t46 + t168 * t49;
t18 = qJD(5) * t51 + t166 * t73 - t169 * t179;
t17 = -qJD(5) * t50 + t166 * t179 + t169 * t73;
t5 = -qJD(6) * t27 - t165 * t17 + t168 * t28;
t4 = qJD(6) * t26 + t165 * t28 + t168 * t17;
t1 = [m(5) * (t108 * t75 + t109 * t74 + t86 * t94 + t87 * t93) + (t8 * t259 + t9 * t258 + Ifges(6,1) * t63 - Ifges(6,4) * t64 + t111 * mrSges(6,2) + t203 * t271 + t201 * t270 + t199 * t269 + (mrSges(6,3) + t204) * t13 + t207 * mrSges(7,3) + (t36 * t259 - t168 * t35 / 0.2e1 + t41 * t205 + t200 * t268 + t202 * t267 + t198 * t265 - t274 * mrSges(7,3)) * qJD(6)) * t107 + (-t108 * t124 + t109 * t125 - t181) * mrSges(5,3) + (-t145 * t223 - t164 * t209) * t245 + (t133 * t260 - t227) * Ifges(5,1) + 0.2e1 * qJD(3) * mrSges(4,3) * t209 - (-t12 * mrSges(6,3) + t60 / 0.2e1 - Ifges(6,4) * t63 + t111 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t64 + t273) * t191 + (m(5) * (qJD(1) * t149 + t144) + m(4) * (t155 + t161) + t224 + (-mrSges(5,1) * t139 - mrSges(5,2) * t190 + (2 * mrSges(3,3)) + t206 + 0.2e1 * t257) * qJD(1)) * qJD(2) + (t139 * t124 - t125 * t190 + t132 * t260 + t133 * t261) * Ifges(5,4) + (t172 - t275) * t69 + t174 * t180 + t149 * t213 + t117 * t214 + t98 * t260 + t97 * t261 + t246 * t18 + m(6) * (t110 * t119 + t111 * t117 + t12 * t51 + t17 * t44 - t18 * t43 + t252) + m(7) * (t15 * t5 + t16 * t4 + t18 * t41 + t2 * t27 + t26 * t3 + t252) + t26 * t21 + t27 * t22 + (t132 * t261 + t229) * Ifges(5,2) + t50 * t14 + t4 * t53 + t5 * t54 + t17 * t89 + t86 * t115 + t87 * t116 + t119 * t66 + qJD(4) * (Ifges(5,5) * t135 - Ifges(5,6) * t134) / 0.2e1 + t144 * (mrSges(5,1) * t134 + mrSges(5,2) * t135) + (t50 * t63 - t51 * t64 - t280) * mrSges(6,3); t134 * t115 + t135 * t116 - t246 * t69 + (-mrSges(3,3) - t257) * t171 - (t63 * mrSges(6,3) + t14) * t107 + (t227 - t229) * mrSges(5,3) - t189 * t180 - (-t64 * mrSges(6,3) + qJD(6) * t193 + t195) * t191 + m(5) * t181 + m(6) * (-t12 * t191 - t236 + t280) + m(7) * (-t236 - t41 * t69 + t274 * t180 - (t182 + t251) * t191) + (-m(4) * t155 - m(5) * t144 - m(6) * t110 - m(7) * t197 - t223 * t245 + t193 - t224 - t66) * qJD(1); -t132 * t115 + t133 * t116 + t165 * t22 + t168 * t21 - t246 * t192 + t194 * qJD(6) + (m(4) + m(5)) * t160 - t223 * t171 * mrSges(4,3) + t189 * t211 - m(5) * (t132 * t94 - t133 * t93) + m(4) * t145 * t209 + t213 + t214 + (-t41 * t192 + t274 * t96 - t207) * m(7) + (t192 * t43 - t211 * t44 + t111) * m(6); t178 + (-t96 * t234 + (-t16 * t96 - t3) * t165) * mrSges(7,3) + (t132 * t93 + t133 * t94) * mrSges(5,3) + t279 * (-pkin(4) * t169 - pkin(5)) - m(6) * (-t43 * t45 + t44 * t46) - m(7) * (t15 * t19 + t16 * t20 + t41 * t45) - (-Ifges(5,2) * t133 + t128 + t98) * t132 / 0.2e1 + t175 * (pkin(4) * t166 + pkin(9)) + (-t174 + t248) * t192 + (-t172 + t249) * t211 + (Ifges(5,1) * t132 - t243) * t262 - t246 * t45 + (-t133 * t66 + (-t166 * t64 - t169 * t63) * mrSges(6,3) + (t246 * t166 - t189 * t169 + m(7) * (t166 * t41 + t169 * t274)) * qJD(5) + (0.2e1 * t110 * t262 + t12 * t166 - t13 * t169 + (-t166 * t43 + t169 * t44) * qJD(5)) * m(6)) * pkin(4) - t20 * t53 - t19 * t54 - t74 * mrSges(5,2) + t75 * mrSges(5,1) - t46 * t89 - t93 * t115 + t94 * t116 + Ifges(5,5) * t124 + Ifges(5,6) * t125 + t133 * t97 / 0.2e1 - qJD(4) * (Ifges(5,5) * t132 - Ifges(5,6) * t133) / 0.2e1 - t144 * (mrSges(5,1) * t133 + mrSges(5,2) * t132); t178 - m(7) * (t15 * t23 + t16 * t24 + t41 * t44) + (t249 - t95 / 0.2e1 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t192 + t275 - t173) * t211 + t175 * pkin(9) + (-t176 + t248) * t192 + t182 * mrSges(7,3) - t246 * t44 - t24 * t53 - t23 * t54 - t43 * t89 - t279 * pkin(5); t60 - t41 * (mrSges(7,1) * t85 + mrSges(7,2) * t84) + (Ifges(7,1) * t84 - t256) * t267 + t35 * t266 + (Ifges(7,5) * t84 - Ifges(7,6) * t85) * t265 - t15 * t53 + t16 * t54 + (t15 * t84 + t16 * t85) * mrSges(7,3) + (-Ifges(7,2) * t85 + t36 + t83) * t268 + t273;];
tauc  = t1(:);
