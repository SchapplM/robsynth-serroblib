% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:05
% EndTime: 2020-01-03 12:11:12
% DurationCPUTime: 2.73s
% Computational Cost: add. (3683->303), mult. (6435->412), div. (0->0), fcn. (3887->6), ass. (0->167)
t279 = Ifges(5,4) + Ifges(6,4);
t280 = Ifges(5,1) + Ifges(6,1);
t278 = Ifges(6,5) + Ifges(5,5);
t277 = Ifges(5,2) + Ifges(6,2);
t276 = Ifges(6,6) + Ifges(5,6);
t188 = sin(qJ(4));
t189 = sin(qJ(3));
t191 = cos(qJ(4));
t192 = cos(qJ(3));
t158 = -t188 * t189 + t191 * t192;
t185 = qJD(1) + qJD(2);
t136 = t158 * t185;
t282 = t279 * t136;
t159 = t188 * t192 + t189 * t191;
t137 = t159 * t185;
t281 = t279 * t137;
t184 = qJD(3) + qJD(4);
t275 = t277 * t136 + t276 * t184 + t281;
t274 = t280 * t137 + t278 * t184 + t282;
t264 = -pkin(8) - pkin(7);
t172 = t264 * t189;
t183 = t192 * pkin(8);
t173 = pkin(7) * t192 + t183;
t117 = t188 * t172 + t191 * t173;
t218 = qJD(3) * t264;
t163 = t189 * t218;
t164 = t192 * t218;
t193 = cos(qJ(2));
t239 = qJD(1) * pkin(1);
t222 = t193 * t239;
t271 = -qJD(4) * t117 + t159 * t222 - t163 * t188 + t191 * t164;
t225 = qJD(4) * t191;
t226 = qJD(4) * t188;
t269 = -t158 * t222 + t191 * t163 + t188 * t164 + t172 * t225 - t173 * t226;
t190 = sin(qJ(2));
t220 = t190 * t239;
t168 = pkin(7) * t185 + t220;
t219 = qJD(2) * t239;
t212 = t193 * t219;
t171 = t192 * t212;
t228 = qJD(3) * t189;
t108 = -t168 * t228 + t171;
t205 = t189 * t212;
t227 = qJD(3) * t192;
t109 = -t168 * t227 - t205;
t203 = t108 * t192 - t109 * t189;
t234 = qJ(5) * t136;
t215 = pkin(8) * t185 + t168;
t125 = t215 * t192;
t112 = t191 * t125;
t124 = t215 * t189;
t113 = qJD(3) * pkin(3) - t124;
t59 = t113 * t188 + t112;
t37 = t59 + t234;
t273 = mrSges(5,3) * t59 + mrSges(6,3) * t37;
t97 = t184 * t158;
t204 = -qJ(5) * t97 - qJD(5) * t159;
t272 = t204 + t271;
t98 = t184 * t159;
t235 = -t98 * qJ(5) + t158 * qJD(5);
t270 = t235 + t269;
t128 = t137 * qJ(5);
t110 = t188 * t125;
t58 = t191 * t113 - t110;
t36 = -t128 + t58;
t177 = pkin(1) * t190 + pkin(7);
t249 = -pkin(8) - t177;
t155 = t249 * t189;
t156 = t177 * t192 + t183;
t93 = t188 * t155 + t191 * t156;
t76 = t97 * t185;
t268 = t76 / 0.2e1;
t77 = t98 * t185;
t267 = -t77 / 0.2e1;
t263 = -t136 / 0.2e1;
t260 = t137 / 0.2e1;
t255 = t189 / 0.2e1;
t251 = pkin(1) * t193;
t250 = pkin(4) * t137;
t248 = mrSges(5,3) * t136;
t247 = mrSges(6,3) * t136;
t246 = Ifges(4,4) * t189;
t242 = Ifges(4,5) * t192;
t241 = Ifges(4,6) * t189;
t240 = pkin(1) * qJD(2);
t238 = t137 * mrSges(5,3);
t237 = t188 * t77;
t236 = t192 * Ifges(4,2);
t233 = qJ(5) * t159;
t230 = t185 * t189;
t229 = t185 * t192;
t61 = -t191 * t124 - t110;
t176 = t190 * t219;
t181 = pkin(3) * t228;
t146 = t185 * t181 + t176;
t224 = -qJD(1) - t185;
t223 = -qJD(2) + t185;
t221 = t193 * t240;
t180 = -t192 * pkin(3) - pkin(2);
t23 = t77 * mrSges(6,1) + t76 * mrSges(6,2);
t89 = pkin(4) * t98 + t181;
t216 = t227 / 0.2e1;
t214 = qJD(3) * t249;
t213 = t168 * (t189 ^ 2 + t192 ^ 2);
t60 = t124 * t188 - t112;
t92 = t191 * t155 - t156 * t188;
t116 = t191 * t172 - t173 * t188;
t31 = pkin(4) * t184 + t36;
t206 = qJD(3) * t215;
t90 = -t189 * t206 + t171;
t91 = -t192 * t206 - t205;
t10 = -qJD(4) * t59 - t188 * t90 + t191 * t91;
t4 = -qJ(5) * t76 - qJD(5) * t137 + t10;
t210 = -t4 * t159 - t31 * t97;
t209 = -t10 * t159 - t58 * t97;
t208 = -mrSges(4,1) * t192 + mrSges(4,2) * t189;
t207 = mrSges(4,1) * t189 + mrSges(4,2) * t192;
t166 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t230;
t167 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t229;
t202 = t166 * t192 + t167 * t189;
t201 = t189 * t166 - t192 * t167;
t133 = -t158 * pkin(4) + t180;
t9 = t113 * t225 - t125 * t226 + t188 * t91 + t191 * t90;
t200 = t189 * (Ifges(4,1) * t192 - t246);
t199 = (t236 + t246) * t185;
t198 = t207 * qJD(3);
t118 = t189 * t214 + t192 * t221;
t119 = -t189 * t221 + t192 * t214;
t20 = t191 * t118 + t188 * t119 + t155 * t225 - t156 * t226;
t139 = t180 * t185 - t222;
t21 = -qJD(4) * t93 - t118 * t188 + t191 * t119;
t3 = -qJ(5) * t77 + qJD(5) * t136 + t9;
t83 = -t136 * pkin(4) + qJD(5) + t139;
t195 = t10 * mrSges(5,1) + t4 * mrSges(6,1) - t9 * mrSges(5,2) - t3 * mrSges(6,2) - t139 * (mrSges(5,1) * t137 + mrSges(5,2) * t136) + t31 * t247 + t58 * t248 - t83 * (mrSges(6,1) * t137 + mrSges(6,2) * t136) - t276 * t77 + t278 * t76 - (t280 * t136 - t281) * t137 / 0.2e1 + t275 * t260 - (t278 * t136 - t276 * t137) * t184 / 0.2e1 + (-t277 * t137 + t274 + t282) * t263;
t134 = Ifges(4,6) * qJD(3) + t199;
t175 = Ifges(4,4) * t229;
t135 = Ifges(4,1) * t230 + Ifges(4,5) * qJD(3) + t175;
t169 = -t185 * pkin(2) - t222;
t55 = pkin(4) * t77 + t146;
t194 = qJD(3) ^ 2 * (-t241 + t242) / 0.2e1 + t169 * t198 + t208 * t176 + t135 * t216 + (t139 * mrSges(5,2) + t83 * mrSges(6,2)) * t97 - (-t139 * mrSges(5,1) - t83 * mrSges(6,1) + t273) * t98 + t274 * t97 / 0.2e1 - t275 * t98 / 0.2e1 + (-t277 * t98 + t279 * t97) * t136 / 0.2e1 + (-t279 * t98 + t280 * t97) * t260 + (-t277 * t77 + t279 * t76) * t158 / 0.2e1 + (-t279 * t77 + t280 * t76) * t159 / 0.2e1 + (-t276 * t98 + t278 * t97) * t184 / 0.2e1 - (t199 + t134) * t228 / 0.2e1 + (t146 * mrSges(5,2) + t55 * mrSges(6,2) + t279 * t267 + t280 * t268) * t159 + (-mrSges(5,1) * t146 - mrSges(6,1) * t55 + mrSges(5,3) * t9 + mrSges(6,3) * t3 + t277 * t267 + t279 * t268) * t158 + t203 * mrSges(4,3) + (qJD(3) * t200 + (0.3e1 * Ifges(4,4) * t192 + (Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t189) * t216) * t185;
t182 = t190 * t240;
t179 = -pkin(2) - t251;
t178 = pkin(3) * t191 + pkin(4);
t170 = t180 - t251;
t165 = t182 + t181;
t154 = t158 * qJ(5);
t147 = t208 * t185;
t138 = t185 * t198;
t121 = t133 - t251;
t107 = mrSges(5,1) * t184 - t238;
t106 = mrSges(6,1) * t184 - t137 * mrSges(6,3);
t105 = -mrSges(5,2) * t184 + t248;
t104 = -mrSges(6,2) * t184 + t247;
t99 = pkin(3) * t230 + t250;
t85 = t154 + t117;
t84 = t116 - t233;
t82 = -mrSges(5,1) * t136 + mrSges(5,2) * t137;
t81 = -mrSges(6,1) * t136 + mrSges(6,2) * t137;
t78 = t182 + t89;
t63 = t154 + t93;
t62 = t92 - t233;
t44 = -t128 + t61;
t43 = t60 - t234;
t24 = mrSges(5,1) * t77 + mrSges(5,2) * t76;
t7 = t204 + t21;
t6 = t20 + t235;
t1 = [t194 + t165 * t82 + t170 * t24 + t179 * t138 + t121 * t23 + t6 * t104 + t20 * t105 + t7 * t106 + t21 * t107 + t78 * t81 + ((t147 + m(4) * (qJD(1) * t179 + t169) + t224 * mrSges(3,1)) * t190 + (m(4) * t213 + mrSges(3,2) * t224 - t201) * t193) * t240 + (m(4) * t203 - t202 * qJD(3)) * t177 + m(5) * (t10 * t92 + t139 * t165 + t146 * t170 + t20 * t59 + t21 * t58 + t9 * t93) + m(6) * (t121 * t55 + t3 * t63 + t31 * t7 + t37 * t6 + t4 * t62 + t78 * t83) + (-t62 * t76 - t63 * t77 + t210) * mrSges(6,3) + (-t76 * t92 - t77 * t93 + t209) * mrSges(5,3); t194 + t180 * t24 - pkin(2) * t138 + t133 * t23 + t89 * t81 + (t189 * pkin(3) * t82 - pkin(7) * t202) * qJD(3) + t271 * t107 + t272 * t106 + t269 * t105 + t270 * t104 + (-t84 * t76 - t85 * t77 + t210) * mrSges(6,3) + (-t116 * t76 - t117 * t77 + t209) * mrSges(5,3) + ((mrSges(3,2) * t223 + t201) * t193 + (mrSges(3,1) * t223 - t147 - t81 - t82) * t190) * t239 + (t133 * t55 + t3 * t85 + t4 * t84 + (-t220 + t89) * t83 + t270 * t37 + t272 * t31) * m(6) + (t10 * t116 + t117 * t9 + t146 * t180 + t269 * t59 + t271 * t58 + (t181 - t220) * t139) * m(5) + (-pkin(2) * t176 + pkin(7) * t203 - (t169 * t190 + t193 * t213) * t239) * m(4); t195 + t202 * t168 + (-t169 * t207 + t134 * t255 + (t236 * t255 - t200 / 0.2e1) * t185 + (t242 / 0.2e1 - t241 / 0.2e1) * qJD(3) - (t175 + t135) * t192 / 0.2e1) * t185 - m(5) * (t58 * t60 + t59 * t61) - t99 * t81 - t44 * t104 - t61 * t105 - t43 * t106 - t60 * t107 - t108 * mrSges(4,2) + t109 * mrSges(4,1) + t59 * t238 + (-mrSges(6,3) * t237 - t82 * t230 + (-t191 * t76 - t237) * mrSges(5,3) + ((t104 + t105) * t191 + (-t106 - t107) * t188) * qJD(4) + (t10 * t191 - t139 * t230 + t188 * t9 + t225 * t59 - t226 * t58) * m(5)) * pkin(3) + (t137 * t37 - t178 * t76) * mrSges(6,3) + ((t188 * t3 + t225 * t37 - t226 * t31) * pkin(3) + t178 * t4 - t31 * t43 - t37 * t44 - t83 * t99) * m(6); t195 + (-pkin(4) * t81 + t273) * t137 - t36 * t104 - t58 * t105 + t37 * t106 + t59 * t107 - pkin(4) * t76 * mrSges(6,3) + (-t83 * t250 - (-t31 + t36) * t37 + t4 * pkin(4)) * m(6); -t136 * t104 + t137 * t106 + 0.2e1 * (t55 / 0.2e1 + t37 * t263 + t31 * t260) * m(6) + t23;];
tauc = t1(:);
