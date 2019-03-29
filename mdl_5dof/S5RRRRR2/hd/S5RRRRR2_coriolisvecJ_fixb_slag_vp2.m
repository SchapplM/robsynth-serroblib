% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:25:52
% EndTime: 2019-03-29 15:25:56
% DurationCPUTime: 2.43s
% Computational Cost: add. (3540->330), mult. (7694->491), div. (0->0), fcn. (5196->8), ass. (0->175)
t163 = cos(qJ(4));
t164 = cos(qJ(3));
t159 = sin(qJ(4));
t160 = sin(qJ(3));
t222 = t159 * t160;
t143 = -t163 * t164 + t222;
t155 = qJD(1) + qJD(2);
t133 = t143 * t155;
t283 = t133 / 0.2e1;
t161 = sin(qJ(2));
t242 = pkin(1) * qJD(1);
t207 = t161 * t242;
t172 = qJD(3) * pkin(2) - t160 * t207;
t194 = t164 * t207;
t112 = t159 * t172 + t163 * t194;
t165 = cos(qJ(2));
t220 = qJD(1) * t165;
t206 = pkin(1) * t220;
t223 = t155 * t164;
t142 = -pkin(2) * t223 - t206;
t158 = sin(qJ(5));
t162 = cos(qJ(5));
t81 = -t112 * t158 + t142 * t162;
t282 = t81 * mrSges(6,1);
t82 = t112 * t162 + t142 * t158;
t281 = t82 * mrSges(6,2);
t154 = qJD(3) + qJD(4);
t280 = Ifges(5,5) * t154;
t279 = Ifges(5,2) * t133;
t278 = Ifges(5,6) * t154;
t144 = t159 * t164 + t160 * t163;
t134 = t144 * t155;
t105 = -t134 * t158 + t154 * t162;
t277 = Ifges(6,6) * t105;
t124 = qJD(5) + t133;
t276 = Ifges(6,3) * t124;
t275 = t154 * t161;
t107 = t154 * t143;
t212 = qJD(5) * t162;
t178 = -t107 * t158 + t144 * t212;
t219 = qJD(2) * t165;
t201 = t164 * t219;
t217 = qJD(3) * t160;
t127 = (-t161 * t217 + t201) * t242;
t202 = t160 * t219;
t216 = qJD(3) * t164;
t128 = (-t161 * t216 - t202) * t242;
t274 = t127 * t164 - t128 * t160;
t254 = pkin(1) * t161;
t205 = qJD(2) * t254;
t193 = qJD(1) * t205;
t204 = pkin(2) * t217;
t137 = t155 * t204 + t193;
t111 = t159 * t194 - t163 * t172;
t215 = qJD(4) * t111;
t61 = t163 * t127 + t159 * t128 - t215;
t15 = qJD(5) * t81 + t137 * t158 + t162 * t61;
t16 = -qJD(5) * t82 + t137 * t162 - t158 * t61;
t273 = t15 * t162 - t16 * t158;
t272 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t106 = t134 * t162 + t154 * t158;
t247 = Ifges(6,4) * t106;
t48 = Ifges(6,2) * t105 + Ifges(6,6) * t124 + t247;
t269 = -t48 / 0.2e1;
t92 = t107 * t155;
t57 = t105 * qJD(5) - t162 * t92;
t268 = t57 / 0.2e1;
t58 = -t106 * qJD(5) + t158 * t92;
t267 = t58 / 0.2e1;
t108 = t154 * t144;
t93 = t108 * t155;
t266 = t93 / 0.2e1;
t271 = -t158 * t81 + t162 * t82;
t270 = Ifges(6,1) * t268 + Ifges(6,4) * t267 + Ifges(6,5) * t266;
t265 = -t105 / 0.2e1;
t264 = -t106 / 0.2e1;
t263 = t106 / 0.2e1;
t262 = -t124 / 0.2e1;
t259 = t134 / 0.2e1;
t258 = t160 / 0.2e1;
t257 = t162 / 0.2e1;
t255 = m(4) * pkin(1) ^ 2;
t253 = pkin(1) * t165;
t251 = mrSges(5,3) * t133;
t250 = Ifges(4,4) * t160;
t248 = Ifges(5,4) * t134;
t246 = Ifges(6,4) * t158;
t245 = Ifges(6,4) * t162;
t244 = Ifges(4,5) * t164;
t243 = Ifges(4,6) * t160;
t241 = t134 * mrSges(5,3);
t236 = t164 * Ifges(4,2);
t214 = qJD(4) * t112;
t62 = t159 * t127 - t163 * t128 + t214;
t235 = t62 * t144;
t234 = mrSges(5,1) * t154 + mrSges(6,1) * t105 - mrSges(6,2) * t106 - t241;
t119 = t143 * t207;
t232 = t111 * t119;
t176 = pkin(1) * t144;
t120 = t176 * t220;
t231 = t111 * t120;
t228 = t133 * t158;
t227 = t133 * t162;
t226 = t144 * t158;
t225 = t144 * t162;
t224 = t155 * t160;
t221 = t160 ^ 2 + t164 ^ 2;
t218 = qJD(3) * t155;
t213 = qJD(5) * t158;
t211 = -qJD(1) - t155;
t210 = -qJD(2) + t155;
t209 = Ifges(6,5) * t57 + Ifges(6,6) * t58 + Ifges(6,3) * t93;
t208 = pkin(2) * t224;
t104 = Ifges(6,4) * t105;
t49 = Ifges(6,1) * t106 + Ifges(6,5) * t124 + t104;
t203 = t49 * t257;
t95 = mrSges(5,1) * t133 + mrSges(5,2) * t134;
t199 = -m(5) * t142 - t95;
t196 = -t213 / 0.2e1;
t191 = -mrSges(4,1) * t164 + mrSges(4,2) * t160;
t190 = mrSges(4,1) * t160 + mrSges(4,2) * t164;
t189 = mrSges(6,1) * t158 + mrSges(6,2) * t162;
t188 = Ifges(6,1) * t162 - t246;
t187 = -Ifges(6,2) * t158 + t245;
t186 = Ifges(6,5) * t162 - Ifges(6,6) * t158;
t131 = t161 * t176;
t69 = (t159 * t201 + (t164 * t275 + t202) * t163) * pkin(1) - t154 * t222 * t254;
t185 = t111 * t69 + t131 * t62;
t74 = -mrSges(6,2) * t124 + mrSges(6,3) * t105;
t75 = mrSges(6,1) * t124 - mrSges(6,3) * t106;
t184 = -t158 * t74 - t162 * t75;
t183 = t158 * t82 + t162 * t81;
t182 = -t111 * t107 + t235;
t132 = t143 * t254;
t151 = -pkin(2) * t164 - t253;
t100 = -t132 * t162 + t151 * t158;
t99 = t132 * t158 + t151 * t162;
t147 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t224;
t148 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t223;
t181 = t147 * t164 + t148 * t160;
t180 = -t147 * t160 + t148 * t164;
t109 = -mrSges(5,2) * t154 - t251;
t179 = -t158 * t75 + t162 * t74 + t109;
t177 = t107 * t162 + t144 * t213;
t175 = t160 * (Ifges(4,1) * t164 - t250);
t174 = (t236 + t250) * t155;
t173 = t143 * t165;
t129 = Ifges(4,6) * qJD(3) + t174;
t152 = Ifges(4,4) * t223;
t130 = Ifges(4,1) * t224 + Ifges(4,5) * qJD(3) + t152;
t47 = Ifges(6,5) * t106 + t276 + t277;
t6 = Ifges(6,4) * t57 + Ifges(6,2) * t58 + Ifges(6,6) * t93;
t78 = t248 + t278 - t279;
t123 = Ifges(5,4) * t133;
t79 = Ifges(5,1) * t134 - t123 + t280;
t169 = (-t177 * Ifges(6,1) - t178 * Ifges(6,4)) * t263 + t225 * t270 + t189 * t235 + t191 * t193 + t111 * (t178 * mrSges(6,1) - t177 * mrSges(6,2)) + t124 * (-t177 * Ifges(6,5) - t178 * Ifges(6,6)) / 0.2e1 + t105 * (-t177 * Ifges(6,4) - t178 * Ifges(6,2)) / 0.2e1 + t175 * t218 - t6 * t226 / 0.2e1 + qJD(3) ^ 2 * (-t243 + t244) / 0.2e1 + t178 * t269 - (t174 + t129) * t217 / 0.2e1 - (Ifges(5,1) * t259 + t280 / 0.2e1 + t142 * mrSges(5,2) + t79 / 0.2e1 + t203) * t107 + (t279 / 0.2e1 + Ifges(6,5) * t263 - t112 * mrSges(5,3) - t278 / 0.2e1 + t142 * mrSges(5,1) + t47 / 0.2e1 - t78 / 0.2e1 + t282 - t281 + t276 / 0.2e1 + t277 / 0.2e1) * t108 + (-t15 * t226 - t16 * t225 + t177 * t81 - t178 * t82) * mrSges(6,3) + t274 * mrSges(4,3) + (mrSges(5,2) * t137 - Ifges(5,1) * t92 + t186 * t266 + t187 * t267 + t188 * t268 + t196 * t49) * t144 + (Ifges(5,2) * t93 + Ifges(6,3) * t266 + Ifges(6,6) * t267 + Ifges(6,5) * t268 - t61 * mrSges(5,3) + t137 * mrSges(5,1) + t209 / 0.2e1 + t272) * t143 + (t107 * t283 - t108 * t259 + t143 * t92 - t93 * t144) * Ifges(5,4) + (t130 + (0.3e1 * Ifges(4,4) * t164 + (Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t160) * t155) * t216 / 0.2e1;
t168 = (-mrSges(6,1) * t162 + mrSges(6,2) * t158 - mrSges(5,1)) * t62 + (t111 * t189 + t203) * qJD(5) + (-t183 * qJD(5) - t227 * t81 - t228 * t82 + t273) * mrSges(6,3) + t78 * t259 + (Ifges(6,5) * t158 + Ifges(6,6) * t162) * t266 + (Ifges(6,2) * t162 + t246) * t267 + (Ifges(6,1) * t158 + t245) * t268 + t158 * t270 + t111 * t251 + t6 * t257 - Ifges(5,5) * t92 - Ifges(5,6) * t93 - t61 * mrSges(5,2) + (Ifges(6,3) * t134 - t186 * t133) * t262 + (Ifges(6,5) * t134 - t188 * t133) * t264 + (Ifges(6,6) * t134 - t187 * t133) * t265 - t154 * (-Ifges(5,5) * t133 - Ifges(5,6) * t134) / 0.2e1 - t142 * (mrSges(5,1) * t134 - mrSges(5,2) * t133) - t134 * t282 + t48 * t196 + t228 * t269 + t134 * t281 + (-Ifges(5,2) * t134 - t123 + t79) * t283 + t49 * t227 / 0.2e1 + (t105 * t187 + t106 * t188 + t124 * t186) * qJD(5) / 0.2e1 - (-Ifges(5,1) * t133 - t248 + t47) * t134 / 0.2e1;
t146 = t204 + t205;
t145 = t190 * qJD(3);
t138 = t191 * t155;
t122 = t173 * t242;
t121 = t144 * t207;
t103 = -t122 * t162 + t158 * t207;
t102 = t122 * t158 + t162 * t207;
t98 = -t121 * t162 + t158 * t208;
t97 = t121 * t158 + t162 * t208;
t80 = t189 * t133;
t68 = (-qJD(2) * t173 - t144 * t275) * pkin(1);
t43 = mrSges(5,1) * t93 - mrSges(5,2) * t92;
t25 = -t100 * qJD(5) + t146 * t162 - t158 * t68;
t24 = t99 * qJD(5) + t146 * t158 + t162 * t68;
t23 = -mrSges(6,2) * t93 + mrSges(6,3) * t58;
t22 = mrSges(6,1) * t93 - mrSges(6,3) * t57;
t12 = -mrSges(6,1) * t58 + mrSges(6,2) * t57;
t1 = [m(5) * (t112 * t68 - t132 * t61 + t137 * t151 + t142 * t146 + t185) + m(6) * (t100 * t15 + t16 * t99 + t24 * t82 + t25 * t81 + t185) + (-qJD(1) * t145 - t190 * t218 + (t211 * mrSges(3,2) + t180) * qJD(2)) * t253 + t146 * t95 + t151 * t43 + t131 * t12 + t68 * t109 + t99 * t22 + t100 * t23 + t24 * t74 + t25 * t75 + ((-0.2e1 + t221) * qJD(1) * t219 * t255 + (m(4) * t274 - t181 * qJD(3) + (t211 * mrSges(3,1) + t138) * qJD(2)) * pkin(1)) * t161 - t234 * t69 + t169 + (-t131 * t92 + t132 * t93 + t182) * mrSges(5,3); t122 * t109 - t102 * t75 - t103 * t74 + t182 * mrSges(5,3) + (-(-0.1e1 + t221) * qJD(1) ^ 2 * t161 * t255 + (t210 * mrSges(3,2) - t145 - t180) * t242) * t165 + ((m(6) * t183 - t184 - t199) * t217 + (-t162 * t22 - t158 * t23 - m(5) * t137 + m(6) * (-t15 * t158 - t16 * t162 - t82 * t212 + t81 * t213) - t43 - t74 * t212 + t75 * t213) * t164) * pkin(2) - m(5) * (-t112 * t122 + t142 * t207 + t231) - m(6) * (t102 * t81 + t103 * t82 + t231) + (t210 * mrSges(3,1) - t138 - t95) * t207 + t169 + t234 * t120; -m(6) * (t81 * t97 + t82 * t98 - t232) + t112 * t241 + t121 * t109 - t127 * mrSges(4,2) + t128 * mrSges(4,1) + t111 * t80 - t97 * t75 - t98 * t74 + (t199 * t224 + (t92 * mrSges(5,3) - t12 - m(6) * t62 + m(5) * (-t62 + t214) + (m(6) * t271 + t179) * qJD(4)) * t163 + (-t93 * mrSges(5,3) - t158 * t22 + t162 * t23 + t184 * qJD(5) - t234 * qJD(4) + m(6) * (-t81 * t212 - t82 * t213 + t215 + t273) + m(5) * (t61 + t215)) * t159) * pkin(2) + t168 + (t190 * t206 + t129 * t258 + (-t175 / 0.2e1 + t236 * t258) * t155 + (t244 / 0.2e1 - t243 / 0.2e1) * qJD(3) - (t130 + t152) * t164 / 0.2e1) * t155 - t234 * t119 - m(5) * (-t112 * t121 - t232) + t181 * t207; (t234 + t241) * t112 + (-m(6) * (t112 - t271) + t80 + t179) * t111 + t168; -t111 * (mrSges(6,1) * t106 + mrSges(6,2) * t105) + (Ifges(6,1) * t105 - t247) * t264 + t48 * t263 + (Ifges(6,5) * t105 - Ifges(6,6) * t106) * t262 - t81 * t74 + t82 * t75 + (t105 * t81 + t106 * t82) * mrSges(6,3) + t209 + (-Ifges(6,2) * t106 + t104 + t49) * t265 + t272;];
tauc  = t1(:);
