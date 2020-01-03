% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRP4
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:40
% EndTime: 2019-12-31 21:50:46
% DurationCPUTime: 2.40s
% Computational Cost: add. (3546->299), mult. (6204->401), div. (0->0), fcn. (3663->6), ass. (0->156)
t250 = Ifges(5,1) + Ifges(6,1);
t248 = Ifges(5,5) + Ifges(6,4);
t249 = -Ifges(5,4) + Ifges(6,5);
t157 = qJD(1) + qJD(2);
t160 = sin(qJ(4));
t161 = sin(qJ(3));
t163 = cos(qJ(4));
t164 = cos(qJ(3));
t129 = t160 * t161 - t163 * t164;
t156 = qJD(3) + qJD(4);
t86 = t156 * t129;
t69 = t86 * t157;
t241 = -t69 / 0.2e1;
t251 = mrSges(5,1) + mrSges(6,1);
t247 = Ifges(6,6) - Ifges(5,6);
t115 = t129 * t157;
t110 = Ifges(5,4) * t115;
t130 = t160 * t164 + t161 * t163;
t116 = t130 * t157;
t219 = Ifges(6,5) * t115;
t246 = t250 * t116 + t248 * t156 - t110 + t219;
t162 = sin(qJ(2));
t218 = pkin(1) * qJD(1);
t199 = t162 * t218;
t136 = pkin(7) * t157 + t199;
t165 = cos(qJ(2));
t217 = pkin(1) * qJD(2);
t196 = qJD(1) * t217;
t188 = t165 * t196;
t139 = t164 * t188;
t204 = qJD(3) * t161;
t93 = -t136 * t204 + t139;
t178 = t161 * t188;
t203 = qJD(3) * t164;
t94 = -t136 * t203 - t178;
t180 = -t161 * t94 + t164 * t93;
t239 = -pkin(8) - pkin(7);
t141 = t239 * t161;
t155 = t164 * pkin(8);
t142 = pkin(7) * t164 + t155;
t100 = t141 * t160 + t142 * t163;
t195 = qJD(3) * t239;
t132 = t161 * t195;
t187 = t164 * t195;
t198 = t165 * t218;
t244 = -qJD(4) * t100 + t130 * t198 - t132 * t160 + t163 * t187;
t174 = t163 * t141 - t142 * t160;
t243 = qJD(4) * t174 + t129 * t198 + t163 * t132 + t160 * t187;
t191 = pkin(8) * t157 + t136;
t106 = t191 * t164;
t210 = t106 * t160;
t105 = t191 * t161;
t96 = qJD(3) * pkin(3) - t105;
t53 = t163 * t96 - t210;
t242 = -t53 + qJD(5);
t149 = pkin(1) * t162 + pkin(7);
t227 = -pkin(8) - t149;
t127 = t227 * t161;
t128 = t149 * t164 + t155;
t177 = t163 * t127 - t128 * t160;
t179 = qJD(3) * t191;
t168 = -t164 * t179 - t178;
t206 = t163 * t106;
t54 = t160 * t96 + t206;
t80 = -t161 * t179 + t139;
t8 = qJD(4) * t54 + t160 * t80 - t163 * t168;
t238 = t8 * t177;
t237 = t8 * t174;
t236 = -t115 / 0.2e1;
t235 = t115 / 0.2e1;
t233 = t116 / 0.2e1;
t231 = -t156 / 0.2e1;
t230 = t156 / 0.2e1;
t229 = pkin(1) * t165;
t228 = t8 * t130;
t223 = mrSges(5,3) * t115;
t89 = -mrSges(5,2) * t156 - t223;
t92 = -t115 * mrSges(6,2) + mrSges(6,3) * t156;
t226 = t89 + t92;
t215 = t116 * mrSges(5,3);
t224 = mrSges(6,2) * t116;
t225 = t251 * t156 - t215 - t224;
t222 = Ifges(4,4) * t161;
t220 = Ifges(5,4) * t116;
t48 = -pkin(4) * t156 + t242;
t216 = t115 * t48;
t212 = Ifges(4,5) * qJD(3);
t211 = Ifges(4,6) * qJD(3);
t209 = t157 * t161;
t208 = t157 * t164;
t146 = t162 * t196;
t194 = t157 * t204;
t122 = pkin(3) * t194 + t146;
t202 = qJD(4) * t163;
t201 = -qJD(1) - t157;
t200 = -qJD(2) + t157;
t197 = t165 * t217;
t153 = pkin(3) * t204;
t152 = -t164 * pkin(3) - pkin(2);
t190 = qJD(3) * t227;
t189 = t136 * (t161 ^ 2 + t164 ^ 2);
t49 = qJ(5) * t156 + t54;
t7 = -qJD(4) * t210 + t160 * t168 + t163 * t80 + t96 * t202;
t5 = qJD(5) * t156 + t7;
t87 = t156 * t130;
t185 = -t5 * t129 - t49 * t87;
t184 = t53 * t86 + t228;
t70 = t87 * t157;
t82 = t127 * t160 + t128 * t163;
t183 = t177 * t69 - t70 * t82;
t182 = -t100 * t70 + t174 * t69;
t181 = -mrSges(4,1) * t164 + mrSges(4,2) * t161;
t71 = pkin(4) * t116 + qJ(5) * t115;
t134 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t209;
t135 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t208;
t176 = t134 * t164 + t135 * t161;
t175 = t161 * t134 - t164 * t135;
t173 = (t164 * Ifges(4,2) + t222) * t157;
t172 = (mrSges(4,1) * t161 + mrSges(4,2) * t164) * qJD(3);
t83 = t129 * pkin(4) - t130 * qJ(5) + t152;
t26 = pkin(4) * t87 + qJ(5) * t86 - qJD(5) * t130 + t153;
t118 = t152 * t157 - t198;
t171 = -t161 * t197 + t164 * t190;
t109 = Ifges(6,5) * t116;
t50 = t115 * pkin(4) - t116 * qJ(5) + t118;
t60 = Ifges(6,6) * t156 + Ifges(6,3) * t115 + t109;
t61 = -Ifges(5,2) * t115 + Ifges(5,6) * t156 + t220;
t167 = -t7 * mrSges(5,2) + t5 * mrSges(6,3) - t118 * (mrSges(5,1) * t116 - mrSges(5,2) * t115) - t53 * t223 + t49 * t224 - t50 * (mrSges(6,1) * t116 + mrSges(6,3) * t115) + t61 * t233 + (Ifges(6,3) * t116 - t219) * t236 - t251 * t8 + t247 * t70 - t248 * t69 + (-t248 * t115 + t247 * t116) * t231 + (-Ifges(5,2) * t116 - t110 + t246) * t235 - (-t250 * t115 + t109 - t220 + t60) * t116 / 0.2e1;
t113 = t173 + t211;
t144 = Ifges(4,4) * t208;
t114 = Ifges(4,1) * t209 + t144 + t212;
t137 = -t157 * pkin(2) - t198;
t9 = pkin(4) * t70 + qJ(5) * t69 - qJD(5) * t116 + t122;
t166 = mrSges(6,2) * t228 + t181 * t146 + t137 * t172 + (Ifges(4,1) * t164 - t222) * t194 + qJD(3) ^ 2 * (Ifges(4,5) * t164 - Ifges(4,6) * t161) / 0.2e1 + (t249 * t70 - t250 * t69) * t130 / 0.2e1 - (t173 + t113) * t204 / 0.2e1 + t180 * mrSges(4,3) + (t118 * mrSges(5,1) + t60 / 0.2e1 + t50 * mrSges(6,1) - t61 / 0.2e1 - t54 * mrSges(5,3) + Ifges(6,3) * t235 - Ifges(5,2) * t236 + t249 * t233 + t247 * t230) * t87 + (t122 * mrSges(5,2) - t9 * mrSges(6,3) + (-Ifges(5,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t70 + t250 * t241) * t130 + (t122 * mrSges(5,1) + t9 * mrSges(6,1) - t7 * mrSges(5,3) + (Ifges(6,3) + Ifges(5,2)) * t70 + 0.2e1 * t249 * t241) * t129 + (-t246 / 0.2e1 - t118 * mrSges(5,2) - t48 * mrSges(6,2) + t50 * mrSges(6,3) - Ifges(5,4) * t236 - Ifges(6,5) * t235 - t248 * t230 - t250 * t233) * t86 + (t114 + (0.3e1 * Ifges(4,4) * t164 + (Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t161) * t157) * t203 / 0.2e1;
t154 = t162 * t217;
t151 = -pkin(2) - t229;
t150 = -pkin(3) * t163 - pkin(4);
t147 = pkin(3) * t160 + qJ(5);
t145 = pkin(3) * t202 + qJD(5);
t138 = t152 - t229;
t133 = t154 + t153;
t123 = t181 * t157;
t117 = t157 * t172;
t101 = t161 * t190 + t164 * t197;
t79 = t83 - t229;
t75 = mrSges(5,1) * t115 + mrSges(5,2) * t116;
t74 = mrSges(6,1) * t115 - mrSges(6,3) * t116;
t57 = pkin(3) * t209 + t71;
t56 = -t105 * t163 - t210;
t55 = -t105 * t160 + t206;
t20 = t154 + t26;
t19 = mrSges(5,1) * t70 - mrSges(5,2) * t69;
t18 = mrSges(6,1) * t70 + mrSges(6,3) * t69;
t16 = qJD(4) * t82 + t101 * t160 - t163 * t171;
t15 = qJD(4) * t177 + t163 * t101 + t160 * t171;
t1 = [((t123 + m(4) * (qJD(1) * t151 + t137) + t201 * mrSges(3,1)) * t162 + (m(4) * t189 + mrSges(3,2) * t201 - t175) * t165) * t217 + (m(4) * t180 - qJD(3) * t176) * t149 + t226 * t15 - t225 * t16 + t151 * t117 + t133 * t75 + t138 * t19 + t20 * t74 + t79 * t18 + t166 + m(5) * (t118 * t133 + t122 * t138 + t15 * t54 - t16 * t53 + t7 * t82 - t238) + m(6) * (t15 * t49 + t16 * t48 + t20 * t50 + t5 * t82 + t79 * t9 - t238) + (t183 + t184) * mrSges(5,3) + (t183 + t185) * mrSges(6,2); t152 * t19 - pkin(2) * t117 + t26 * t74 + t83 * t18 + t166 + (t182 + t184) * mrSges(5,3) + (t161 * pkin(3) * t75 - pkin(7) * t176) * qJD(3) + ((mrSges(3,2) * t200 + t175) * t165 + (mrSges(3,1) * t200 - t123 - t74 - t75) * t162) * t218 + (t182 + t185) * mrSges(6,2) + t243 * t226 + t244 * t225 + (t100 * t5 + t83 * t9 - t237 + (-t199 + t26) * t50 + t243 * t49 - t244 * t48) * m(6) + (t100 * t7 + t122 * t152 - t237 + t243 * t54 + t244 * t53 + (t153 - t199) * t118) * m(5) + (-(t137 * t162 + t165 * t189) * t218 - pkin(2) * t146 + t180 * pkin(7)) * m(4); t176 * t136 - t226 * t56 + t225 * t55 + (-t147 * t70 - t150 * t69 + t216) * mrSges(6,2) + m(6) * (t145 * t49 + t147 * t5 + t150 * t8) + ((-t137 * mrSges(4,2) + t212 / 0.2e1 - t144 / 0.2e1 - t114 / 0.2e1) * t164 + (-t137 * mrSges(4,1) - t211 / 0.2e1 + t113 / 0.2e1 + (t222 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t164) * t157 + (-m(5) * t118 - t75) * pkin(3)) * t161) * t157 + (m(5) * (t160 * t7 - t163 * t8) + (-t160 * t70 + t163 * t69) * mrSges(5,3) + ((m(5) * t54 + t89) * t163 + (-m(5) * t53 + m(6) * t48 - t225) * t160) * qJD(4)) * pkin(3) - m(5) * (-t53 * t55 + t54 * t56) + t145 * t92 + t94 * mrSges(4,1) - t93 * mrSges(4,2) - t57 * t74 + t54 * t215 + t167 - m(6) * (t48 * t55 + t49 * t56 + t50 * t57); (t215 + t225) * t54 - t226 * t53 + (pkin(4) * t69 - qJ(5) * t70 + t216) * mrSges(6,2) + qJD(5) * t92 - t71 * t74 + t167 + (-t8 * pkin(4) + t5 * qJ(5) + t242 * t49 - t48 * t54 - t50 * t71) * m(6); -t69 * mrSges(6,2) + t116 * t74 - t156 * t92 + 0.2e1 * (t8 / 0.2e1 + t50 * t233 + t49 * t231) * m(6);];
tauc = t1(:);
