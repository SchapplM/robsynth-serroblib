% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:10
% EndTime: 2019-03-09 02:05:25
% DurationCPUTime: 7.81s
% Computational Cost: add. (3441->444), mult. (7053->614), div. (0->0), fcn. (3729->6), ass. (0->200)
t278 = Ifges(6,1) + Ifges(7,1);
t268 = Ifges(7,4) + Ifges(6,5);
t277 = Ifges(6,6) - Ifges(7,6);
t120 = sin(qJ(4));
t188 = qJD(1) * qJD(4);
t168 = t120 * t188;
t121 = cos(qJ(5));
t119 = sin(qJ(5));
t193 = qJD(5) * t120;
t122 = cos(qJ(4));
t195 = qJD(4) * t122;
t130 = -t119 * t193 + t121 * t195;
t187 = qJD(4) * qJD(5);
t65 = -qJD(1) * t130 + t121 * t187;
t190 = t119 * qJD(4);
t172 = t122 * t190;
t192 = qJD(5) * t121;
t129 = t120 * t192 + t172;
t66 = qJD(1) * t129 - t119 * t187;
t276 = (Ifges(6,4) - Ifges(7,5)) * t66 + t278 * t65 - t268 * t168;
t116 = sin(pkin(9));
t117 = cos(pkin(9));
t202 = qJD(1) * t122;
t175 = t117 * t202;
t207 = t121 * t122;
t179 = t116 * t207;
t194 = qJD(5) * t119;
t203 = qJD(1) * t121;
t275 = -qJD(5) * t179 + t117 * t194 + t119 * t175 + (t120 * t190 - t203) * t116;
t196 = qJD(4) * t121;
t171 = t120 * t196;
t209 = t119 * t122;
t84 = t116 * t209 + t117 * t121;
t274 = -qJD(5) * t84 - t116 * t171 - (t116 * t119 + t117 * t207) * qJD(1);
t204 = qJD(1) * t120;
t178 = mrSges(5,3) * t204;
t92 = t119 * t204 + t196;
t93 = t120 * t203 - t190;
t214 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t92 + mrSges(6,2) * t93 + t178;
t123 = -pkin(1) - pkin(2);
t104 = qJD(1) * t123 + qJD(2);
t205 = qJD(1) * t117;
t80 = qJ(2) * t205 + t116 * t104;
t73 = -qJD(1) * pkin(7) + t80;
t50 = qJD(3) * t122 - t120 * t73;
t273 = qJD(4) * t50;
t189 = qJD(1) * qJD(2);
t169 = t117 * t189;
t199 = qJD(3) * t120;
t51 = t122 * t73 + t199;
t211 = qJD(4) * t51;
t33 = t120 * t169 + t211;
t217 = t120 * t33;
t32 = t122 * t169 + t273;
t135 = t122 * t32 + t217;
t197 = qJD(4) * t120;
t272 = t50 * t195 + t51 * t197 - t135;
t219 = Ifges(7,5) * t121;
t221 = Ifges(6,4) * t121;
t271 = t119 * t278 - t219 + t221;
t252 = t66 / 0.2e1;
t108 = qJD(5) + t202;
t47 = qJD(4) * pkin(8) + t51;
t158 = pkin(4) * t122 + pkin(8) * t120;
t113 = t116 * qJ(2);
t79 = -qJD(1) * t113 + t104 * t117;
t72 = qJD(1) * pkin(3) - t79;
t48 = qJD(1) * t158 + t72;
t157 = -pkin(4) * t120 + pkin(8) * t122;
t201 = qJD(2) * t116;
t83 = qJD(4) * t157 + t201;
t74 = t83 * qJD(1);
t3 = t119 * t74 + t121 * t32 + t48 * t192 - t194 * t47;
t1 = -qJ(6) * t168 + qJD(6) * t108 + t3;
t270 = t1 * mrSges(7,3);
t269 = t3 * mrSges(6,2);
t267 = Ifges(7,2) + Ifges(6,3);
t238 = Ifges(7,5) * t92;
t91 = Ifges(6,4) * t92;
t266 = t268 * t108 - t278 * t93 - t238 + t91;
t265 = -m(3) * qJ(2) - mrSges(3,3);
t12 = -t119 * t47 + t121 * t48;
t264 = qJD(6) - t12;
t263 = t268 * t119 + t277 * t121;
t262 = -t277 * t119 + t268 * t121;
t220 = Ifges(7,5) * t119;
t222 = Ifges(6,4) * t119;
t261 = t278 * t121 + t220 - t222;
t200 = qJD(2) * t117;
t173 = t122 * t200;
t206 = t117 * qJ(2) + t116 * t123;
t95 = -pkin(7) + t206;
t260 = -t95 * t197 + t173;
t13 = t119 * t48 + t121 * t47;
t4 = -qJD(5) * t13 - t119 * t32 + t121 * t74;
t155 = -t119 * t4 + t121 * t3;
t2 = pkin(5) * t168 - t4;
t156 = t1 * t121 + t119 * t2;
t191 = qJD(5) * t122;
t163 = t117 * t123 - t113;
t94 = pkin(3) - t163;
t75 = t158 + t94;
t9 = -t119 * (qJD(5) * t75 + t260) - t121 * (t191 * t95 - t83);
t257 = -t168 * t267 + t268 * t65 + t277 * t66;
t137 = -t13 * t119 - t12 * t121;
t10 = -pkin(5) * t108 + t264;
t11 = qJ(6) * t108 + t13;
t138 = t10 * t121 - t11 * t119;
t142 = Ifges(7,3) * t119 + t219;
t148 = -Ifges(6,2) * t119 + t221;
t46 = -qJD(4) * pkin(4) - t50;
t15 = -pkin(5) * t92 + qJ(6) * t93 + t46;
t153 = mrSges(7,1) * t119 - mrSges(7,3) * t121;
t154 = mrSges(6,1) * t119 + mrSges(6,2) * t121;
t90 = Ifges(7,5) * t93;
t34 = t108 * Ifges(7,6) - t92 * Ifges(7,3) - t90;
t218 = t119 * t34;
t244 = t121 / 0.2e1;
t246 = t108 / 0.2e1;
t249 = -t93 / 0.2e1;
t250 = t92 / 0.2e1;
t251 = -t92 / 0.2e1;
t239 = Ifges(6,4) * t93;
t37 = t92 * Ifges(6,2) + t108 * Ifges(6,6) - t239;
t256 = t142 * t251 + t148 * t250 + t15 * t153 + t46 * t154 + t246 * t262 + t249 * t261 + t138 * mrSges(7,2) + t137 * mrSges(6,3) + t218 / 0.2e1 - t119 * t37 / 0.2e1 + t266 * t244;
t255 = -t65 * Ifges(7,5) / 0.2e1 + Ifges(7,6) * t168 / 0.2e1 + Ifges(7,3) * t252;
t254 = t65 / 0.2e1;
t253 = -t66 / 0.2e1;
t247 = -t108 / 0.2e1;
t241 = mrSges(6,3) * t92;
t240 = mrSges(6,3) * t93;
t5 = -pkin(5) * t66 - qJ(6) * t65 + qJD(6) * t93 + t33;
t234 = t120 * t5;
t22 = -mrSges(7,1) * t66 - mrSges(7,3) * t65;
t23 = -mrSges(6,1) * t66 + mrSges(6,2) * t65;
t231 = t22 + t23;
t41 = mrSges(7,2) * t66 - mrSges(7,3) * t168;
t44 = mrSges(6,2) * t168 + mrSges(6,3) * t66;
t230 = t41 + t44;
t42 = -mrSges(6,1) * t168 - mrSges(6,3) * t65;
t43 = mrSges(7,1) * t168 + t65 * mrSges(7,2);
t229 = -t42 + t43;
t97 = t157 * qJD(1);
t25 = t119 * t97 + t121 * t50;
t68 = -mrSges(6,2) * t108 + t241;
t71 = mrSges(7,2) * t92 + mrSges(7,3) * t108;
t226 = t68 + t71;
t69 = mrSges(6,1) * t108 + t240;
t70 = -mrSges(7,1) * t108 - mrSges(7,2) * t93;
t225 = -t70 + t69;
t31 = t119 * t75 + t95 * t207;
t224 = Ifges(5,4) * t120;
t223 = Ifges(5,4) * t122;
t216 = t121 * t75;
t213 = Ifges(5,5) * qJD(4);
t212 = Ifges(5,6) * qJD(4);
t208 = t120 * t121;
t186 = t119 * t83 + t121 * t173 + t75 * t192;
t53 = -mrSges(7,1) * t92 + mrSges(7,3) * t93;
t185 = -t53 + t214;
t180 = t95 * t195;
t177 = mrSges(5,3) * t202;
t176 = t117 * t204;
t164 = t193 / 0.2e1;
t147 = Ifges(6,2) * t121 + t222;
t141 = -Ifges(7,3) * t121 + t220;
t140 = -pkin(5) * t121 - qJ(6) * t119;
t139 = pkin(5) * t119 - qJ(6) * t121;
t24 = -t119 * t50 + t121 * t97;
t134 = -t139 + t95;
t96 = (mrSges(5,1) * t122 - mrSges(5,2) * t120) * qJD(1);
t133 = (-mrSges(5,1) * t120 - mrSges(5,2) * t122) * qJD(4);
t132 = (-t120 * Ifges(5,1) - t223) * qJD(1);
t131 = (-t122 * Ifges(5,2) - t224) * qJD(1);
t128 = t119 * t229 + t121 * t230;
t127 = -t119 * t226 - t121 * t225;
t124 = qJD(1) ^ 2;
t103 = -qJD(4) * mrSges(5,2) - t177;
t99 = -pkin(4) + t140;
t89 = qJD(1) * t133;
t87 = t132 + t213;
t86 = t131 + t212;
t85 = -t117 * t119 + t179;
t82 = qJD(5) * t139 - qJD(6) * t119;
t52 = -pkin(5) * t93 - qJ(6) * t92;
t49 = t134 * t120;
t36 = -t93 * Ifges(7,4) + t108 * Ifges(7,2) - t92 * Ifges(7,6);
t35 = -t93 * Ifges(6,5) + t92 * Ifges(6,6) + t108 * Ifges(6,3);
t30 = -t209 * t95 + t216;
t28 = t199 + (-qJD(1) * t139 + t73) * t122;
t27 = -t216 + (t119 * t95 - pkin(5)) * t122;
t26 = qJ(6) * t122 + t31;
t21 = pkin(5) * t204 - t24;
t20 = -qJ(6) * t204 + t25;
t17 = t65 * Ifges(6,4) + t66 * Ifges(6,2) - Ifges(6,6) * t168;
t14 = t134 * t195 + (qJD(5) * t140 + qJD(6) * t121 + t200) * t120;
t8 = (-t119 * t191 - t171) * t95 + t186;
t7 = pkin(5) * t197 - t9;
t6 = (-t194 * t95 + qJD(6)) * t122 + (-t121 * t95 - qJ(6)) * t197 + t186;
t16 = [t257 * t122 / 0.2e1 + 0.2e1 * (mrSges(4,1) * t116 - t265) * t189 + qJD(4) ^ 2 * (-Ifges(5,5) * t122 + Ifges(5,6) * t120) / 0.2e1 + t120 * t95 * t23 - t153 * t234 - (-Ifges(5,1) * t122 + t224) * t168 - t154 * t217 + t94 * t89 + t2 * (-mrSges(7,1) * t122 - mrSges(7,2) * t208) + t4 * (mrSges(6,1) * t122 + mrSges(6,3) * t208) + t11 * (mrSges(7,2) * t129 - mrSges(7,3) * t197) + t12 * (-mrSges(6,1) * t197 + mrSges(6,3) * t130) + t13 * (mrSges(6,2) * t197 + mrSges(6,3) * t129) + t10 * (mrSges(7,1) * t197 - mrSges(7,2) * t130) + m(6) * (t46 * t180 + t12 * t9 + t13 * t8 + t3 * t31 + t30 * t4 + (t200 * t46 + t33 * t95) * t120) + 0.2e1 * t96 * t201 + (t263 * t193 + (-t120 * t267 - t122 * t262) * qJD(4)) * t246 + (-t120 * t261 + t122 * t268) * t254 + t8 * t68 + t9 * t69 + t7 * t70 + t6 * t71 + t26 * t41 + t30 * t42 + t27 * t43 + t31 * t44 + t49 * t22 + t14 * t53 + m(7) * (t1 * t26 + t10 * t7 + t11 * t6 + t14 * t15 + t2 * t27 + t49 * t5) - t122 * t269 + (t131 + t86) * t197 / 0.2e1 + t15 * (-mrSges(7,1) * t129 + mrSges(7,3) * t130) + t46 * (-mrSges(6,1) * t129 - mrSges(6,2) * t130) + (t271 * t193 + (-t120 * t268 - t122 * t261) * qJD(4)) * t249 + t272 * mrSges(5,3) + ((t17 / 0.2e1 + t3 * mrSges(6,3) + t1 * mrSges(7,2) + t255) * t120 + t266 * t164) * t119 + 0.2e1 * mrSges(4,2) * t169 + t122 * t270 + (t147 * t193 + (-Ifges(6,6) * t120 - t122 * t148) * qJD(4)) * t250 + (t141 * t193 + (-Ifges(7,6) * t120 - t122 * t142) * qJD(4)) * t251 + (Ifges(6,6) * t122 - t120 * t148) * t252 + (Ifges(7,6) * t122 - t120 * t142) * t253 + m(4) * (-t79 * t116 + t80 * t117 + (-t116 * t163 + t117 * t206) * qJD(1)) * qJD(2) - t122 * (Ifges(5,2) * t120 - t223) * t188 + t214 * (-t120 * t200 - t180) - (qJD(5) * t34 + t276) * t208 / 0.2e1 + m(5) * (((-t51 * t120 - t50 * t122) * qJD(4) + t135) * t95 + ((-t120 * t50 + t122 * t51) * t117 + (qJD(1) * t94 + t72) * t116) * qJD(2)) + t72 * t133 - (t121 * t266 + t132 + t218 + t87) * t195 / 0.2e1 - (t36 + t35 + (-t120 * t262 + t122 * t267) * qJD(1)) * t197 / 0.2e1 + t260 * t103 + (t121 * t164 + t172 / 0.2e1) * t37; t230 * t85 + t229 * t84 + t265 * t124 + (-t124 * mrSges(4,2) - t89 + (-t122 * t103 + t120 * t185) * qJD(1)) * t117 - m(5) * (t175 * t51 - t176 * t50) - m(4) * t80 * t205 + t274 * t226 + t275 * t225 + (t1 * t85 - t275 * t10 + t274 * t11 - t15 * t176 + t2 * t84) * m(7) + (t275 * t12 + t274 * t13 - t176 * t46 + t3 * t85 - t4 * t84) * m(6) + (-t124 * mrSges(4,1) + t231 * t120 + (-t103 * t120 - t122 * t185) * qJD(4) + m(5) * (-t169 - t272) + m(7) * (t15 * t195 + t234) + m(6) * (t195 * t46 + t217) + (m(4) * t79 - m(5) * t72 - t96) * qJD(1)) * t116; ((-t119 * t225 + t121 * t226 + t103 + t177) * qJD(4) + m(7) * (t10 * t190 + t11 * t196 - t5) + m(5) * (-t33 + t211) + m(6) * (-t12 * t190 + t13 * t196 - t33) - t231) * t122 + (t127 * qJD(5) + (t178 - t185) * qJD(4) + m(7) * (qJD(4) * t15 + t10 * t192 - t11 * t194 + t156) + m(5) * (t32 - t273) + m(6) * (qJD(4) * t46 - t12 * t192 - t13 * t194 + t155) + t128) * t120; ((m(6) * t137 + m(7) * t138 + t127) * pkin(8) + t256) * qJD(5) + m(7) * (pkin(8) * t156 + t15 * t82 + t5 * t99) + m(6) * (-pkin(4) * t33 + pkin(8) * t155) + t5 * (-mrSges(7,1) * t121 - mrSges(7,3) * t119) + t99 * t22 - t50 * t103 + t214 * t51 + t128 * pkin(8) + (-mrSges(6,1) * t121 + mrSges(6,2) * t119 - mrSges(5,1)) * t33 - t25 * t68 - t24 * t69 - t21 * t70 - t20 * t71 + (t82 - t28) * t53 - t32 * mrSges(5,2) - pkin(4) * t23 - m(6) * (t12 * t24 + t13 * t25 + t46 * t51) - m(7) * (t10 * t21 + t11 * t20 + t15 * t28) + t271 * t254 + t121 * t255 + t147 * t252 + t141 * t253 + t17 * t244 + t276 * t119 / 0.2e1 + t155 * mrSges(6,3) + t156 * mrSges(7,2) + ((t212 / 0.2e1 + t35 / 0.2e1 + t36 / 0.2e1 - t86 / 0.2e1 + t12 * mrSges(6,1) - t13 * mrSges(6,2) - t10 * mrSges(7,1) + t11 * mrSges(7,3) + t72 * mrSges(5,1) - t51 * mrSges(5,3) + Ifges(5,4) * t204 / 0.2e1 + (-Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1) * t93 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t92 + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t108 - t263 * qJD(4) / 0.2e1) * t120 + (t87 / 0.2e1 + (-t223 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t120) * qJD(1) - t50 * mrSges(5,3) + t72 * mrSges(5,2) - t213 / 0.2e1 + t256) * t122) * qJD(1); (-t226 + t241) * t12 + (t225 - t240) * t13 - t46 * (-mrSges(6,1) * t93 + mrSges(6,2) * t92) - t15 * (-mrSges(7,1) * t93 - mrSges(7,3) * t92) + t257 + qJD(6) * t71 + qJ(6) * t41 - pkin(5) * t43 - t52 * t53 + t270 - t2 * mrSges(7,1) - t269 + t4 * mrSges(6,1) + (t278 * t92 + t239 + t34 - t90) * t93 / 0.2e1 + (-t10 * t92 - t11 * t93) * mrSges(7,2) + (Ifges(6,2) * t93 + t266 + t91) * t251 + (-pkin(5) * t2 + qJ(6) * t1 - t10 * t13 + t11 * t264 - t15 * t52) * m(7) + (t268 * t92 + t277 * t93) * t247 + t37 * t249 + (-Ifges(7,3) * t93 + t238) * t250; -t108 * t71 - t93 * t53 + 0.2e1 * (t2 / 0.2e1 + t11 * t247 + t15 * t249) * m(7) + t43;];
tauc  = t16(:);
