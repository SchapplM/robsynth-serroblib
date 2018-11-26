% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2018-11-23 15:48
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:48:37
% EndTime: 2018-11-23 15:48:41
% DurationCPUTime: 4.95s
% Computational Cost: add. (5029->433), mult. (10717->606), div. (0->0), fcn. (6391->8), ass. (0->208)
t149 = sin(qJ(4));
t144 = t149 * qJD(1);
t148 = sin(qJ(5));
t253 = -pkin(9) - pkin(8);
t188 = qJD(5) * t253;
t152 = cos(qJ(4));
t179 = pkin(4) * t152 + pkin(8) * t149;
t128 = t179 * qJD(1);
t151 = cos(qJ(5));
t139 = -cos(pkin(10)) * pkin(1) - pkin(2) - pkin(7);
t120 = qJD(1) * t139 + qJD(3);
t92 = -t149 * qJD(2) + t120 * t152;
t56 = t148 * t128 + t151 * t92;
t291 = -t56 + (-pkin(9) * t144 + t188) * t148;
t207 = t149 * t151;
t194 = pkin(9) * t207;
t55 = t151 * t128 - t148 * t92;
t290 = t151 * t188 - (pkin(5) * t152 + t194) * qJD(1) - t55;
t198 = t151 * qJD(4);
t206 = qJD(1) * t152;
t123 = -t148 * t206 + t198;
t140 = t144 + qJD(5);
t205 = qJD(2) * t152;
t93 = t149 * t120 + t205;
t82 = qJD(4) * pkin(8) + t93;
t141 = sin(pkin(10)) * pkin(1) + qJ(3);
t117 = t149 * pkin(4) - pkin(8) * t152 + t141;
t96 = t117 * qJD(1);
t47 = -t148 * t82 + t151 * t96;
t48 = t148 * t96 + t151 * t82;
t168 = t48 * t148 + t47 * t151;
t225 = Ifges(6,4) * t151;
t173 = -Ifges(6,2) * t148 + t225;
t226 = Ifges(6,4) * t148;
t175 = Ifges(6,1) * t151 - t226;
t176 = mrSges(6,1) * t148 + mrSges(6,2) * t151;
t223 = Ifges(6,6) * t148;
t224 = Ifges(6,5) * t151;
t238 = t151 / 0.2e1;
t239 = -t148 / 0.2e1;
t199 = t148 * qJD(4);
t124 = t151 * t206 + t199;
t244 = t124 / 0.2e1;
t227 = Ifges(6,4) * t124;
t60 = t123 * Ifges(6,2) + t140 * Ifges(6,6) + t227;
t118 = Ifges(6,4) * t123;
t61 = t124 * Ifges(6,1) + t140 * Ifges(6,5) + t118;
t81 = -qJD(4) * pkin(4) - t92;
t289 = -t168 * mrSges(6,3) + (-t223 + t224) * t140 / 0.2e1 + t173 * t123 / 0.2e1 + t175 * t244 + t81 * t176 + t238 * t61 + t239 * t60;
t196 = qJD(4) * qJD(5);
t200 = qJD(5) * t152;
t88 = t151 * t196 + (-t148 * t200 - t149 * t198) * qJD(1);
t159 = t149 * t199 - t151 * t200;
t89 = qJD(1) * t159 - t148 * t196;
t49 = -mrSges(6,1) * t89 + mrSges(6,2) * t88;
t84 = qJD(4) * t93;
t50 = -t89 * pkin(5) + t84;
t147 = sin(qJ(6));
t150 = cos(qJ(6));
t181 = t150 * t123 - t124 * t147;
t23 = qJD(6) * t181 + t147 * t89 + t150 * t88;
t68 = t123 * t147 + t124 * t150;
t24 = -qJD(6) * t68 - t147 * t88 + t150 * t89;
t8 = -mrSges(7,1) * t24 + mrSges(7,2) * t23;
t288 = m(6) * (-t198 * t48 + t199 * t47 + t84) + m(7) * t50 + t49 + t8;
t35 = pkin(9) * t123 + t48;
t216 = t150 * t35;
t34 = -pkin(9) * t124 + t47;
t32 = pkin(5) * t140 + t34;
t10 = t147 * t32 + t216;
t203 = qJD(4) * t152;
t184 = qJD(1) * t203;
t137 = Ifges(7,3) * t184;
t236 = Ifges(7,4) * t68;
t136 = qJD(6) + t140;
t243 = -t136 / 0.2e1;
t257 = -t68 / 0.2e1;
t259 = -t181 / 0.2e1;
t62 = Ifges(7,4) * t181;
t31 = Ifges(7,1) * t68 + Ifges(7,5) * t136 + t62;
t57 = -t123 * pkin(5) + t81;
t218 = t147 * t35;
t9 = t150 * t32 - t218;
t287 = t137 + (Ifges(7,5) * t181 - Ifges(7,6) * t68) * t243 + (t10 * t68 + t181 * t9) * mrSges(7,3) + (-Ifges(7,2) * t68 + t31 + t62) * t259 - t57 * (mrSges(7,1) * t68 + mrSges(7,2) * t181) + (Ifges(7,1) * t181 - t236) * t257;
t134 = t253 * t148;
t135 = t253 * t151;
t86 = t134 * t147 - t135 * t150;
t286 = -qJD(6) * t86 - t291 * t147 + t290 * t150;
t85 = t134 * t150 + t135 * t147;
t285 = qJD(6) * t85 + t290 * t147 + t291 * t150;
t165 = t147 * t148 - t150 * t151;
t268 = qJD(5) + qJD(6);
t284 = t268 * t165;
t212 = Ifges(5,5) * qJD(4);
t283 = t212 / 0.2e1 + (t152 * Ifges(5,1) - Ifges(5,4) * t149) * qJD(1) / 0.2e1 - t92 * mrSges(5,3) + t289;
t121 = qJD(4) * t179 + qJD(3);
t100 = t121 * qJD(1);
t83 = qJD(4) * t92;
t18 = -qJD(5) * t48 + t151 * t100 - t148 * t83;
t13 = pkin(5) * t184 - t88 * pkin(9) + t18;
t201 = qJD(5) * t151;
t202 = qJD(5) * t148;
t17 = t148 * t100 + t151 * t83 + t96 * t201 - t202 * t82;
t14 = pkin(9) * t89 + t17;
t2 = qJD(6) * t9 + t13 * t147 + t14 * t150;
t3 = -qJD(6) * t10 + t13 * t150 - t14 * t147;
t282 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t23 + Ifges(7,6) * t24;
t30 = Ifges(7,2) * t181 + Ifges(7,6) * t136 + t236;
t279 = t30 / 0.2e1;
t278 = -m(5) - m(4);
t185 = -Ifges(5,6) * qJD(4) / 0.2e1;
t251 = m(7) * t57;
t33 = -mrSges(7,1) * t181 + mrSges(7,2) * t68;
t273 = t33 + t251;
t106 = t165 * t152;
t126 = t147 * t151 + t148 * t150;
t114 = t126 * qJD(1);
t73 = t268 * t126;
t272 = -qJD(4) * t106 - t149 * t73 - t114;
t104 = t126 * t152;
t271 = t165 * qJD(1) - qJD(4) * t104 + t149 * t284;
t270 = t165 * t149;
t122 = t139 * t207;
t70 = t148 * t117 + t122;
t169 = -t148 * t18 + t151 * t17;
t90 = -mrSges(6,2) * t140 + mrSges(6,3) * t123;
t91 = mrSges(6,1) * t140 - mrSges(6,3) * t124;
t166 = -t148 * t90 - t151 * t91;
t266 = -m(6) * t168 + t166;
t63 = mrSges(6,1) * t184 - t88 * mrSges(6,3);
t64 = -mrSges(6,2) * t184 + t89 * mrSges(6,3);
t167 = -t148 * t63 + t151 * t64;
t265 = m(6) * (qJD(4) * t81 - t201 * t47 - t202 * t48 + t169) + qJD(4) * t251 + t166 * qJD(5) + t167;
t264 = t18 * mrSges(6,1) - t17 * mrSges(6,2) + Ifges(6,5) * t88 + Ifges(6,6) * t89 + t282;
t261 = t23 / 0.2e1;
t260 = t24 / 0.2e1;
t258 = t181 / 0.2e1;
t256 = t68 / 0.2e1;
t255 = t88 / 0.2e1;
t254 = t89 / 0.2e1;
t249 = -t104 / 0.2e1;
t248 = -t106 / 0.2e1;
t247 = -t123 / 0.2e1;
t245 = -t124 / 0.2e1;
t242 = t136 / 0.2e1;
t241 = -t140 / 0.2e1;
t235 = pkin(5) * t148;
t234 = pkin(9) * t152;
t98 = qJD(1) * t270;
t230 = -t284 - t98;
t97 = t149 * t114;
t229 = -t73 - t97;
t228 = Ifges(5,4) * t152;
t189 = mrSges(5,3) * t206;
t214 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t123 + mrSges(6,2) * t124 + t189;
t213 = mrSges(4,3) * qJD(1);
t208 = t139 * t148;
t131 = qJD(1) * t141;
t204 = qJD(4) * t149;
t197 = qJD(1) * qJD(3);
t193 = t33 + t214;
t192 = mrSges(5,3) * t144;
t191 = Ifges(5,4) * t144;
t190 = t149 * t208;
t187 = t139 * t203;
t186 = -t212 / 0.2e1;
t183 = pkin(5) - t208;
t180 = 0.2e1 * t131;
t177 = mrSges(6,1) * t151 - mrSges(6,2) * t148;
t174 = Ifges(6,1) * t148 + t225;
t172 = Ifges(6,2) * t151 + t226;
t170 = Ifges(6,5) * t148 + Ifges(6,6) * t151;
t102 = t151 * t117;
t54 = t149 * t183 - t151 * t234 + t102;
t58 = -t148 * t234 + t70;
t25 = -t147 * t58 + t150 * t54;
t26 = t147 * t54 + t150 * t58;
t132 = -qJD(4) * mrSges(5,2) - t192;
t164 = -t148 * t91 + t151 * t90 + t132;
t36 = -qJD(5) * t190 + t117 * t201 + t148 * t121 + t151 * t187;
t154 = t131 * mrSges(5,1) + t47 * mrSges(6,1) + t9 * mrSges(7,1) + t185 - (-Ifges(5,2) * t149 + t228) * qJD(1) / 0.2e1 + t136 * Ifges(7,3) + t68 * Ifges(7,5) + t181 * Ifges(7,6) + t140 * Ifges(6,3) + t124 * Ifges(6,5) + t123 * Ifges(6,6) - t10 * mrSges(7,2) - t48 * mrSges(6,2) - t93 * mrSges(5,3);
t142 = -pkin(5) * t151 - pkin(4);
t138 = Ifges(6,3) * t184;
t127 = (t149 * mrSges(5,1) + mrSges(5,2) * t152) * qJD(1);
t109 = t151 * t121;
t107 = (-t139 + t235) * t152;
t103 = t126 * t149;
t75 = -pkin(5) * t159 + t139 * t204;
t71 = t205 + (-qJD(1) * t235 + t120) * t149;
t69 = t102 - t190;
t53 = mrSges(7,1) * t136 - mrSges(7,3) * t68;
t52 = -mrSges(7,2) * t136 + mrSges(7,3) * t181;
t43 = Ifges(6,1) * t88 + Ifges(6,4) * t89 + Ifges(6,5) * t184;
t42 = Ifges(6,4) * t88 + Ifges(6,2) * t89 + Ifges(6,6) * t184;
t41 = t126 * t204 + t152 * t284;
t39 = qJD(4) * t270 - t152 * t73;
t37 = -qJD(5) * t70 - t148 * t187 + t109;
t28 = pkin(9) * t159 + t36;
t27 = t109 + (-t122 + (-t117 + t234) * t148) * qJD(5) + (t152 * t183 + t194) * qJD(4);
t20 = -mrSges(7,2) * t184 + t24 * mrSges(7,3);
t19 = mrSges(7,1) * t184 - t23 * mrSges(7,3);
t12 = t150 * t34 - t218;
t11 = -t147 * t34 - t216;
t7 = Ifges(7,1) * t23 + Ifges(7,4) * t24 + Ifges(7,5) * t184;
t6 = Ifges(7,4) * t23 + Ifges(7,2) * t24 + Ifges(7,6) * t184;
t5 = -qJD(6) * t26 - t147 * t28 + t150 * t27;
t4 = qJD(6) * t25 + t147 * t27 + t150 * t28;
t1 = [(mrSges(5,1) * t197 + t137 / 0.2e1 + t138 / 0.2e1 + (m(5) * t139 - mrSges(5,3)) * t83 + (0.3e1 / 0.2e1 * t191 + t186 - t180 * mrSges(5,2) + (-m(5) * t92 + m(6) * t81 + t214) * t139 - t283) * qJD(4) + t264) * t149 + m(7) * (t10 * t4 + t107 * t50 + t2 * t26 + t25 * t3 + t5 * t9 + t57 * t75) + (t42 * t239 + t43 * t238 + t175 * t255 + t173 * t254 + mrSges(5,2) * t197 + (-t148 * t17 - t151 * t18) * mrSges(6,3) + (mrSges(5,3) + t176) * t84 + (t172 * t247 + t174 * t245 + t170 * t241 + t81 * t177 + t61 * t239 - t151 * t60 / 0.2e1 + (t148 * t47 - t151 * t48) * mrSges(6,3)) * qJD(5) + (((-0.3e1 / 0.2e1 * Ifges(5,4) + t224 / 0.2e1 - t223 / 0.2e1) * t152 + Ifges(7,5) * t248 + Ifges(7,6) * t249 + t141 * mrSges(5,1) + (Ifges(7,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(5,2) + Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,1)) * t149) * qJD(1) + t154 + t185) * qJD(4) + (-t49 + (-m(5) - m(6)) * t84 + (m(5) * t93 + t132) * qJD(4)) * t139) * t152 + (-Ifges(7,4) * t106 - Ifges(7,2) * t104) * t260 + (-Ifges(7,1) * t106 - Ifges(7,4) * t104) * t261 + (t10 * t41 - t104 * t2 + t106 * t3 - t39 * t9) * mrSges(7,3) + t50 * (mrSges(7,1) * t104 - mrSges(7,2) * t106) + m(6) * (t17 * t70 + t18 * t69 + t48 * t36 + t47 * t37) + t41 * t279 + (-t180 * t278 + t127 + 0.2e1 * t213) * qJD(3) + t4 * t52 + t5 * t53 + t57 * (-mrSges(7,1) * t41 + mrSges(7,2) * t39) + t69 * t63 + t70 * t64 + t75 * t33 + t36 * t90 + t37 * t91 + t107 * t8 + (Ifges(7,5) * t39 + Ifges(7,6) * t41) * t242 + t7 * t248 + t6 * t249 + (Ifges(7,1) * t39 + Ifges(7,4) * t41) * t256 + (Ifges(7,4) * t39 + Ifges(7,2) * t41) * t258 + t25 * t19 + t26 * t20 + t39 * t31 / 0.2e1; m(7) * (t10 * t39 - t3 * t104 - t2 * t106 + t9 * t41) - t106 * t20 - t104 * t19 + t39 * t52 + t41 * t53 + ((-t164 - t192) * qJD(4) + t288) * t149 + ((-t189 + t193) * qJD(4) + t265) * t152; -t103 * t19 - t270 * t20 + t271 * t53 + t272 * t52 + (t131 * t278 - t127 - t213 + t266) * qJD(1) + (t10 * t272 - t3 * t103 - t2 * t270 + t271 * t9) * m(7) + (qJD(4) * t164 - t288) * t152 + (qJD(4) * t193 + t265) * t149; ((t131 * mrSges(5,2) + t186 - t191 / 0.2e1 + t283) * t149 + (-t154 + (t228 / 0.2e1 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t149) * qJD(1) + t185 + (Ifges(7,5) * t126 - Ifges(7,6) * t165 + t170) * qJD(4) / 0.2e1) * t152) * qJD(1) + (t10 * t229 - t126 * t3 - t165 * t2 - t230 * t9) * mrSges(7,3) + (Ifges(7,4) * t126 - Ifges(7,2) * t165) * t260 + (Ifges(7,1) * t126 - Ifges(7,4) * t165) * t261 + t50 * (mrSges(7,1) * t165 + mrSges(7,2) * t126) - t165 * t6 / 0.2e1 + (-t73 / 0.2e1 - t97 / 0.2e1) * t30 + (-Ifges(7,4) * t284 - Ifges(7,2) * t73) * t258 + (-t284 / 0.2e1 - t98 / 0.2e1) * t31 + (-Ifges(7,5) * t284 - Ifges(7,6) * t73) * t242 + (-Ifges(7,1) * t284 - Ifges(7,4) * t73) * t256 + t285 * t52 + t286 * t53 + (t285 * t10 + t142 * t50 + t2 * t86 + t286 * t9 + t3 * t85 - t57 * t71) * m(7) + (-mrSges(5,1) - t177) * t84 + (-mrSges(7,1) * t229 + mrSges(7,2) * t230) * t57 + t167 * pkin(8) + t169 * mrSges(6,3) - t214 * t93 + (pkin(8) * t266 + t235 * t273 + t289) * qJD(5) + (-pkin(4) * t84 + pkin(8) * t169 - t47 * t55 - t48 * t56 - t81 * t93) * m(6) - pkin(4) * t49 - t71 * t33 - t83 * mrSges(5,2) + t85 * t19 + t86 * t20 - t56 * t90 - t55 * t91 + t42 * t238 + (Ifges(7,5) * t98 + Ifges(7,6) * t97) * t243 + t172 * t254 + t174 * t255 + (Ifges(7,1) * t98 + Ifges(7,4) * t97) * t257 + (Ifges(7,4) * t98 + Ifges(7,2) * t97) * t259 + t126 * t7 / 0.2e1 - t92 * t132 + t142 * t8 + t148 * t43 / 0.2e1; (t147 * t20 + t150 * t19 + m(7) * (t147 * t2 + t150 * t3) - t273 * t124 + (-t147 * t53 + t150 * t52 + m(7) * (t10 * t150 - t147 * t9)) * qJD(6)) * pkin(5) + (-Ifges(6,2) * t124 + t118 + t61) * t247 - m(7) * (t10 * t12 + t11 * t9) + t68 * t279 + t264 + t138 + (t123 * t47 + t124 * t48) * mrSges(6,3) - t12 * t52 - t11 * t53 - t47 * t90 + t48 * t91 + (Ifges(6,5) * t123 - Ifges(6,6) * t124) * t241 + t60 * t244 + (Ifges(6,1) * t123 - t227) * t245 - t81 * (mrSges(6,1) * t124 + mrSges(6,2) * t123) + t287; t10 * t53 + t30 * t256 - t9 * t52 + t282 + t287;];
tauc  = t1(:);
