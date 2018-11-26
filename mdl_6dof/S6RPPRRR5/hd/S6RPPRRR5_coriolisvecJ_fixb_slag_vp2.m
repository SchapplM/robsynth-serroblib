% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2018-11-23 15:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:49:48
% EndTime: 2018-11-23 15:49:51
% DurationCPUTime: 2.93s
% Computational Cost: add. (4386->388), mult. (9187->541), div. (0->0), fcn. (5555->6), ass. (0->188)
t148 = sin(qJ(6));
t150 = cos(qJ(6));
t149 = sin(qJ(4));
t151 = cos(qJ(4));
t253 = sin(qJ(5));
t254 = cos(qJ(5));
t111 = t149 * t254 + t151 * t253;
t142 = qJD(4) + qJD(5);
t75 = t142 * t111 * qJD(1);
t191 = qJD(1) * t253;
t182 = t149 * t191;
t192 = qJD(1) * t254;
t106 = t151 * t192 - t182;
t89 = -t106 * t148 + t142 * t150;
t44 = qJD(6) * t89 - t150 * t75;
t187 = t254 * qJD(5);
t272 = t254 * qJD(4) + t187;
t157 = t272 * t151;
t76 = qJD(1) * t157 - t142 * t182;
t15 = mrSges(7,1) * t76 - mrSges(7,3) * t44;
t90 = t106 * t150 + t142 * t148;
t45 = -qJD(6) * t90 + t148 * t75;
t16 = -mrSges(7,2) * t76 + mrSges(7,3) * t45;
t168 = -t148 * t15 + t150 * t16;
t209 = qJD(6) * t150;
t210 = qJD(6) * t148;
t105 = -t149 * t192 - t151 * t191;
t102 = qJD(6) - t105;
t56 = -mrSges(7,2) * t102 + mrSges(7,3) * t89;
t57 = mrSges(7,1) * t102 - mrSges(7,3) * t90;
t274 = -t57 * t209 - t56 * t210 + t168;
t175 = mrSges(7,1) * t148 + mrSges(7,2) * t150;
t140 = qJD(1) * qJ(2) + qJD(3);
t121 = -qJD(1) * pkin(7) + t140;
t99 = (-pkin(8) * qJD(1) + t121) * t149;
t201 = t253 * t99;
t113 = t151 * t121;
t216 = qJD(1) * t151;
t100 = -pkin(8) * t216 + t113;
t95 = qJD(4) * pkin(4) + t100;
t59 = t254 * t95 - t201;
t53 = -t142 * pkin(5) - t59;
t162 = t53 * t175;
t170 = Ifges(7,5) * t150 - Ifges(7,6) * t148;
t238 = Ifges(7,4) * t150;
t172 = -Ifges(7,2) * t148 + t238;
t239 = Ifges(7,4) * t148;
t174 = Ifges(7,1) * t150 - t239;
t255 = t150 / 0.2e1;
t259 = t102 / 0.2e1;
t261 = t90 / 0.2e1;
t263 = t89 / 0.2e1;
t252 = Ifges(7,4) * t90;
t37 = Ifges(7,2) * t89 + Ifges(7,6) * t102 + t252;
t87 = Ifges(7,4) * t89;
t38 = Ifges(7,1) * t90 + Ifges(7,5) * t102 + t87;
t273 = t170 * t259 + t172 * t263 + t174 * t261 + t38 * t255 + t162 - t148 * t37 / 0.2e1;
t234 = t106 * mrSges(6,3);
t243 = -mrSges(6,1) * t142 - mrSges(7,1) * t89 + mrSges(7,2) * t90 + t234;
t165 = t148 * t57 - t150 * t56;
t242 = mrSges(6,3) * t105;
t91 = -mrSges(6,2) * t142 + t242;
t163 = -t165 + t91;
t206 = qJD(1) * qJD(4);
t186 = t151 * t206;
t189 = qJD(5) * t253;
t271 = -qJD(4) * t253 - t189;
t202 = t254 * t99;
t60 = t253 * t95 + t202;
t54 = t142 * pkin(9) + t60;
t147 = pkin(1) + qJ(3);
t268 = qJD(1) * t147;
t122 = -qJD(2) + t268;
t217 = qJD(1) * t149;
t108 = pkin(4) * t217 + t122;
t55 = -pkin(5) * t105 - pkin(9) * t106 + t108;
t17 = -t148 * t54 + t150 * t55;
t18 = t148 * t55 + t150 * t54;
t269 = -t148 * t17 + t150 * t18;
t212 = qJD(4) * t149;
t194 = t121 * t212;
t214 = qJD(2) * t151;
t154 = -t194 + (pkin(8) * t212 + t214) * qJD(1);
t208 = qJD(1) * qJD(2);
t211 = qJD(4) * t151;
t97 = t121 * t211 + t149 * t208;
t88 = -pkin(8) * t186 + t97;
t13 = t59 * qJD(5) + t154 * t253 + t254 * t88;
t123 = pkin(4) * t211 + qJD(3);
t114 = t123 * qJD(1);
t26 = pkin(5) * t76 + pkin(9) * t75 + t114;
t2 = qJD(6) * t17 + t13 * t150 + t148 * t26;
t3 = -qJD(6) * t18 - t13 * t148 + t150 * t26;
t267 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t44 + Ifges(7,6) * t45;
t250 = t148 * t3;
t266 = -t17 * t209 - t18 * t210 - t250;
t265 = (t149 ^ 2 + t151 ^ 2) * t121;
t152 = qJD(1) ^ 2;
t264 = -t89 / 0.2e1;
t262 = -t90 / 0.2e1;
t260 = -t102 / 0.2e1;
t256 = t148 / 0.2e1;
t14 = t60 * qJD(5) - t154 * t254 + t253 * t88;
t146 = qJ(2) - pkin(7);
t244 = pkin(8) - t146;
t115 = t244 * t149;
t116 = t244 * t151;
t85 = -t115 * t253 + t116 * t254;
t251 = t14 * t85;
t249 = t150 * t2;
t248 = t75 * mrSges(6,3);
t247 = t76 * mrSges(6,3);
t246 = t89 * Ifges(7,6);
t245 = t90 * Ifges(7,5);
t241 = Ifges(5,4) * t149;
t240 = Ifges(5,4) * t151;
t101 = Ifges(6,4) * t105;
t236 = t102 * Ifges(7,3);
t235 = t105 * Ifges(6,2);
t233 = t106 * Ifges(6,1);
t232 = t106 * Ifges(6,4);
t110 = t149 * t253 - t151 * t254;
t231 = t110 * t14;
t230 = t142 * Ifges(6,5);
t229 = t142 * Ifges(6,6);
t222 = Ifges(5,5) * qJD(4);
t221 = Ifges(5,6) * qJD(4);
t220 = t105 * t148;
t219 = t105 * t150;
t118 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t216;
t218 = t118 * t149;
t213 = qJD(3) * t122;
t207 = qJD(1) * qJD(3);
t133 = t149 * pkin(4) + t147;
t205 = t254 * pkin(4);
t204 = t253 * pkin(4);
t203 = pkin(4) * t216;
t193 = m(5) * t146 - mrSges(5,3);
t183 = m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3);
t180 = t76 * mrSges(6,1) - t75 * mrSges(6,2);
t80 = pkin(5) * t106 - pkin(9) * t105;
t179 = t148 * t2 + t150 * t3;
t83 = -t149 * t272 + t151 * t271;
t84 = t149 * t271 + t157;
t178 = t59 * t83 + t60 * t84;
t177 = mrSges(5,1) * t151 - mrSges(5,2) * t149;
t176 = -mrSges(7,1) * t150 + mrSges(7,2) * t148;
t173 = Ifges(7,1) * t148 + t238;
t171 = Ifges(7,2) * t150 + t239;
t169 = Ifges(7,5) * t148 + Ifges(7,6) * t150;
t167 = t148 * t18 + t150 * t17;
t164 = -t148 * t56 - t150 * t57;
t77 = pkin(5) * t111 + pkin(9) * t110 + t133;
t86 = -t115 * t254 - t116 * t253;
t39 = -t148 * t86 + t150 * t77;
t40 = t148 * t77 + t150 * t86;
t158 = t212 * t244 + t214;
t156 = -qJD(6) * t167 - t250;
t155 = t156 + t249;
t36 = t236 + t245 + t246;
t67 = t229 + t232 + t235;
t68 = t101 + t230 + t233;
t8 = t44 * Ifges(7,4) + t45 * Ifges(7,2) + t76 * Ifges(7,6);
t9 = t44 * Ifges(7,1) + t45 * Ifges(7,4) + t76 * Ifges(7,5);
t153 = -t13 * mrSges(6,2) - t38 * t219 / 0.2e1 + t37 * t220 / 0.2e1 + (Ifges(7,3) * t106 + t105 * t170) * t260 + (Ifges(7,5) * t106 + t105 * t174) * t262 + (Ifges(7,6) * t106 + t105 * t172) * t264 + t8 * t255 + t9 * t256 + t45 * t171 / 0.2e1 + t44 * t173 / 0.2e1 - t105 * t162 + t59 * t242 + mrSges(7,3) * t249 - Ifges(6,5) * t75 + t106 * t67 / 0.2e1 - t142 * (Ifges(6,5) * t105 - Ifges(6,6) * t106) / 0.2e1 - t18 * (-mrSges(7,2) * t106 - mrSges(7,3) * t220) - t17 * (mrSges(7,1) * t106 - mrSges(7,3) * t219) - t108 * (mrSges(6,1) * t106 + mrSges(6,2) * t105) + (t169 / 0.2e1 - Ifges(6,6)) * t76 - (-Ifges(6,2) * t106 + t101 + t68) * t105 / 0.2e1 - (Ifges(6,1) * t105 - t232 + t36) * t106 / 0.2e1 + (-mrSges(6,1) + t176) * t14 + t273 * qJD(6);
t137 = -t205 - pkin(5);
t117 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t217;
t112 = qJD(1) * (mrSges(5,1) * t149 + mrSges(5,2) * t151);
t104 = t222 + (t151 * Ifges(5,1) - t241) * qJD(1);
t103 = t221 + (-t149 * Ifges(5,2) + t240) * qJD(1);
t98 = qJD(2) * t149 - qJD(4) * t116;
t96 = t151 * t208 - t194;
t79 = -mrSges(6,1) * t105 + mrSges(6,2) * t106;
t72 = Ifges(7,3) * t76;
t66 = t80 + t203;
t62 = t100 * t254 - t201;
t61 = t100 * t253 + t202;
t41 = pkin(5) * t84 - pkin(9) * t83 + t123;
t35 = qJD(5) * t86 - t158 * t254 + t253 * t98;
t34 = -qJD(5) * t85 + t158 * t253 + t254 * t98;
t25 = t148 * t80 + t150 * t59;
t24 = -t148 * t59 + t150 * t80;
t23 = t148 * t66 + t150 * t62;
t22 = -t148 * t62 + t150 * t66;
t10 = -mrSges(7,1) * t45 + mrSges(7,2) * t44;
t5 = -qJD(6) * t40 - t148 * t34 + t150 * t41;
t4 = qJD(6) * t39 + t148 * t41 + t150 * t34;
t1 = [t39 * t15 + t40 * t16 + t4 * t56 + t5 * t57 + t85 * t10 + t34 * t91 + qJD(3) * t112 + t123 * t79 + t133 * t180 + (t17 * mrSges(7,1) - t18 * mrSges(7,2) + t246 / 0.2e1 + t245 / 0.2e1 + t236 / 0.2e1 + t36 / 0.2e1 - t67 / 0.2e1 - t235 / 0.2e1 - t232 / 0.2e1 + t108 * mrSges(6,1) - t229 / 0.2e1) * t84 + (t68 / 0.2e1 + t101 / 0.2e1 + t233 / 0.2e1 + t108 * mrSges(6,2) + t230 / 0.2e1 - t167 * mrSges(7,3) + t273) * t83 + t243 * t35 + 0.2e1 * (qJD(3) * mrSges(4,3) + qJD(2) * t183) * qJD(1) + (-t75 * t85 - t76 * t86 - t178) * mrSges(6,3) + m(5) * (qJD(2) * t265 + t147 * t207 + t213) + m(7) * (t17 * t5 + t18 * t4 + t2 * t40 + t3 * t39 + t35 * t53 + t251) + m(6) * (t108 * t123 + t114 * t133 + t13 * t86 + t34 * t60 - t35 * t59 + t251) + m(4) * (qJD(2) * t140 + t213 + (qJ(2) * qJD(2) + qJD(3) * t147) * qJD(1)) + (mrSges(5,2) * t207 + qJD(2) * t118 + t193 * t96 + (t146 * t117 - t103 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,4) * t216 - t221 / 0.2e1 + (t122 + t268) * mrSges(5,1)) * qJD(4)) * t151 + (-t13 * mrSges(6,3) + t72 / 0.2e1 + Ifges(6,4) * t75 + t114 * mrSges(6,1) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t76 + t267) * t111 + (t8 * t256 - t150 * t9 / 0.2e1 + Ifges(6,1) * t75 - t114 * mrSges(6,2) - t44 * t174 / 0.2e1 - t45 * t172 / 0.2e1 + (-mrSges(6,3) - t175) * t14 + t179 * mrSges(7,3) + (mrSges(7,3) * t269 + t169 * t259 + t171 * t263 + t173 * t261 + t176 * t53 + t255 * t37 + t256 * t38) * qJD(6) + (Ifges(6,4) - t170 / 0.2e1) * t76) * t110 + (mrSges(5,1) * t207 + qJD(2) * t117 + t193 * t97 + (-t146 * t118 - t104 / 0.2e1 - t122 * mrSges(5,2) - t222 / 0.2e1 + (-t147 * mrSges(5,2) + 0.3e1 / 0.2e1 * t241 + (-0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2)) * t151) * qJD(1)) * qJD(4)) * t149; -t148 * t16 - t150 * t15 + t243 * t106 + t165 * qJD(6) - t183 * t152 + t163 * t105 + (-t149 * t117 - t151 * t118 - t177 * qJD(4) + (-qJD(3) - t265) * m(5) + (-qJD(3) - t140) * m(4)) * qJD(1) - t180 + (-t102 * t269 + t106 * t53 - t179) * m(7) + (t105 * t60 - t106 * t59 - t114) * m(6); -t152 * mrSges(4,3) - t243 * t83 + (t10 - t248) * t110 + (t117 * t151 - t218) * qJD(4) + t163 * t84 + (qJD(6) * t164 + t168 - t247) * t111 + m(6) * (t111 * t13 + t178 + t231) + m(7) * (t111 * t155 + t269 * t84 - t53 * t83 + t231) + m(5) * (t149 * t97 + t151 * t96) + (-m(6) * t108 - m(5) * t122 - m(7) * t167 - t79 - t112 + (qJD(2) - t122) * m(4) + t164) * qJD(1); t153 - t204 * t247 + t205 * t248 + (-t151 * (-Ifges(5,1) * t149 - t240) / 0.2e1 + t149 * (-Ifges(5,2) * t151 - t241) / 0.2e1) * t152 - t122 * t177 * qJD(1) + t266 * mrSges(7,3) + (m(7) * t155 + t274) * (t204 + pkin(9)) - t117 * t113 + t103 * t216 / 0.2e1 + t104 * t217 / 0.2e1 - t79 * t203 - Ifges(5,6) * t186 / 0.2e1 + (t14 * t137 - t17 * t22 - t18 * t23 - t53 * t61) * m(7) + (-t108 * t203 + t59 * t61 - t60 * t62) * m(6) + (t163 * t187 + (t253 * t53 + t254 * t269) * qJD(5) * m(7) + t243 * t189 + m(6) * (-t254 * t14 + t253 * t13 + (-t253 * t59 + t254 * t60) * qJD(5))) * pkin(4) - t243 * t61 + t121 * t218 + t60 * t234 - t206 * Ifges(5,5) * t149 / 0.2e1 - t23 * t56 - t22 * t57 - t62 * t91 + t96 * mrSges(5,1) - t97 * mrSges(5,2) + t137 * t10; t153 + t274 * pkin(9) - pkin(5) * t10 + (t234 - t243) * t60 + t156 * mrSges(7,3) - t25 * t56 - t24 * t57 - t59 * t91 + ((t249 + t266) * pkin(9) - t17 * t24 - t18 * t25 - t53 * t60 - t14 * pkin(5)) * m(7); t72 - t53 * (mrSges(7,1) * t90 + mrSges(7,2) * t89) + (Ifges(7,1) * t89 - t252) * t262 + t37 * t261 + (Ifges(7,5) * t89 - Ifges(7,6) * t90) * t260 - t17 * t56 + t18 * t57 + (t17 * t89 + t18 * t90) * mrSges(7,3) + (-Ifges(7,2) * t90 + t38 + t87) * t264 + t267;];
tauc  = t1(:);
