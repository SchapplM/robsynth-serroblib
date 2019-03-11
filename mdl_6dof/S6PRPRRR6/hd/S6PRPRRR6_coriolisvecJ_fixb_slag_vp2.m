% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:23
% EndTime: 2019-03-08 20:45:37
% DurationCPUTime: 7.71s
% Computational Cost: add. (5345->506), mult. (12187->720), div. (0->0), fcn. (8168->10), ass. (0->249)
t168 = sin(qJ(5));
t169 = sin(qJ(4));
t174 = cos(qJ(2));
t165 = sin(pkin(6));
t237 = qJD(1) * t165;
t170 = sin(qJ(2));
t172 = cos(qJ(5));
t238 = t170 * t172;
t106 = (t168 * t174 + t169 * t238) * t237;
t173 = cos(qJ(4));
t201 = pkin(4) * t173 + pkin(9) * t169;
t135 = qJD(4) * t201 + qJD(3);
t149 = pkin(4) * t169 - pkin(9) * t173 + qJ(3);
t175 = -pkin(2) - pkin(8);
t228 = qJD(4) * t175;
t209 = t173 * t228;
t239 = t169 * t175;
t216 = t168 * t239;
t226 = qJD(5) * t172;
t58 = -qJD(5) * t216 + t168 * t135 + t149 * t226 + t172 * t209;
t334 = t58 - t106;
t240 = t168 * t170;
t105 = (-t169 * t240 + t172 * t174) * t237;
t119 = t172 * t135;
t157 = t172 * t239;
t205 = -t168 * t175 + pkin(5);
t222 = pkin(10) * t169 * t172;
t270 = pkin(10) * t173;
t333 = -t105 + t119 + (-t157 + (-t149 + t270) * t168) * qJD(5) + (t173 * t205 + t222) * qJD(4);
t225 = qJD(5) * t173;
t231 = qJD(4) * t169;
t180 = t168 * t231 - t172 * t225;
t332 = -pkin(10) * t180 - t334;
t164 = t169 * qJD(2);
t291 = -pkin(10) - pkin(9);
t215 = qJD(5) * t291;
t144 = t201 * qJD(2);
t213 = t174 * t237;
t192 = qJD(3) - t213;
t123 = qJD(2) * t175 + t192;
t166 = cos(pkin(6));
t236 = qJD(1) * t169;
t94 = t123 * t173 - t166 * t236;
t61 = t168 * t144 + t172 * t94;
t331 = -t61 + (-pkin(10) * t164 + t215) * t168;
t60 = t172 * t144 - t168 * t94;
t330 = t172 * t215 - (pkin(5) * t173 + t222) * qJD(2) - t60;
t219 = mrSges(5,3) * t164;
t153 = -qJD(4) * mrSges(5,2) - t219;
t241 = t166 * t173;
t212 = qJD(1) * t241;
t95 = t123 * t169 + t212;
t329 = m(5) * t95 + t153;
t167 = sin(qJ(6));
t171 = cos(qJ(6));
t134 = t172 * t149;
t75 = t169 * t205 - t172 * t270 + t134;
t111 = t168 * t149 + t157;
t88 = -t168 * t270 + t111;
t36 = -t167 * t88 + t171 * t75;
t328 = qJD(6) * t36 + t167 * t333 - t332 * t171;
t37 = t167 * t75 + t171 * t88;
t327 = -qJD(6) * t37 + t332 * t167 + t171 * t333;
t230 = qJD(4) * t172;
t234 = qJD(2) * t173;
t139 = -t168 * t234 + t230;
t162 = t164 + qJD(5);
t214 = t170 * t237;
t109 = qJD(2) * t149 + t214;
t84 = qJD(4) * pkin(9) + t95;
t43 = t172 * t109 - t168 * t84;
t44 = t109 * t168 + t172 * t84;
t190 = t168 * t44 + t172 * t43;
t263 = Ifges(6,4) * t172;
t196 = -Ifges(6,2) * t168 + t263;
t264 = Ifges(6,4) * t168;
t198 = Ifges(6,1) * t172 - t264;
t199 = mrSges(6,1) * t168 + mrSges(6,2) * t172;
t261 = Ifges(6,6) * t168;
t262 = Ifges(6,5) * t172;
t273 = t172 / 0.2e1;
t274 = -t168 / 0.2e1;
t232 = qJD(4) * t168;
t140 = t172 * t234 + t232;
t279 = t140 / 0.2e1;
t258 = t140 * Ifges(6,4);
t69 = t139 * Ifges(6,2) + t162 * Ifges(6,6) + t258;
t132 = Ifges(6,4) * t139;
t70 = Ifges(6,1) * t140 + Ifges(6,5) * t162 + t132;
t83 = -qJD(4) * pkin(4) - t94;
t326 = -t190 * mrSges(6,3) + (-t261 + t262) * t162 / 0.2e1 + t196 * t139 / 0.2e1 + t198 * t279 + t83 * t199 + t273 * t70 + t274 * t69;
t32 = pkin(10) * t139 + t44;
t254 = t167 * t32;
t31 = -pkin(10) * t140 + t43;
t27 = pkin(5) * t162 + t31;
t12 = t171 * t27 - t254;
t251 = t171 * t32;
t13 = t167 * t27 + t251;
t224 = qJD(2) * qJD(4);
t206 = t173 * t224;
t160 = Ifges(7,3) * t206;
t204 = t171 * t139 - t140 * t167;
t79 = t139 * t167 + t140 * t171;
t272 = Ifges(7,4) * t79;
t158 = qJD(6) + t162;
t278 = -t158 / 0.2e1;
t293 = -t79 / 0.2e1;
t295 = -t204 / 0.2e1;
t72 = Ifges(7,4) * t204;
t30 = Ifges(7,1) * t79 + Ifges(7,5) * t158 + t72;
t65 = -pkin(5) * t139 + t83;
t325 = t160 + (Ifges(7,5) * t204 - Ifges(7,6) * t79) * t278 + (t12 * t204 + t13 * t79) * mrSges(7,3) + (-Ifges(7,2) * t79 + t30 + t72) * t295 - t65 * (mrSges(7,1) * t79 + mrSges(7,2) * t204) + (Ifges(7,1) * t204 - t272) * t293;
t155 = t291 * t168;
t156 = t291 * t172;
t101 = t155 * t167 - t156 * t171;
t324 = -qJD(6) * t101 - t167 * t331 + t330 * t171;
t100 = t155 * t171 + t156 * t167;
t323 = qJD(6) * t100 + t330 * t167 + t171 * t331;
t187 = t167 * t168 - t171 * t172;
t306 = qJD(5) + qJD(6);
t322 = t306 * t187;
t147 = qJD(2) * qJ(3) + t214;
t247 = Ifges(5,5) * qJD(4);
t321 = t147 * mrSges(5,2) + t247 / 0.2e1 + (Ifges(5,1) * t173 - Ifges(5,4) * t169) * qJD(2) / 0.2e1 - t94 * mrSges(5,3) + t326;
t223 = qJD(4) * qJD(5);
t102 = t172 * t223 + (-t168 * t225 - t169 * t230) * qJD(2);
t235 = qJD(2) * t170;
t211 = t165 * t235;
t229 = qJD(4) * t173;
t66 = t123 * t229 + (-qJD(4) * t166 + t211) * t236;
t96 = (t135 + t213) * qJD(2);
t18 = -qJD(5) * t44 - t168 * t66 + t172 * t96;
t11 = pkin(5) * t206 - pkin(10) * t102 + t18;
t103 = qJD(2) * t180 - t168 * t223;
t227 = qJD(5) * t168;
t17 = t109 * t226 + t168 * t96 + t172 * t66 - t227 * t84;
t14 = pkin(10) * t103 + t17;
t2 = qJD(6) * t12 + t11 * t167 + t14 * t171;
t25 = qJD(6) * t204 + t102 * t171 + t103 * t167;
t26 = -qJD(6) * t79 - t102 * t167 + t103 * t171;
t3 = -qJD(6) * t13 + t11 * t171 - t14 * t167;
t320 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t25 + Ifges(7,6) * t26;
t29 = Ifges(7,2) * t204 + Ifges(7,6) * t158 + t272;
t317 = t29 / 0.2e1;
t207 = -Ifges(5,6) * qJD(4) / 0.2e1;
t288 = m(7) * t65;
t38 = -mrSges(7,1) * t204 + mrSges(7,2) * t79;
t311 = t38 + t288;
t117 = t187 * t173;
t142 = t167 * t172 + t168 * t171;
t127 = t142 * qJD(2);
t86 = t306 * t142;
t310 = -qJD(4) * t117 - t169 * t86 - t127;
t115 = t142 * t173;
t309 = t187 * qJD(2) - qJD(4) * t115 + t169 * t322;
t308 = t187 * t169;
t191 = -t168 * t18 + t17 * t172;
t307 = -m(5) * t94 + m(6) * t83;
t107 = -mrSges(6,2) * t162 + mrSges(6,3) * t139;
t108 = mrSges(6,1) * t162 - mrSges(6,3) * t140;
t188 = -t168 * t107 - t172 * t108;
t304 = -m(6) * t190 + t188;
t303 = t18 * mrSges(6,1) - t17 * mrSges(6,2) + Ifges(6,5) * t102 + Ifges(6,6) * t103 + t320;
t300 = -m(5) / 0.2e1;
t299 = m(5) / 0.2e1;
t298 = m(6) / 0.2e1;
t297 = t25 / 0.2e1;
t296 = t26 / 0.2e1;
t294 = t204 / 0.2e1;
t292 = t79 / 0.2e1;
t286 = t102 / 0.2e1;
t285 = t103 / 0.2e1;
t284 = -t115 / 0.2e1;
t283 = -t117 / 0.2e1;
t282 = -t139 / 0.2e1;
t280 = -t140 / 0.2e1;
t277 = t158 / 0.2e1;
t276 = -t162 / 0.2e1;
t271 = pkin(5) * t168;
t10 = -mrSges(7,1) * t26 + mrSges(7,2) * t25;
t57 = -mrSges(6,1) * t103 + mrSges(6,2) * t102;
t266 = -t10 - t57;
t265 = Ifges(5,4) * t173;
t242 = t165 * t174;
t121 = t166 * t169 + t173 * t242;
t203 = t173 * t211;
t245 = qJD(4) * t95;
t67 = -qJD(1) * t203 + t245;
t260 = t121 * t67;
t112 = t169 * t127;
t250 = t112 + t86;
t113 = qJD(2) * t308;
t249 = t113 + t322;
t248 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t139 + mrSges(6,2) * t140 + mrSges(5,3) * t234;
t243 = t147 * t174;
t233 = qJD(2) * t174;
t221 = t288 / 0.2e1;
t220 = t38 + t248;
t218 = Ifges(5,4) * t164;
t217 = t169 * t242;
t210 = t165 * t233;
t208 = -t247 / 0.2e1;
t200 = mrSges(6,1) * t172 - mrSges(6,2) * t168;
t197 = Ifges(6,1) * t168 + t263;
t195 = Ifges(6,2) * t172 + t264;
t193 = Ifges(6,5) * t168 + Ifges(6,6) * t172;
t122 = -t217 + t241;
t91 = -t122 * t168 + t165 * t238;
t92 = t122 * t172 + t165 * t240;
t41 = -t167 * t92 + t171 * t91;
t42 = t167 * t91 + t171 * t92;
t73 = mrSges(6,1) * t206 - mrSges(6,3) * t102;
t74 = -mrSges(6,2) * t206 + mrSges(6,3) * t103;
t189 = -t168 * t73 + t172 * t74;
t136 = (qJD(3) + t213) * qJD(2);
t183 = t136 * t170 + t147 * t233;
t178 = t12 * mrSges(7,1) + t147 * mrSges(5,1) + t43 * mrSges(6,1) + t207 - (-t169 * Ifges(5,2) + t265) * qJD(2) / 0.2e1 + t158 * Ifges(7,3) + t79 * Ifges(7,5) + t204 * Ifges(7,6) + t162 * Ifges(6,3) + t140 * Ifges(6,5) + t139 * Ifges(6,6) - t13 * mrSges(7,2) - t44 * mrSges(6,2) - t95 * mrSges(5,3);
t176 = qJD(2) ^ 2;
t163 = -pkin(5) * t172 - pkin(4);
t161 = Ifges(6,3) * t206;
t143 = (mrSges(5,1) * t169 + mrSges(5,2) * t173) * qJD(2);
t138 = -qJD(2) * pkin(2) + t192;
t137 = (-t175 + t271) * t173;
t130 = (mrSges(5,1) * t173 - mrSges(5,2) * t169) * t224;
t114 = t142 * t169;
t110 = t134 - t216;
t104 = -pkin(5) * t180 + t169 * t228;
t90 = -qJD(4) * t217 + t166 * t229 - t203;
t89 = -qJD(4) * t121 + t169 * t211;
t71 = t212 + (-qJD(2) * t271 + t123) * t169;
t63 = mrSges(7,1) * t158 - mrSges(7,3) * t79;
t62 = -mrSges(7,2) * t158 + mrSges(7,3) * t204;
t59 = -qJD(5) * t111 - t168 * t209 + t119;
t51 = t102 * Ifges(6,1) + t103 * Ifges(6,4) + Ifges(6,5) * t206;
t50 = t102 * Ifges(6,4) + t103 * Ifges(6,2) + Ifges(6,6) * t206;
t48 = t142 * t231 + t173 * t322;
t46 = qJD(4) * t308 - t173 * t86;
t40 = -pkin(5) * t103 + t67;
t34 = qJD(5) * t91 + t168 * t210 + t172 * t89;
t33 = -qJD(5) * t92 - t168 * t89 + t172 * t210;
t22 = -mrSges(7,2) * t206 + mrSges(7,3) * t26;
t21 = mrSges(7,1) * t206 - mrSges(7,3) * t25;
t16 = t171 * t31 - t254;
t15 = -t167 * t31 - t251;
t9 = t25 * Ifges(7,1) + t26 * Ifges(7,4) + Ifges(7,5) * t206;
t8 = t25 * Ifges(7,4) + t26 * Ifges(7,2) + Ifges(7,6) * t206;
t5 = -qJD(6) * t42 - t167 * t34 + t171 * t33;
t4 = qJD(6) * t41 + t167 * t33 + t171 * t34;
t1 = [-t122 * mrSges(5,3) * t206 + t34 * t107 + t33 * t108 + t89 * t153 + t41 * t21 + t42 * t22 + t4 * t62 + t5 * t63 + t91 * t73 + t92 * t74 + t220 * t90 + (-qJD(4) * t219 - t266) * t121 + m(7) * (t12 * t5 + t121 * t40 + t13 * t4 + t2 * t42 + t3 * t41 + t65 * t90) + m(6) * (t17 * t92 + t18 * t91 + t33 * t43 + t34 * t44 + t83 * t90 + t260) + m(5) * (t122 * t66 + t89 * t95 - t90 * t94 + t260) + (t143 * t233 + t170 * t130 + ((-mrSges(3,2) + mrSges(4,3)) * t174 + (-mrSges(3,1) + mrSges(4,2)) * t170) * t176 + 0.2e1 * t183 * t299 + (t138 * t235 - t214 * t233 + t183) * m(4)) * t165; (-Ifges(7,4) * t117 - Ifges(7,2) * t115) * t296 + (-Ifges(7,1) * t117 - Ifges(7,4) * t115) * t297 + t40 * (mrSges(7,1) * t115 - mrSges(7,2) * t117) + (-t115 * t2 + t117 * t3 - t12 * t46 + t13 * t48) * mrSges(7,3) + 0.2e1 * (t243 * t300 - (pkin(2) * t235 + t138 * t170 + t243) * m(4) / 0.2e1) * t237 + t334 * t107 - m(6) * (t105 * t43 + t106 * t44) + m(6) * (t110 * t18 + t111 * t17 + t43 * t59 + t44 * t58) + (t136 * mrSges(5,2) + t198 * t286 + t196 * t285 + t50 * t274 + t51 * t273 + (-t168 * t17 - t172 * t18) * mrSges(6,3) + (t83 * t200 + t195 * t282 + t197 * t280 + t193 * t276 - t172 * t69 / 0.2e1 + t70 * t274 + (t43 * t168 - t44 * t172) * mrSges(6,3)) * qJD(5) + (t178 + ((t262 / 0.2e1 - t261 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,4)) * t173 + Ifges(7,5) * t283 + Ifges(7,6) * t284 + (Ifges(6,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(5,1) + Ifges(7,3) / 0.2e1) * t169) * qJD(2) + t207) * qJD(4) + (t329 * qJD(4) - t57) * t175 + (t220 + 0.2e1 * t221 + t307) * t214 + (mrSges(5,3) + t199 + (0.2e1 * t300 - m(6)) * t175) * t67) * t173 + (t160 / 0.2e1 + t161 / 0.2e1 + t136 * mrSges(5,1) + (0.3e1 / 0.2e1 * t218 + t208 + (t248 + t307) * t175 - t321) * qJD(4) + (m(5) * t175 - mrSges(5,3)) * t66 + t303 - t329 * t214) * t169 + (m(4) + m(5)) * (qJ(3) * t136 + qJD(3) * t147) + (t59 - t105) * t108 + (qJD(2) * t192 + t136) * mrSges(4,3) + t192 * t143 + (Ifges(7,5) * t46 + Ifges(7,6) * t48) * t277 + t9 * t283 + t8 * t284 + (Ifges(7,1) * t46 + Ifges(7,4) * t48) * t292 + (Ifges(7,4) * t46 + Ifges(7,2) * t48) * t294 + t36 * t21 + t37 * t22 + t46 * t30 / 0.2e1 + t65 * (-mrSges(7,1) * t48 + mrSges(7,2) * t46) + t104 * t38 + t110 * t73 + t111 * t74 + qJ(3) * t130 + t137 * t10 + t327 * t63 + t328 * t62 + (t104 * t65 + t327 * t12 + t13 * t328 + t137 * t40 + t2 * t37 + t3 * t36) * m(7) + t48 * t317; -t176 * mrSges(4,3) - t114 * t21 - t308 * t22 + t309 * t63 + t310 * t62 + ((t107 * t172 - t108 * t168 + t153) * qJD(4) + t266) * t173 + (qJD(4) * t220 + qJD(5) * t188 + t189) * t169 + 0.2e1 * ((-t67 + t245) * t299 + (t230 * t44 - t232 * t43 - t67) * t298) * t173 + 0.2e1 * ((-qJD(4) * t94 + t66) * t299 + qJD(4) * t221 + (qJD(4) * t83 - t226 * t43 - t227 * t44 + t191) * t298) * t169 + (-m(5) * t147 - t143 + (-t147 + t214) * m(4) + t304) * qJD(2) + (-t114 * t3 + t12 * t309 + t13 * t310 - t40 * t173 - t2 * t308) * m(7); (-pkin(4) * t67 + pkin(9) * t191 - t43 * t60 - t44 * t61 - t83 * t95) * m(6) + (pkin(9) * t304 + t271 * t311 + t326) * qJD(5) + ((t208 - t218 / 0.2e1 + t321) * t169 + (-t178 + (t265 / 0.2e1 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t169) * qJD(2) + t207 + (Ifges(7,5) * t142 - Ifges(7,6) * t187 + t193) * qJD(4) / 0.2e1) * t173) * qJD(2) + (-t322 / 0.2e1 - t113 / 0.2e1) * t30 + (-Ifges(7,5) * t322 - Ifges(7,6) * t86) * t277 + (-Ifges(7,1) * t322 - Ifges(7,4) * t86) * t292 + (-Ifges(7,4) * t322 - Ifges(7,2) * t86) * t294 + (mrSges(7,1) * t250 - mrSges(7,2) * t249) * t65 - t248 * t95 + t323 * t62 + t324 * t63 + (t100 * t3 + t101 * t2 + t12 * t324 + t13 * t323 + t163 * t40 - t65 * t71) * m(7) + t189 * pkin(9) + t191 * mrSges(6,3) + (-mrSges(5,1) - t200) * t67 + t50 * t273 + (Ifges(7,5) * t113 + Ifges(7,6) * t112) * t278 + t195 * t285 + t197 * t286 + (Ifges(7,1) * t113 + Ifges(7,4) * t112) * t293 + (Ifges(7,4) * t113 + Ifges(7,2) * t112) * t295 - pkin(4) * t57 - t66 * mrSges(5,2) - t71 * t38 + t100 * t21 + t101 * t22 - t61 * t107 - t60 * t108 + t142 * t9 / 0.2e1 - t94 * t153 + t163 * t10 + t168 * t51 / 0.2e1 + (t12 * t249 - t13 * t250 - t142 * t3 - t187 * t2) * mrSges(7,3) + (Ifges(7,4) * t142 - Ifges(7,2) * t187) * t296 + (Ifges(7,1) * t142 - Ifges(7,4) * t187) * t297 + t40 * (mrSges(7,1) * t187 + mrSges(7,2) * t142) - t187 * t8 / 0.2e1 + (-t86 / 0.2e1 - t112 / 0.2e1) * t29; (t139 * t43 + t140 * t44) * mrSges(6,3) + (-Ifges(6,2) * t140 + t132 + t70) * t282 + (t167 * t22 + t171 * t21 + m(7) * (t167 * t2 + t171 * t3) - t311 * t140 + (-t167 * t63 + t171 * t62 + m(7) * (-t12 * t167 + t13 * t171)) * qJD(6)) * pkin(5) - m(7) * (t12 * t15 + t13 * t16) + t303 + t161 + (Ifges(6,5) * t139 - Ifges(6,6) * t140) * t276 + t69 * t279 + (Ifges(6,1) * t139 - t258) * t280 - t16 * t62 - t15 * t63 - t43 * t107 + t44 * t108 - t83 * (mrSges(6,1) * t140 + mrSges(6,2) * t139) + t79 * t317 + t325; -t12 * t62 + t13 * t63 + t29 * t292 + t320 + t325;];
tauc  = t1(:);
