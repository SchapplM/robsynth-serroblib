% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRRR6
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

function tauc = S6RPPRRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:50:33
% EndTime: 2018-11-23 15:50:38
% DurationCPUTime: 4.76s
% Computational Cost: add. (4638->461), mult. (9549->642), div. (0->0), fcn. (5528->6), ass. (0->217)
t154 = sin(qJ(5));
t155 = sin(qJ(4));
t147 = t155 * qJD(1);
t195 = t154 * t147;
t260 = -pkin(9) - pkin(8);
t196 = qJD(5) * t260;
t158 = cos(qJ(4));
t187 = pkin(4) * t158 + pkin(8) * t155;
t122 = t187 * qJD(1);
t146 = qJD(1) * qJ(2) + qJD(3);
t134 = -qJD(1) * pkin(7) + t146;
t157 = cos(qJ(5));
t213 = t157 * t158;
t67 = t154 * t122 + t134 * t213;
t298 = -pkin(9) * t195 + t154 * t196 - t67;
t215 = t155 * t157;
t171 = pkin(5) * t158 + pkin(9) * t215;
t216 = t154 * t158;
t66 = t157 * t122 - t134 * t216;
t297 = -qJD(1) * t171 + t157 * t196 - t66;
t108 = -qJD(4) * pkin(4) - t134 * t158;
t207 = qJD(4) * t157;
t212 = qJD(1) * t158;
t117 = -t154 * t212 + t207;
t142 = t147 + qJD(5);
t125 = t155 * t134;
t107 = qJD(4) * pkin(8) + t125;
t220 = t107 * t154;
t152 = pkin(1) + qJ(3);
t126 = pkin(4) * t155 - pkin(8) * t158 + t152;
t93 = qJD(1) * t126 - qJD(2);
t51 = t157 * t93 - t220;
t52 = t107 * t157 + t154 * t93;
t176 = t154 * t52 + t157 * t51;
t237 = Ifges(6,4) * t157;
t182 = -Ifges(6,2) * t154 + t237;
t238 = Ifges(6,4) * t154;
t184 = Ifges(6,1) * t157 - t238;
t185 = mrSges(6,1) * t154 + mrSges(6,2) * t157;
t235 = Ifges(6,6) * t154;
t236 = Ifges(6,5) * t157;
t248 = t157 / 0.2e1;
t249 = -t154 / 0.2e1;
t118 = qJD(4) * t154 + t157 * t212;
t254 = t118 / 0.2e1;
t239 = Ifges(6,4) * t118;
t55 = Ifges(6,2) * t117 + Ifges(6,6) * t142 + t239;
t113 = Ifges(6,4) * t117;
t56 = Ifges(6,1) * t118 + Ifges(6,5) * t142 + t113;
t296 = -t176 * mrSges(6,3) + t108 * t185 + (-t235 + t236) * t142 / 0.2e1 + t182 * t117 / 0.2e1 + t184 * t254 + t248 * t56 + t249 * t55;
t153 = sin(qJ(6));
t156 = cos(qJ(6));
t45 = pkin(9) * t117 + t52;
t227 = t156 * t45;
t44 = -pkin(9) * t118 + t51;
t33 = pkin(5) * t142 + t44;
t10 = t153 * t33 + t227;
t206 = qJD(4) * t158;
t190 = qJD(1) * t206;
t140 = Ifges(7,3) * t190;
t189 = t156 * t117 - t118 * t153;
t64 = t117 * t153 + t118 * t156;
t247 = Ifges(7,4) * t64;
t199 = qJD(5) + qJD(6);
t136 = t147 + t199;
t253 = -t136 / 0.2e1;
t266 = -t64 / 0.2e1;
t268 = -t189 / 0.2e1;
t58 = Ifges(7,4) * t189;
t28 = Ifges(7,1) * t64 + Ifges(7,5) * t136 + t58;
t70 = -pkin(5) * t117 + t108;
t229 = t153 * t45;
t9 = t156 * t33 - t229;
t295 = t140 + (Ifges(7,5) * t189 - Ifges(7,6) * t64) * t253 + (t10 * t64 + t189 * t9) * mrSges(7,3) + (-Ifges(7,2) * t64 + t28 + t58) * t268 - t70 * (mrSges(7,1) * t64 + mrSges(7,2) * t189) + (Ifges(7,1) * t189 - t247) * t266;
t131 = t260 * t154;
t132 = t260 * t157;
t75 = t131 * t156 + t132 * t153;
t294 = qJD(6) * t75 + t297 * t153 + t298 * t156;
t76 = t131 * t153 - t132 * t156;
t293 = -qJD(6) * t76 - t298 * t153 + t297 * t156;
t203 = qJD(6) * t156;
t205 = qJD(5) * t157;
t218 = t153 * t154;
t68 = -t156 * t205 - t157 * t203 + t199 * t218;
t214 = t156 * t157;
t89 = t147 * t214 - t153 * t195;
t242 = t68 - t89;
t120 = t153 * t157 + t154 * t156;
t292 = t199 * t120;
t106 = t120 * qJD(1);
t88 = t155 * t106;
t241 = t292 + t88;
t116 = qJD(4) * t187 + qJD(3);
t151 = qJ(2) - pkin(7);
t291 = -qJD(5) * t151 * t155 + t116;
t224 = Ifges(5,5) * qJD(4);
t290 = t224 / 0.2e1 + (t158 * Ifges(5,1) - Ifges(5,4) * t155) * qJD(1) / 0.2e1 + t296;
t289 = qJD(2) * t155 + t151 * t206;
t200 = qJD(4) * qJD(5);
t204 = qJD(5) * t158;
t77 = t157 * t200 + (-t154 * t204 - t155 * t207) * qJD(1);
t208 = qJD(4) * t155;
t164 = t154 * t208 - t157 * t204;
t78 = qJD(1) * t164 - t154 * t200;
t19 = qJD(6) * t189 + t153 * t78 + t156 * t77;
t202 = qJD(1) * qJD(2);
t92 = t134 * t206 + t155 * t202;
t94 = t116 * qJD(1);
t22 = -qJD(5) * t52 - t154 * t92 + t157 * t94;
t11 = pkin(5) * t190 - pkin(9) * t77 + t22;
t21 = -qJD(5) * t220 + t154 * t94 + t157 * t92 + t93 * t205;
t14 = pkin(9) * t78 + t21;
t280 = qJD(6) * t9;
t2 = t11 * t153 + t14 * t156 + t280;
t20 = -qJD(6) * t64 - t153 * t77 + t156 * t78;
t3 = -qJD(6) * t10 + t11 * t156 - t14 * t153;
t288 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t19 + Ifges(7,6) * t20;
t27 = Ifges(7,2) * t189 + Ifges(7,6) * t136 + t247;
t286 = t27 / 0.2e1;
t191 = -Ifges(5,6) * qJD(4) / 0.2e1;
t285 = m(6) * t108;
t259 = m(7) * t70;
t32 = -mrSges(7,1) * t189 + mrSges(7,2) * t64;
t279 = t32 + t259;
t119 = -t214 + t218;
t98 = t119 * t158;
t278 = -qJD(4) * t98 - t155 * t292 - t106;
t162 = t199 * t119;
t96 = t120 * t158;
t277 = t119 * qJD(1) - qJD(4) * t96 + t155 * t162;
t178 = -t154 * t22 + t157 * t21;
t276 = qJD(1) * t152;
t82 = -mrSges(6,2) * t142 + mrSges(6,3) * t117;
t83 = mrSges(6,1) * t142 - mrSges(6,3) * t118;
t172 = -t154 * t82 - t157 * t83;
t274 = -m(6) * t176 + t172;
t273 = (t155 ^ 2 + t158 ^ 2) * t134;
t272 = t22 * mrSges(6,1) - t21 * mrSges(6,2) + Ifges(6,5) * t77 + Ifges(6,6) * t78 + t288;
t270 = t19 / 0.2e1;
t269 = t20 / 0.2e1;
t267 = t189 / 0.2e1;
t265 = t64 / 0.2e1;
t264 = t77 / 0.2e1;
t263 = t78 / 0.2e1;
t262 = -t96 / 0.2e1;
t261 = -t98 / 0.2e1;
t257 = -t117 / 0.2e1;
t255 = -t118 / 0.2e1;
t252 = t136 / 0.2e1;
t251 = -t142 / 0.2e1;
t246 = pkin(5) * t154;
t240 = Ifges(5,4) * t158;
t91 = t134 * t208 - t158 * t202;
t230 = t151 * t91;
t222 = qJD(4) * mrSges(5,1);
t225 = -mrSges(6,1) * t117 + mrSges(6,2) * t118 + mrSges(5,3) * t212 - t222;
t221 = qJD(4) * mrSges(5,2);
t217 = t154 * t155;
t81 = t154 * t126 + t151 * t215;
t135 = -qJD(2) + t276;
t209 = qJD(3) * t135;
t201 = qJD(1) * qJD(3);
t198 = t32 + t225;
t197 = Ifges(5,4) * t147;
t192 = -t224 / 0.2e1;
t188 = m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3);
t186 = mrSges(6,1) * t157 - mrSges(6,2) * t154;
t183 = Ifges(6,1) * t154 + t237;
t181 = Ifges(6,2) * t157 + t238;
t179 = Ifges(6,5) * t154 + Ifges(6,6) * t157;
t112 = t157 * t126;
t57 = -pkin(9) * t213 + t112 + (-t151 * t154 + pkin(5)) * t155;
t65 = -pkin(9) * t216 + t81;
t30 = -t153 * t65 + t156 * t57;
t31 = t153 * t57 + t156 * t65;
t177 = -t154 * t21 - t157 * t22;
t175 = t154 * t51 - t157 * t52;
t59 = mrSges(6,1) * t190 - mrSges(6,3) * t77;
t60 = -mrSges(6,2) * t190 + mrSges(6,3) * t78;
t174 = -t154 * t59 + t157 * t60;
t173 = t154 * t83 - t157 * t82;
t170 = t291 * t157;
t129 = -mrSges(5,3) * t147 - t221;
t169 = t129 - t173;
t34 = t126 * t205 + t154 * t291 + t157 * t289;
t161 = t135 * mrSges(5,1) + t51 * mrSges(6,1) + t9 * mrSges(7,1) + t191 - (-t155 * Ifges(5,2) + t240) * qJD(1) / 0.2e1 + t136 * Ifges(7,3) + t64 * Ifges(7,5) + t189 * Ifges(7,6) + t142 * Ifges(6,3) + t118 * Ifges(6,5) + t117 * Ifges(6,6) - t10 * mrSges(7,2) - t52 * mrSges(6,2);
t159 = qJD(1) ^ 2;
t145 = -pkin(5) * t157 - pkin(4);
t141 = Ifges(6,3) * t190;
t121 = (mrSges(5,1) * t155 + mrSges(5,2) * t158) * qJD(1);
t115 = (-t151 + t246) * t158;
t97 = t119 * t155;
t95 = t120 * t155;
t90 = -pkin(5) * t195 + t125;
t80 = -t151 * t217 + t112;
t72 = -pkin(5) * t164 - qJD(2) * t158 + t151 * t208;
t49 = mrSges(7,1) * t136 - mrSges(7,3) * t64;
t48 = -mrSges(7,2) * t136 + mrSges(7,3) * t189;
t47 = -pkin(5) * t78 + t91;
t46 = -mrSges(6,1) * t78 + mrSges(6,2) * t77;
t41 = t77 * Ifges(6,1) + t78 * Ifges(6,4) + Ifges(6,5) * t190;
t40 = t77 * Ifges(6,4) + t78 * Ifges(6,2) + Ifges(6,6) * t190;
t39 = t120 * t208 + t158 * t162;
t37 = t119 * t208 - t158 * t292;
t35 = (-qJD(5) * t126 - t289) * t154 + t170;
t29 = pkin(9) * t164 + t34;
t25 = t171 * qJD(4) + ((pkin(9) * t158 - t126) * qJD(5) - t289) * t154 + t170;
t16 = -mrSges(7,2) * t190 + mrSges(7,3) * t20;
t15 = mrSges(7,1) * t190 - mrSges(7,3) * t19;
t13 = t156 * t44 - t229;
t12 = -t153 * t44 - t227;
t8 = -mrSges(7,1) * t20 + mrSges(7,2) * t19;
t7 = t19 * Ifges(7,1) + t20 * Ifges(7,4) + Ifges(7,5) * t190;
t6 = t19 * Ifges(7,4) + t20 * Ifges(7,2) + Ifges(7,6) * t190;
t5 = -qJD(6) * t31 - t153 * t29 + t156 * t25;
t4 = qJD(6) * t30 + t153 * t25 + t156 * t29;
t1 = [(mrSges(5,2) * t201 + t184 * t264 + t182 * t263 - t151 * t46 + t40 * t249 + t41 * t248 + (mrSges(5,3) + t185) * t91 - t225 * qJD(2) + t177 * mrSges(6,3) + m(6) * (-qJD(2) * t108 - t230) - m(5) * t230 + (t56 * t249 - t157 * t55 / 0.2e1 + t108 * t186 + t181 * t257 + t183 * t255 + t179 * t251 + t175 * mrSges(6,3)) * qJD(5) + (t151 * t129 + (Ifges(7,5) * t261 + Ifges(7,6) * t262 + (-0.3e1 / 0.2e1 * Ifges(5,4) + t236 / 0.2e1 - t235 / 0.2e1) * t158 + t152 * mrSges(5,1) + (-0.3e1 / 0.2e1 * Ifges(5,1) + Ifges(7,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(5,2) + Ifges(6,3) / 0.2e1) * t155) * qJD(1) + t191 + t161) * qJD(4)) * t158 + t5 * t49 + t4 * t48 + t37 * t28 / 0.2e1 + t30 * t15 + t31 * t16 + t39 * t286 + m(7) * (t10 * t4 + t115 * t47 + t2 * t31 + t3 * t30 + t5 * t9 + t70 * t72) + m(5) * (qJD(2) * t273 + t152 * t201 + t209) + (t10 * t39 - t2 * t96 + t3 * t98 - t37 * t9) * mrSges(7,3) + t47 * (mrSges(7,1) * t96 - mrSges(7,2) * t98) + (-Ifges(7,4) * t98 - Ifges(7,2) * t96) * t269 + (-Ifges(7,1) * t98 - Ifges(7,4) * t96) * t270 + m(4) * (qJD(2) * t146 + t209 + (qJ(2) * qJD(2) + qJD(3) * t152) * qJD(1)) + m(6) * (t21 * t81 + t22 * t80 + t34 * t52 + t35 * t51) + (mrSges(5,1) * t201 + t141 / 0.2e1 + t140 / 0.2e1 + qJD(2) * t129 + (m(5) * t151 - mrSges(5,3)) * t92 + (0.3e1 / 0.2e1 * t197 + t192 + (-t135 - t276) * mrSges(5,2) + (t225 + t285) * t151 - t290) * qJD(4) + t272) * t155 + 0.2e1 * (qJD(3) * mrSges(4,3) + qJD(2) * t188) * qJD(1) + t70 * (-mrSges(7,1) * t39 + mrSges(7,2) * t37) + t72 * t32 + t80 * t59 + t81 * t60 + t34 * t82 + t35 * t83 + t115 * t8 + qJD(3) * t121 + (Ifges(7,5) * t37 + Ifges(7,6) * t39) * t252 + t7 * t261 + t6 * t262 + (Ifges(7,1) * t37 + Ifges(7,4) * t39) * t265 + (Ifges(7,4) * t37 + Ifges(7,2) * t39) * t267; t119 * t15 - t120 * t16 - t154 * t60 - t157 * t59 + t241 * t49 + t242 * t48 + t173 * qJD(5) - t188 * t159 + m(6) * (qJD(5) * t175 + t177) + ((-qJD(3) - t146) * m(4) + (-t169 + t221) * t155 - m(6) * (t215 * t52 - t217 * t51) + (t198 - t222 + t259 + t285) * t158 + (-qJD(3) - t273) * m(5)) * qJD(1) + (t10 * t242 + t119 * t3 - t120 * t2 + t241 * t9) * m(7); -t159 * mrSges(4,3) - t95 * t15 - t97 * t16 + t277 * t49 + t278 * t48 + (qJD(4) * t169 - t46 - t8) * t158 + (qJD(4) * t198 + qJD(5) * t172 + t174) * t155 + m(5) * (t155 * t92 - t158 * t91) + m(6) * ((-qJD(4) * t175 - t91) * t158 + (qJD(4) * t108 - qJD(5) * t176 + t178) * t155) + (-m(5) * t135 - t121 + (qJD(2) - t135) * m(4) + t274) * qJD(1) + (t10 * t278 - t158 * t47 - t2 * t97 + t208 * t70 + t277 * t9 - t3 * t95) * m(7); (pkin(8) * t274 + t246 * t279 + t296) * qJD(5) + (-t10 * t241 - t119 * t2 - t120 * t3 + t242 * t9) * mrSges(7,3) + (mrSges(7,1) * t241 - mrSges(7,2) * t242) * t70 + (-t158 * t129 - t155 * t225) * t134 + t178 * mrSges(6,3) - pkin(4) * t46 + t293 * t49 + t174 * pkin(8) + (-t68 / 0.2e1 + t89 / 0.2e1) * t28 + ((t192 + t135 * mrSges(5,2) - t197 / 0.2e1 + t290) * t155 + ((t240 / 0.2e1 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t155) * qJD(1) + t191 - t161 + (Ifges(7,5) * t120 - Ifges(7,6) * t119 + t179) * qJD(4) / 0.2e1) * t158) * qJD(1) + (-mrSges(5,1) - t186) * t91 + (-t292 / 0.2e1 - t88 / 0.2e1) * t27 + (-Ifges(7,5) * t68 - Ifges(7,6) * t292) * t252 + (-Ifges(7,1) * t68 - Ifges(7,4) * t292) * t265 + (-Ifges(7,4) * t68 - Ifges(7,2) * t292) * t267 + t294 * t48 + (t10 * t294 + t145 * t47 + t2 * t76 + t293 * t9 + t3 * t75 - t70 * t90) * m(7) + t75 * t15 + t76 * t16 - t67 * t82 - t66 * t83 - t90 * t32 - t92 * mrSges(5,2) - t119 * t6 / 0.2e1 + t120 * t7 / 0.2e1 + t47 * (mrSges(7,1) * t119 + mrSges(7,2) * t120) + t145 * t8 + t40 * t248 + (-Ifges(7,5) * t89 + Ifges(7,6) * t88) * t253 + t181 * t263 + t183 * t264 + (-Ifges(7,1) * t89 + Ifges(7,4) * t88) * t266 + (-Ifges(7,4) * t89 + Ifges(7,2) * t88) * t268 + (Ifges(7,4) * t120 - Ifges(7,2) * t119) * t269 + (Ifges(7,1) * t120 - Ifges(7,4) * t119) * t270 + t154 * t41 / 0.2e1 + (-pkin(4) * t91 + pkin(8) * t178 - t108 * t125 - t51 * t66 - t52 * t67) * m(6); -t12 * t49 - t13 * t48 + (-Ifges(6,2) * t118 + t113 + t56) * t257 + (m(7) * t10 * t203 + (m(7) * t3 + qJD(6) * t48 + t15) * t156 + (t16 - t49 * qJD(6) + m(7) * (t2 - t280)) * t153 - t279 * t118) * pkin(5) - m(7) * (t10 * t13 + t12 * t9) + t272 + (t117 * t51 + t118 * t52) * mrSges(6,3) + t64 * t286 + t141 - t51 * t82 + t52 * t83 - t108 * (mrSges(6,1) * t118 + mrSges(6,2) * t117) + (Ifges(6,5) * t117 - Ifges(6,6) * t118) * t251 + t55 * t254 + (Ifges(6,1) * t117 - t239) * t255 + t295; t10 * t49 + t27 * t265 - t9 * t48 + t288 + t295;];
tauc  = t1(:);
