% Calculate time derivative of joint inertia matrix for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:55:49
% EndTime: 2019-03-08 23:56:02
% DurationCPUTime: 5.52s
% Computational Cost: add. (4135->417), mult. (10230->591), div. (0->0), fcn. (9387->10), ass. (0->181)
t308 = Ifges(6,4) + Ifges(7,4);
t307 = Ifges(6,1) + Ifges(7,1);
t306 = Ifges(6,2) + Ifges(7,2);
t166 = cos(qJ(5));
t305 = t308 * t166;
t162 = sin(qJ(5));
t304 = t308 * t162;
t258 = Ifges(6,5) + Ifges(7,5);
t295 = Ifges(6,6) + Ifges(7,6);
t302 = -t306 * t162 + t305;
t301 = t307 * t166 - t304;
t300 = t162 * t258 + t295 * t166;
t298 = t166 * t306 + t304;
t297 = t162 * t307 + t305;
t160 = sin(pkin(6));
t165 = sin(qJ(2));
t225 = qJD(2) * t165;
t211 = t160 * t225;
t163 = sin(qJ(4));
t167 = cos(qJ(4));
t161 = cos(pkin(6));
t164 = sin(qJ(3));
t168 = cos(qJ(3));
t232 = t160 * t165;
t112 = t161 * t168 - t164 * t232;
t113 = t161 * t164 + t168 * t232;
t181 = t167 * t112 - t113 * t163;
t169 = cos(qJ(2));
t224 = qJD(2) * t169;
t210 = t160 * t224;
t97 = -qJD(3) * t113 - t164 * t210;
t98 = qJD(3) * t112 + t168 * t210;
t25 = qJD(4) * t181 + t163 * t97 + t167 * t98;
t231 = t160 * t169;
t67 = t112 * t163 + t113 * t167;
t51 = -t162 * t67 - t166 * t231;
t10 = qJD(5) * t51 + t162 * t211 + t166 * t25;
t222 = qJD(5) * t166;
t223 = qJD(5) * t162;
t126 = mrSges(7,1) * t223 + mrSges(7,2) * t222;
t190 = mrSges(6,1) * t162 + mrSges(6,2) * t166;
t127 = t190 * qJD(5);
t137 = -mrSges(7,1) * t166 + mrSges(7,2) * t162;
t179 = t162 * t231 - t166 * t67;
t11 = qJD(5) * t179 - t162 * t25 + t166 * t211;
t174 = -t11 * t162 + qJD(5) * (t162 * t179 - t166 * t51);
t242 = t166 * mrSges(7,3);
t243 = t166 * mrSges(6,3);
t26 = qJD(4) * t67 + t163 * t98 - t167 * t97;
t138 = -mrSges(6,1) * t166 + mrSges(6,2) * t162;
t290 = t138 - mrSges(5,1);
t296 = -t25 * mrSges(5,2) + (-t126 - t127) * t181 + (t242 + t243) * t10 + (t137 + t290) * t26 + (mrSges(6,3) + mrSges(7,3)) * t174;
t282 = (t162 ^ 2 + t166 ^ 2) * t167;
t294 = Ifges(6,3) + Ifges(7,3);
t124 = t163 * t168 + t164 * t167;
t123 = t163 * t164 - t167 * t168;
t280 = qJD(3) + qJD(4);
t90 = t280 * t123;
t241 = t166 * t90;
t177 = t124 * t223 + t241;
t208 = t124 * t222;
t246 = t162 * t90;
t178 = t208 - t246;
t91 = t280 * t124;
t293 = -t177 * t308 - t178 * t306 + t295 * t91;
t292 = -t177 * t307 - t178 * t308 + t258 * t91;
t291 = t295 * t123 + t302 * t124;
t254 = t258 * t123 + t301 * t124;
t288 = t302 * qJD(5);
t287 = t301 * qJD(5);
t226 = t258 * t222;
t267 = -pkin(9) - pkin(8);
t144 = t267 * t164;
t145 = t267 * t168;
t283 = t167 * t144 + t145 * t163;
t171 = t287 * t162 + t288 * t166 + t222 * t297 - t223 * t298;
t279 = 2 * m(6);
t278 = 2 * m(7);
t277 = -2 * mrSges(5,3);
t107 = t144 * t163 - t145 * t167;
t212 = qJD(3) * t267;
t134 = t164 * t212;
t198 = t168 * t212;
t45 = qJD(4) * t107 + t134 * t163 - t167 * t198;
t276 = 0.2e1 * t45;
t275 = -0.2e1 * t283;
t274 = 0.2e1 * t126;
t273 = 0.2e1 * t137;
t272 = -0.2e1 * t162;
t271 = 0.2e1 * t166;
t270 = m(5) / 0.2e1;
t269 = m(6) / 0.2e1;
t266 = m(6) * t26;
t265 = mrSges(7,3) * pkin(5);
t261 = pkin(3) * t167;
t260 = t181 * t26;
t218 = pkin(3) * qJD(3) * t164;
t37 = pkin(4) * t91 + pkin(10) * t90 + t218;
t44 = qJD(4) * t283 + t167 * t134 + t163 * t198;
t202 = -t162 * t44 + t166 * t37;
t150 = -pkin(3) * t168 - pkin(2);
t78 = pkin(4) * t123 - pkin(10) * t124 + t150;
t99 = t166 * t107;
t40 = t162 * t78 + t99;
t7 = -qJD(5) * t40 + t202;
t259 = t7 * t162;
t256 = -qJ(6) - pkin(10);
t249 = pkin(3) * qJD(4);
t248 = t283 * t45;
t245 = t163 * mrSges(5,1);
t244 = t163 * t181;
t240 = t167 * mrSges(5,2);
t239 = t283 * t163;
t238 = t124 * t162;
t237 = t124 * t166;
t147 = pkin(3) * t163 + pkin(10);
t235 = t147 * t166;
t230 = t162 * t167;
t229 = t163 * t138;
t156 = t166 * qJ(6);
t228 = t166 * t167;
t227 = -qJ(6) - t147;
t221 = 0.2e1 * t164;
t220 = t162 * t37 + t166 * t44 + t78 * t222;
t219 = 0.2e1 * qJD(5);
t217 = t167 * t249;
t216 = pkin(5) * t223;
t215 = m(7) * pkin(5) + mrSges(7,1);
t149 = -pkin(5) * t166 - pkin(4);
t205 = t295 * t162;
t204 = -t223 / 0.2e1;
t201 = qJD(5) * t256;
t200 = -t241 * t258 + t294 * t91;
t39 = -t107 * t162 + t166 * t78;
t199 = qJD(5) * t227;
t197 = t160 ^ 2 * t165 * t224;
t194 = mrSges(6,3) * t282;
t193 = -(2 * Ifges(5,4)) - t205;
t192 = -t181 * t45 - t26 * t283;
t191 = -mrSges(4,1) * t168 + mrSges(4,2) * t164;
t29 = pkin(5) * t123 - t124 * t156 + t39;
t34 = -qJ(6) * t238 + t40;
t185 = -t162 * t34 - t166 * t29;
t184 = -t162 * t40 - t166 * t39;
t182 = qJ(6) * t90 - qJD(6) * t124;
t27 = mrSges(7,1) * t178 - mrSges(7,2) * t177;
t173 = -t97 * t164 + t98 * t168 + (-t112 * t168 - t113 * t164) * qJD(3);
t172 = m(6) * (t10 * t166 + t174);
t23 = pkin(5) * t178 + t45;
t4 = -qJ(6) * t208 + (-qJD(5) * t107 + t182) * t162 + t220;
t6 = -t107 * t223 + t220;
t63 = pkin(5) * t238 - t283;
t170 = -t44 * mrSges(5,2) - Ifges(5,5) * t90 - t283 * t127 + t63 * t126 + t23 * t137 + t4 * t242 + t6 * t243 + t290 * t45 + (-t223 * t295 + t226) * t123 / 0.2e1 + t292 * t162 / 0.2e1 + t293 * t166 / 0.2e1 - t288 * t238 / 0.2e1 + t287 * t237 / 0.2e1 + t291 * t204 + t254 * t222 / 0.2e1 + t298 * (-t208 / 0.2e1 + t246 / 0.2e1) + t297 * (t124 * t204 - t241 / 0.2e1) + (-Ifges(5,6) + t300 / 0.2e1) * t91;
t155 = t166 * qJD(6);
t148 = -pkin(4) - t261;
t139 = pkin(10) * t166 + t156;
t136 = t256 * t162;
t135 = t149 - t261;
t133 = t163 * t249 + t216;
t128 = (mrSges(4,1) * t164 + mrSges(4,2) * t168) * qJD(3);
t120 = t156 + t235;
t119 = t227 * t162;
t111 = -qJD(6) * t162 + t166 * t201;
t110 = t162 * t201 + t155;
t96 = mrSges(5,1) * t123 + mrSges(5,2) * t124;
t89 = (-qJD(6) - t217) * t162 + t166 * t199;
t88 = t162 * t199 + t166 * t217 + t155;
t83 = mrSges(6,1) * t123 - mrSges(6,3) * t237;
t82 = mrSges(7,1) * t123 - mrSges(7,3) * t237;
t81 = -mrSges(6,2) * t123 - mrSges(6,3) * t238;
t80 = -mrSges(7,2) * t123 - mrSges(7,3) * t238;
t73 = t190 * t124;
t72 = (mrSges(7,1) * t162 + mrSges(7,2) * t166) * t124;
t41 = mrSges(5,1) * t91 - mrSges(5,2) * t90;
t33 = -mrSges(6,2) * t91 - mrSges(6,3) * t178;
t32 = -mrSges(7,2) * t91 - mrSges(7,3) * t178;
t31 = mrSges(6,1) * t91 + mrSges(6,3) * t177;
t30 = mrSges(7,1) * t91 + mrSges(7,3) * t177;
t28 = mrSges(6,1) * t178 - mrSges(6,2) * t177;
t2 = pkin(5) * t91 + t182 * t166 + (-t99 + (qJ(6) * t124 - t78) * t162) * qJD(5) + t202;
t1 = [0.2e1 * m(5) * (t67 * t25 - t197 - t260) + 0.2e1 * m(4) * (t112 * t97 + t113 * t98 - t197) + 0.4e1 * (m(7) / 0.2e1 + t269) * (-t10 * t179 + t11 * t51 - t260); -(t27 + t28) * t181 - (t32 + t33) * t179 + (t30 + t31) * t51 + (t72 + t73) * t26 + (t82 + t83) * t11 + (t80 + t81) * t10 + (-t123 * t25 + t124 * t26 + t181 * t90 - t67 * t91) * mrSges(5,3) + t173 * mrSges(4,3) + ((-t128 - t41) * t169 + (-t169 * mrSges(3,2) + (-mrSges(3,1) + t191 + t96) * t165) * qJD(2)) * t160 + m(5) * (t107 * t25 + t44 * t67 + (t150 * t225 - t169 * t218) * t160 + t192) + m(7) * (t10 * t34 + t11 * t29 - t179 * t4 - t181 * t23 + t2 * t51 + t26 * t63) + m(6) * (t10 * t40 + t11 * t39 - t179 * t6 + t51 * t7 + t192) + (-pkin(2) * t211 + pkin(8) * t173) * m(4); t107 * t91 * t277 - 0.2e1 * pkin(2) * t128 + t28 * t275 + 0.2e1 * t150 * t41 + 0.2e1 * t2 * t82 + 0.2e1 * t23 * t72 + 0.2e1 * t63 * t27 + 0.2e1 * t29 * t30 + 0.2e1 * t39 * t31 + 0.2e1 * t34 * t32 + 0.2e1 * t40 * t33 + 0.2e1 * t4 * t80 + t73 * t276 + 0.2e1 * t6 * t81 + 0.2e1 * t7 * t83 + 0.2e1 * m(5) * (t107 * t44 + t150 * t218 - t248) + (t2 * t29 + t23 * t63 + t34 * t4) * t278 + (t39 * t7 + t40 * t6 - t248) * t279 - (mrSges(5,3) * t275 - t162 * t291 + t254 * t166) * t90 + (t44 * t277 + ((2 * Ifges(5,2)) + t294) * t91 - t193 * t90 + t200) * t123 + (mrSges(5,3) * t276 - 0.2e1 * Ifges(5,1) * t90 + t292 * t166 - t293 * t162 + (t166 * t258 + t193) * t91 + (-t123 * t300 - t162 * t254 - t166 * t291) * qJD(5)) * t124 + ((-Ifges(4,4) * t164 + pkin(3) * t96) * t221 + (0.2e1 * Ifges(4,4) * t168 + (Ifges(4,1) - Ifges(4,2)) * t221) * t168) * qJD(3); 0.2e1 * ((t163 * t25 - t167 * t26) * t270 + ((-t179 * t228 - t230 * t51 - t244) * t269 + (t167 * t67 - t244) * t270) * qJD(4)) * pkin(3) + t147 * t172 + m(7) * (t10 * t120 + t11 * t119 - t133 * t181 + t135 * t26 - t179 * t88 + t51 * t89) + t148 * t266 + t97 * mrSges(4,1) - t98 * mrSges(4,2) + t296; (Ifges(4,5) * t168 - Ifges(4,6) * t164 + pkin(8) * t191) * qJD(3) + (-t7 * mrSges(6,3) - t2 * mrSges(7,3) - t147 * t31) * t162 + m(6) * (-t147 * t259 + t148 * t45 + t6 * t235) + m(7) * (t119 * t2 + t120 * t4 + t133 * t63 + t135 * t23 + t29 * t89 + t34 * t88) + t148 * t28 + t135 * t27 + t133 * t72 + t119 * t30 + t120 * t32 + t88 * t80 + t89 * t82 + t170 + (m(5) * (t163 * t44 - t167 * t45) + (-t163 * t91 + t167 * t90) * mrSges(5,3) + ((t124 * mrSges(5,3) + t73) * t163 + (-t123 * mrSges(5,3) - t162 * t83 + t166 * t81) * t167 + m(6) * (t228 * t40 - t230 * t39 - t239) + m(5) * (t107 * t167 - t239)) * qJD(4)) * pkin(3) + (t185 * mrSges(7,3) + t184 * mrSges(6,3) + (m(6) * t184 - t162 * t81 - t166 * t83) * t147) * qJD(5) + t33 * t235; 0.2e1 * t148 * t127 + t133 * t273 + t135 * t274 + (t119 * t89 + t120 * t88 + t133 * t135) * t278 + (t89 * t272 + t88 * t271 + (-t119 * t166 - t120 * t162) * t219) * mrSges(7,3) + (-0.2e1 * t240 - 0.2e1 * t245 + 0.2e1 * t229 + (t147 * t282 + t148 * t163) * t279 + 0.2e1 * t194) * t249 + t171; m(7) * (t10 * t139 + t11 * t136 - t110 * t179 + t111 * t51 + t149 * t26 - t181 * t216) - pkin(4) * t266 + pkin(10) * t172 + t296; t72 * t216 + (qJD(5) * t184 - t259) * mrSges(6,3) + (qJD(5) * t185 - t2 * t162) * mrSges(7,3) + m(7) * (t110 * t34 + t111 * t29 + t136 * t2 + t139 * t4 + t149 * t23 + t216 * t63) + t149 * t27 + t136 * t30 + t139 * t32 + t110 * t80 + t111 * t82 + t170 + (-t83 * t222 - t81 * t223 + m(6) * (t166 * t6 - t222 * t39 - t223 * t40 - t259) + t166 * t33 - t162 * t31) * pkin(10) + (-m(6) * t45 - t28) * pkin(4); m(7) * (t110 * t120 + t111 * t119 + t133 * t149 + t135 * t216 + t136 * t89 + t139 * t88) + (t133 + t216) * t137 + (t148 - pkin(4)) * t127 + (t149 + t135) * t126 + (-t240 - t245 + t229 + m(6) * (-pkin(4) * t163 + pkin(10) * t282) + t194) * t249 + ((t110 + t88) * t166 + (-t111 - t89) * t162 + ((-t119 - t136) * t166 + (-t120 - t139) * t162) * qJD(5)) * mrSges(7,3) + t171; -0.2e1 * pkin(4) * t127 + t216 * t273 + t149 * t274 + (t110 * t139 + t111 * t136 + t149 * t216) * t278 + (t110 * t271 + t111 * t272 + (-t136 * t166 - t139 * t162) * t219) * mrSges(7,3) + t171; (-mrSges(6,2) - mrSges(7,2)) * t10 + (mrSges(6,1) + t215) * t11; mrSges(6,1) * t7 + mrSges(7,1) * t2 - mrSges(6,2) * t6 - mrSges(7,2) * t4 + t90 * t205 + (m(7) * t2 + t30) * pkin(5) - t300 * t124 * qJD(5) + t200; -mrSges(7,2) * t88 + t215 * t89 - t190 * t217 + ((-mrSges(6,1) * t147 - t265) * t166 + (mrSges(6,2) * t147 - t295) * t162) * qJD(5) + t226; -mrSges(7,2) * t110 + t215 * t111 + ((-mrSges(6,1) * pkin(10) - t265) * t166 + (mrSges(6,2) * pkin(10) - t295) * t162) * qJD(5) + t226; 0; m(7) * t26; m(7) * t23 + t27; m(7) * t133 + t126; m(7) * t216 + t126; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
