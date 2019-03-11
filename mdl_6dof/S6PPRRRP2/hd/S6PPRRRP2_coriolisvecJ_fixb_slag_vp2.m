% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PPRRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:56:13
% EndTime: 2019-03-08 18:56:23
% DurationCPUTime: 5.12s
% Computational Cost: add. (4291->447), mult. (11623->613), div. (0->0), fcn. (9251->12), ass. (0->216)
t277 = Ifges(6,1) + Ifges(7,1);
t266 = Ifges(7,4) + Ifges(6,5);
t278 = qJD(4) / 0.2e1;
t276 = Ifges(6,6) - Ifges(7,6);
t143 = cos(qJ(4));
t201 = qJD(3) * t143;
t131 = Ifges(5,4) * t201;
t181 = Ifges(5,5) * t278;
t139 = sin(qJ(5));
t142 = cos(qJ(5));
t194 = t142 * qJD(4);
t140 = sin(qJ(4));
t203 = qJD(3) * t140;
t113 = t139 * t203 - t194;
t114 = qJD(4) * t139 + t142 * t203;
t138 = cos(pkin(6));
t125 = qJD(1) * t138 + qJD(2);
t133 = sin(pkin(12));
t135 = sin(pkin(6));
t141 = sin(qJ(3));
t144 = cos(qJ(3));
t136 = cos(pkin(12));
t137 = cos(pkin(7));
t207 = t136 * t137;
t151 = (t133 * t144 + t141 * t207) * t135;
t134 = sin(pkin(7));
t209 = t134 * t141;
t59 = qJD(1) * t151 + t125 * t209;
t56 = qJD(3) * pkin(9) + t59;
t204 = qJD(1) * t135;
t184 = t136 * t204;
t95 = t125 * t137 - t134 * t184;
t32 = -t140 * t56 + t143 * t95;
t30 = -qJD(4) * pkin(4) - t32;
t15 = pkin(5) * t113 - qJ(6) * t114 + t30;
t217 = t140 * t95;
t33 = t143 * t56 + t217;
t31 = qJD(4) * pkin(10) + t33;
t122 = -pkin(4) * t143 - pkin(10) * t140 - pkin(3);
t210 = t133 * t141;
t58 = t144 * (t125 * t134 + t137 * t184) - t204 * t210;
t46 = qJD(3) * t122 - t58;
t10 = -t139 * t31 + t142 * t46;
t11 = t139 * t46 + t142 * t31;
t160 = t10 * t142 + t11 * t139;
t219 = Ifges(7,5) * t142;
t164 = Ifges(7,3) * t139 + t219;
t222 = Ifges(6,4) * t142;
t168 = -Ifges(6,2) * t139 + t222;
t173 = mrSges(7,1) * t139 - mrSges(7,3) * t142;
t175 = mrSges(6,1) * t139 + mrSges(6,2) * t142;
t128 = qJD(5) - t201;
t263 = qJD(6) - t10;
t8 = -pkin(5) * t128 + t263;
t9 = qJ(6) * t128 + t11;
t177 = t9 * t139 - t8 * t142;
t241 = t142 / 0.2e1;
t243 = t139 / 0.2e1;
t244 = -t139 / 0.2e1;
t247 = t114 / 0.2e1;
t249 = t113 / 0.2e1;
t250 = -t113 / 0.2e1;
t111 = Ifges(6,4) * t113;
t221 = Ifges(7,5) * t113;
t264 = t277 * t114 + t266 * t128 - t111 + t221;
t267 = t128 / 0.2e1;
t220 = Ifges(7,5) * t139;
t223 = Ifges(6,4) * t139;
t269 = t277 * t142 + t220 - t223;
t110 = Ifges(7,5) * t114;
t60 = Ifges(7,6) * t128 + Ifges(7,3) * t113 + t110;
t224 = Ifges(6,4) * t114;
t63 = -Ifges(6,2) * t113 + Ifges(6,6) * t128 + t224;
t258 = t177 * mrSges(7,2) + t160 * mrSges(6,3) - t15 * t173 - t164 * t249 - t168 * t250 - t30 * t175 - t60 * t243 - t63 * t244 - t269 * t247 - (-t139 * t276 + t142 * t266) * t267 - t264 * t241;
t274 = t203 / 0.2e1;
t55 = -qJD(3) * pkin(3) - t58;
t275 = t55 * mrSges(5,2) + Ifges(5,1) * t274 + t131 / 0.2e1 + t181 - t32 * mrSges(5,3) - t258;
t193 = qJD(3) * qJD(4);
t179 = t140 * t193;
t192 = qJD(4) * qJD(5);
t197 = qJD(5) * t139;
t89 = t142 * t192 + (-t140 * t197 + t143 * t194) * qJD(3);
t196 = qJD(5) * t142;
t198 = qJD(4) * t143;
t90 = t139 * t192 + (t139 * t198 + t140 * t196) * qJD(3);
t273 = (-Ifges(6,4) + Ifges(7,5)) * t90 + t277 * t89 + t266 * t179;
t272 = qJD(4) * t33;
t152 = (t144 * t207 - t210) * t135;
t208 = t134 * t144;
t271 = t138 * t208 + t152;
t270 = t277 * t139 - t219 + t222;
t268 = -m(5) * t32 + m(6) * t30;
t180 = -Ifges(5,6) * qJD(4) / 0.2e1;
t262 = t266 * t139 + t276 * t142;
t52 = (qJD(1) * t152 + t125 * t208) * qJD(3);
t13 = qJD(4) * t32 + t143 * t52;
t178 = pkin(4) * t140 - pkin(10) * t143;
t117 = t178 * qJD(4);
t41 = (t117 + t59) * qJD(3);
t3 = t142 * t13 + t139 * t41 + t46 * t196 - t197 * t31;
t4 = -qJD(5) * t11 - t13 * t139 + t142 * t41;
t261 = -t139 * t4 + t142 * t3;
t1 = qJ(6) * t179 + qJD(6) * t128 + t3;
t2 = -pkin(5) * t179 - t4;
t260 = t1 * t142 + t139 * t2;
t195 = qJD(5) * t143;
t199 = qJD(4) * t140;
t50 = pkin(9) * (t139 * t199 - t142 * t195) + t117 * t142 - t122 * t197;
t188 = -Ifges(6,3) / 0.2e1 - Ifges(7,2) / 0.2e1;
t189 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t190 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t257 = -t189 * t113 - t190 * t114 + t188 * t128 - t10 * mrSges(6,1) - t55 * mrSges(5,1) - t9 * mrSges(7,3) - t180 + (Ifges(5,4) * t140 + t143 * Ifges(5,2)) * qJD(3) / 0.2e1 - Ifges(6,6) * t250 - Ifges(7,6) * t249 + t11 * mrSges(6,2) + t33 * mrSges(5,3) + t8 * mrSges(7,1) - (Ifges(6,3) + Ifges(7,2)) * t267 - t266 * t247;
t256 = 0.2e1 * m(5);
t255 = t89 / 0.2e1;
t254 = -t90 / 0.2e1;
t253 = t90 / 0.2e1;
t248 = -t114 / 0.2e1;
t246 = -t128 / 0.2e1;
t242 = -t142 / 0.2e1;
t14 = t140 * t52 + t272;
t101 = -t134 * t135 * t136 + t137 * t138;
t74 = t138 * t209 + t151;
t158 = t101 * t143 - t140 * t74;
t237 = t14 * t158;
t53 = qJD(3) * t59;
t234 = t53 * t271;
t44 = mrSges(7,1) * t90 - mrSges(7,3) * t89;
t45 = mrSges(6,1) * t90 + mrSges(6,2) * t89;
t231 = t44 + t45;
t69 = -mrSges(7,2) * t90 + mrSges(7,3) * t179;
t72 = -mrSges(6,2) * t179 - mrSges(6,3) * t90;
t230 = t69 + t72;
t70 = mrSges(6,1) * t179 - mrSges(6,3) * t89;
t71 = -mrSges(7,1) * t179 + t89 * mrSges(7,2);
t229 = -t70 + t71;
t226 = mrSges(6,3) * t113;
t91 = -mrSges(6,2) * t128 - t226;
t94 = -t113 * mrSges(7,2) + mrSges(7,3) * t128;
t228 = t91 + t94;
t225 = mrSges(6,3) * t114;
t92 = mrSges(6,1) * t128 - t225;
t93 = -mrSges(7,1) * t128 + t114 * mrSges(7,2);
t227 = t92 - t93;
t102 = -t143 * t137 + t140 * t209;
t218 = t102 * t14;
t215 = t144 * t53;
t116 = t178 * qJD(3);
t23 = t139 * t116 + t142 * t32;
t214 = qJD(4) * mrSges(5,1) - mrSges(6,1) * t113 - mrSges(6,2) * t114 - mrSges(5,3) * t203;
t211 = t122 * t142;
t206 = t139 * t143;
t205 = t142 * t143;
t97 = pkin(9) * t205 + t139 * t122;
t202 = qJD(3) * t141;
t200 = qJD(3) * t144;
t76 = mrSges(7,1) * t113 - mrSges(7,3) * t114;
t191 = -t76 + t214;
t187 = mrSges(5,3) * t201;
t185 = t139 * t208;
t183 = t134 * t202;
t182 = t134 * t200;
t176 = mrSges(6,1) * t142 - mrSges(6,2) * t139;
t174 = mrSges(7,1) * t142 + mrSges(7,3) * t139;
t167 = Ifges(6,2) * t142 + t223;
t163 = -Ifges(7,3) * t142 + t220;
t162 = pkin(5) * t142 + qJ(6) * t139;
t161 = pkin(5) * t139 - qJ(6) * t142;
t43 = t101 * t140 + t143 * t74;
t21 = -t139 * t271 + t142 * t43;
t20 = t139 * t43 + t142 * t271;
t22 = t116 * t142 - t139 * t32;
t157 = pkin(9) + t161;
t103 = t137 * t140 + t143 * t209;
t80 = t103 * t139 + t142 * t208;
t150 = -t4 * mrSges(6,1) + t2 * mrSges(7,1) + t3 * mrSges(6,2) - t1 * mrSges(7,3);
t49 = t139 * t117 + t122 * t196 + (-t139 * t195 - t140 * t194) * pkin(9);
t145 = qJD(3) ^ 2;
t127 = Ifges(7,2) * t179;
t126 = Ifges(6,3) * t179;
t124 = -qJD(4) * mrSges(5,2) + t187;
t119 = -pkin(4) - t162;
t115 = (-mrSges(5,1) * t143 + mrSges(5,2) * t140) * qJD(3);
t108 = (mrSges(5,1) * t140 + mrSges(5,2) * t143) * t193;
t100 = qJD(5) * t161 - qJD(6) * t139;
t99 = t157 * t140;
t96 = -pkin(9) * t206 + t211;
t88 = -t211 + (pkin(9) * t139 + pkin(5)) * t143;
t87 = -qJ(6) * t143 + t97;
t86 = Ifges(7,4) * t89;
t85 = Ifges(6,5) * t89;
t84 = Ifges(6,6) * t90;
t83 = Ifges(7,6) * t90;
t81 = t103 * t142 - t185;
t79 = qJD(4) * t103 + t140 * t182;
t78 = -qJD(4) * t102 + t143 * t182;
t75 = pkin(5) * t114 + qJ(6) * t113;
t67 = t74 * qJD(3);
t66 = t271 * qJD(3);
t48 = (qJD(5) * t162 - qJD(6) * t142) * t140 + t157 * t198;
t47 = -pkin(5) * t199 - t50;
t39 = qJ(6) * t199 - qJD(6) * t143 + t49;
t35 = t89 * Ifges(6,4) - t90 * Ifges(6,2) + Ifges(6,6) * t179;
t34 = t89 * Ifges(7,5) + Ifges(7,6) * t179 + t90 * Ifges(7,3);
t29 = -qJD(5) * t185 + t103 * t196 + t139 * t78 - t142 * t183;
t28 = -qJD(5) * t80 + t139 * t183 + t142 * t78;
t26 = t139 * t59 + t205 * t58;
t25 = -t142 * t59 + t206 * t58;
t24 = t217 + (qJD(3) * t161 + t56) * t143;
t19 = -pkin(5) * t203 - t22;
t18 = qJ(6) * t203 + t23;
t17 = qJD(4) * t158 + t143 * t66;
t16 = qJD(4) * t43 + t140 * t66;
t7 = pkin(5) * t90 - qJ(6) * t89 - qJD(6) * t114 + t14;
t6 = -qJD(5) * t20 + t139 * t67 + t142 * t17;
t5 = qJD(5) * t21 + t139 * t17 - t67 * t142;
t12 = [-t271 * t108 + t67 * t115 + t17 * t124 + t228 * t6 - t227 * t5 - t231 * t158 + t230 * t21 + t229 * t20 - t191 * t16 + (-t67 * mrSges(4,1) - t66 * mrSges(4,2) + (-t140 * t43 - t143 * t158) * qJD(4) * mrSges(5,3)) * qJD(3) + m(5) * (t13 * t43 - t16 * t32 + t17 * t33 + t55 * t67 - t234 - t237) + m(6) * (-t10 * t5 + t11 * t6 + t16 * t30 - t20 * t4 + t21 * t3 - t237) + m(7) * (t1 * t21 + t15 * t16 - t158 * t7 + t2 * t20 + t5 * t8 + t6 * t9) + m(4) * (t52 * t74 - t58 * t67 + t59 * t66 - t234); -t103 * mrSges(5,3) * t179 + t78 * t124 + t230 * t81 + t229 * t80 - t227 * t29 + t228 * t28 - t191 * t79 + (qJD(4) * t187 + t231) * t102 + m(5) * (t103 * t13 - t32 * t79 + t33 * t78 + t218) + m(6) * (-t10 * t29 + t11 * t28 + t3 * t81 + t30 * t79 - t4 * t80 + t218) + m(7) * (t1 * t81 + t102 * t7 + t15 * t79 + t2 * t80 + t28 * t9 + t29 * t8) + ((-mrSges(4,2) * t145 - t108) * t144 + (-mrSges(4,1) * t145 + qJD(3) * t115) * t141 + m(4) * (t141 * t52 + t200 * t59 - t202 * t58 - t215) + m(5) * (t202 * t55 - t215)) * t134; -pkin(3) * t108 - t59 * t115 + t39 * t94 + t99 * t44 + t47 * t93 + t48 * t76 + t49 * t91 + t50 * t92 + t87 * t69 + t96 * t70 + t88 * t71 + t97 * t72 - t228 * t26 + t227 * t25 + (qJD(3) * t58 - t52) * mrSges(4,2) - m(7) * (t25 * t8 + t26 * t9) - m(6) * (-t10 * t25 + t11 * t26) + m(6) * (t10 * t50 + t11 * t49 + t3 * t97 + t4 * t96) + m(7) * (t1 * t87 + t15 * t48 + t2 * t88 + t39 * t9 + t47 * t8 + t7 * t99) + (-pkin(3) * t53 / 0.2e1 - t55 * t59 / 0.2e1) * t256 + (t13 * mrSges(5,3) - t58 * t124 - t53 * mrSges(5,1) - t126 / 0.2e1 - t127 / 0.2e1 - t86 / 0.2e1 + t84 / 0.2e1 - t85 / 0.2e1 - t83 / 0.2e1 - t189 * t90 - t190 * t89 + (pkin(9) * t13 / 0.2e1 - t33 * t58 / 0.2e1) * t256 + (t181 + 0.3e1 / 0.2e1 * t131 + (-t214 + t268) * pkin(9) + t275) * qJD(4) + t150) * t143 + (t34 * t243 + t35 * t244 + t7 * t173 + t53 * mrSges(5,2) + t164 * t253 + t168 * t254 + (mrSges(5,3) + t175) * t14 + (-t3 * t139 - t4 * t142) * mrSges(6,3) + (-t1 * t139 + t2 * t142) * mrSges(7,2) + (((-0.3e1 / 0.2e1 * Ifges(5,4) + t190 * t142 + t189 * t139) * t140 + (-0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(5,1) + t188) * t143) * qJD(3) + t180 - t257) * qJD(4) + (t30 * t176 + t15 * t174 + t163 * t250 + t167 * t249 + t63 * t242 + (t10 * t139 - t11 * t142) * mrSges(6,3) + (-t8 * t139 - t9 * t142) * mrSges(7,2) + t270 * t248 + t262 * t246 + t264 * t244) * qJD(5) + t269 * t255 + (t60 * qJD(5) + t273) * t241 + (m(6) * t14 - qJD(4) * t124 + t45 + (t14 - t272) * m(5)) * pkin(9) + (-m(7) * t15 + t191 - t268) * t58) * t140; t273 * t243 + t270 * t255 + (-t258 + (-m(6) * t160 - m(7) * t177 - t139 * t228 - t142 * t227) * pkin(10)) * qJD(5) + (-t24 + t100) * t76 - t13 * mrSges(5,2) - t7 * t174 + (-mrSges(5,1) - t176) * t14 - m(7) * (t15 * t24 + t18 * t9 + t19 * t8) - m(6) * (t10 * t22 + t11 * t23 + t30 * t33) + m(7) * (pkin(10) * t260 + t100 * t15 + t119 * t7) + t260 * mrSges(7,2) + m(6) * (-pkin(4) * t14 + pkin(10) * t261) + t261 * mrSges(6,3) - pkin(4) * t45 + ((Ifges(5,4) * t274 + t262 * t278 + t180 + t257) * t140 + ((-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t203 - t131 / 0.2e1 + t181 - t275) * t143) * qJD(3) - t23 * t91 - t22 * t92 - t19 * t93 - t18 * t94 + t119 * t44 - t32 * t124 + t35 * t241 + t34 * t242 + t163 * t253 + t167 * t254 + t214 * t33 + (t139 * t229 + t142 * t230) * pkin(10); -t150 + (t113 * t8 + t114 * t9) * mrSges(7,2) + t126 + t127 + t86 - t84 + t85 + t83 + qJ(6) * t69 - pkin(5) * t71 - t75 * t76 + qJD(6) * t94 - t15 * (mrSges(7,1) * t114 + mrSges(7,3) * t113) - t30 * (mrSges(6,1) * t114 - mrSges(6,2) * t113) + t63 * t247 + (Ifges(7,3) * t114 - t221) * t250 + (t225 + t227) * t11 + (-t226 - t228) * t10 + (-t113 * t266 - t114 * t276) * t246 + (-pkin(5) * t2 + qJ(6) * t1 - t11 * t8 - t15 * t75 + t263 * t9) * m(7) + (-Ifges(6,2) * t114 - t111 + t264) * t249 + (-t277 * t113 + t110 - t224 + t60) * t248; t114 * t76 - t128 * t94 + 0.2e1 * (t2 / 0.2e1 + t15 * t247 + t9 * t246) * m(7) + t71;];
tauc  = t12(:);
