% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:13
% EndTime: 2019-03-09 03:11:22
% DurationCPUTime: 4.88s
% Computational Cost: add. (3001->443), mult. (6820->561), div. (0->0), fcn. (3506->6), ass. (0->209)
t269 = Ifges(6,1) + Ifges(7,1);
t264 = Ifges(6,5) + Ifges(7,4);
t134 = sin(qJ(5));
t136 = cos(qJ(5));
t138 = -pkin(3) - pkin(8);
t137 = cos(qJ(3));
t126 = t137 * qJD(2);
t135 = sin(qJ(3));
t119 = sin(pkin(9)) * pkin(1) + pkin(7);
t105 = t119 * qJD(1);
t176 = pkin(4) * qJD(1) + t105;
t151 = t176 * t135;
t64 = t126 - t151;
t45 = qJD(3) * t138 + qJD(4) - t64;
t180 = -cos(pkin(9)) * pkin(1) - pkin(2);
t150 = -qJ(4) * t135 + t180;
t80 = t137 * t138 + t150;
t54 = t80 * qJD(1);
t12 = -t134 * t54 + t136 * t45;
t13 = t134 * t45 + t136 * t54;
t155 = t12 * t134 - t13 * t136;
t198 = qJD(1) * t135;
t118 = qJD(5) + t198;
t249 = qJD(6) - t12;
t10 = -pkin(5) * t118 + t249;
t11 = qJ(6) * t118 + t13;
t156 = t10 * t134 + t11 * t136;
t171 = mrSges(7,1) * t136 + mrSges(7,3) * t134;
t173 = mrSges(6,1) * t136 - mrSges(6,2) * t134;
t196 = qJD(3) * t134;
t197 = qJD(1) * t137;
t100 = t136 * t197 + t196;
t179 = t134 * t197;
t194 = qJD(3) * t136;
t101 = -t179 + t194;
t129 = qJD(3) * qJ(4);
t125 = t135 * qJD(2);
t84 = t137 * t105 + t125;
t65 = pkin(4) * t197 + t84;
t52 = t129 + t65;
t19 = pkin(5) * t100 - qJ(6) * t101 + t52;
t226 = t136 / 0.2e1;
t227 = -t136 / 0.2e1;
t229 = -t134 / 0.2e1;
t213 = Ifges(7,5) * t100;
t96 = Ifges(6,4) * t100;
t251 = t269 * t101 + t264 * t118 + t213 - t96;
t95 = Ifges(7,5) * t101;
t31 = t118 * Ifges(7,6) + t100 * Ifges(7,3) + t95;
t218 = Ifges(6,4) * t101;
t34 = -t100 * Ifges(6,2) + t118 * Ifges(6,6) + t218;
t270 = -t156 * mrSges(7,2) + t155 * mrSges(6,3) + t19 * t171 + t52 * t173 + t226 * t31 + t227 * t34 + t251 * t229;
t107 = t180 * qJD(1);
t212 = Ifges(7,5) * t134;
t160 = -Ifges(7,3) * t136 + t212;
t217 = Ifges(6,4) * t134;
t164 = Ifges(6,2) * t136 + t217;
t230 = t118 / 0.2e1;
t232 = t101 / 0.2e1;
t234 = t100 / 0.2e1;
t235 = -t100 / 0.2e1;
t211 = Ifges(7,5) * t136;
t216 = Ifges(6,4) * t136;
t246 = t269 * t134 - t211 + t216;
t209 = Ifges(7,6) * t136;
t210 = Ifges(6,6) * t136;
t214 = Ifges(6,5) * t134;
t215 = Ifges(7,4) * t134;
t247 = -t209 + t215 + t210 + t214;
t255 = qJD(3) / 0.2e1;
t256 = -qJD(3) / 0.2e1;
t257 = qJD(1) / 0.2e1;
t73 = -t129 - t84;
t94 = -pkin(3) * t137 + t150;
t74 = t94 * qJD(1);
t268 = t107 * mrSges(4,1) + t73 * mrSges(5,1) + Ifges(4,6) * t256 - (Ifges(4,4) * t135 + t137 * Ifges(4,2)) * qJD(1) / 0.2e1 + Ifges(5,5) * t255 + (-Ifges(5,6) * t135 - t137 * Ifges(5,3)) * t257 + t164 * t235 + t160 * t234 - t74 * mrSges(5,2) - t84 * mrSges(4,3) + t246 * t232 + t247 * t230 - t270;
t189 = qJD(1) * qJD(3);
t177 = t137 * t189;
t178 = t135 * t189;
t188 = qJD(3) * qJD(5);
t63 = -qJD(5) * t179 + (-t178 + t188) * t136;
t38 = -mrSges(7,2) * t63 + mrSges(7,3) * t177;
t190 = qJD(5) * t136;
t195 = qJD(3) * t135;
t62 = -t134 * t188 + (t134 * t195 - t137 * t190) * qJD(1);
t39 = mrSges(6,1) * t177 - mrSges(6,3) * t62;
t40 = -mrSges(7,1) * t177 + t62 * mrSges(7,2);
t41 = -mrSges(6,2) * t177 - mrSges(6,3) * t63;
t145 = (t38 + t41) * t134 + (t39 - t40) * t136;
t219 = mrSges(6,3) * t101;
t68 = mrSges(6,1) * t118 - t219;
t69 = -mrSges(7,1) * t118 + mrSges(7,2) * t101;
t221 = t68 - t69;
t208 = t100 * mrSges(6,3);
t66 = -mrSges(6,2) * t118 - t208;
t67 = -mrSges(7,2) * t100 + mrSges(7,3) * t118;
t222 = t66 + t67;
t243 = t134 * t221 - t136 * t222;
t267 = t243 * qJD(5) - t145;
t266 = -t197 / 0.2e1;
t265 = -mrSges(4,1) + mrSges(5,2);
t254 = mrSges(5,1) + mrSges(4,3);
t225 = pkin(4) + t119;
t263 = (-Ifges(6,4) + Ifges(7,5)) * t63 + t269 * t62 + t264 * t177;
t83 = t105 * t135 - t126;
t262 = -qJD(4) - t83;
t261 = t269 * t136 + t212 - t217;
t260 = -Ifges(4,1) / 0.2e1;
t233 = -t101 / 0.2e1;
t231 = -t118 / 0.2e1;
t258 = Ifges(4,4) * t266;
t252 = Ifges(6,6) - Ifges(7,6);
t77 = qJD(3) * t84;
t250 = -qJD(3) * t73 - t77;
t191 = qJD(5) * t134;
t117 = pkin(3) * t178;
t157 = pkin(8) * t135 - qJ(4) * t137;
t192 = qJD(4) * t135;
t142 = qJD(3) * t157 - t192;
t50 = qJD(1) * t142 + t117;
t51 = (t137 * t176 + t125) * qJD(3);
t3 = t134 * t51 + t136 * t50 + t45 * t190 - t191 * t54;
t4 = -qJD(5) * t13 - t134 * t50 + t136 * t51;
t174 = t134 * t3 + t136 * t4;
t245 = -qJD(3) * t52 + t174;
t1 = qJ(6) * t177 + qJD(6) * t118 + t3;
t2 = -pkin(5) * t177 - t4;
t175 = t1 * t134 - t136 * t2;
t244 = -qJD(3) * t19 + t175;
t97 = t225 * t135;
t220 = t134 * t97 + t136 * t80;
t124 = pkin(3) * t195;
t71 = t124 + t142;
t193 = qJD(3) * t137;
t88 = t225 * t193;
t9 = -qJD(5) * t220 - t134 * t71 + t136 * t88;
t181 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t182 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t184 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t70 = -qJD(3) * pkin(3) - t262;
t241 = -t181 * t100 - t184 * t101 - t182 * t118 + t10 * mrSges(7,1) + t13 * mrSges(6,2) + t74 * mrSges(5,3) + Ifges(6,6) * t234 + Ifges(7,6) * t235 + t198 * t260 + Ifges(4,5) * t256 + t258 + Ifges(5,4) * t255 + (-t135 * Ifges(5,2) - Ifges(5,6) * t137) * t257 - t107 * mrSges(4,2) - t11 * mrSges(7,3) - t12 * mrSges(6,1) - t70 * mrSges(5,1) - t83 * mrSges(4,3) + t264 * t233 + (Ifges(6,3) + Ifges(7,2)) * t231;
t240 = 0.2e1 * t119;
t239 = m(4) / 0.2e1;
t237 = -t63 / 0.2e1;
t236 = t63 / 0.2e1;
t228 = t134 / 0.2e1;
t122 = pkin(3) * t198;
t85 = qJD(1) * t157 + t122;
t25 = t134 * t65 + t136 * t85;
t109 = -mrSges(5,1) * t197 - qJD(3) * mrSges(5,3);
t48 = mrSges(6,1) * t100 + mrSges(6,2) * t101;
t207 = -t109 + t48;
t203 = t134 * t138;
t202 = t136 * t138;
t108 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t197;
t201 = t108 - t109;
t200 = t265 * qJD(3) + t254 * t198;
t98 = t225 * t137;
t47 = mrSges(7,1) * t100 - mrSges(7,3) * t101;
t187 = -t47 - t207;
t186 = -0.3e1 / 0.2e1 * Ifges(4,4) - 0.3e1 / 0.2e1 * Ifges(5,6);
t185 = -Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t183 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t120 = qJD(3) * t126;
t76 = -t105 * t195 + t120;
t172 = mrSges(6,1) * t134 + mrSges(6,2) * t136;
t170 = mrSges(7,1) * t134 - mrSges(7,3) * t136;
t165 = -Ifges(6,2) * t134 + t216;
t161 = Ifges(7,3) * t134 + t211;
t159 = pkin(5) * t136 + qJ(6) * t134;
t158 = -pkin(5) * t134 + qJ(6) * t136;
t24 = -t134 * t85 + t136 * t65;
t29 = -t134 * t80 + t136 * t97;
t149 = -pkin(4) - t159;
t8 = t134 * t88 + t136 * t71 + t97 * t190 - t191 * t80;
t146 = -qJ(4) * t193 - t192;
t143 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t128 = qJD(3) * qJD(4);
t42 = -qJD(3) * t151 + t120 + t128;
t116 = Ifges(7,2) * t177;
t115 = Ifges(6,3) * t177;
t104 = qJ(4) - t158;
t103 = -qJ(4) * t197 + t122;
t102 = (mrSges(5,2) * t137 - mrSges(5,3) * t135) * qJD(1);
t89 = t124 + t146;
t87 = t225 * t195;
t82 = qJD(5) * t159 - qJD(6) * t136 + qJD(4);
t75 = qJD(1) * t146 + t117;
t61 = -t128 - t76;
t59 = Ifges(7,4) * t62;
t58 = Ifges(6,5) * t62;
t57 = Ifges(6,6) * t63;
t56 = Ifges(7,6) * t63;
t49 = t137 * t159 + t98;
t46 = pkin(5) * t101 + qJ(6) * t100;
t28 = t126 + (qJD(1) * t149 - t105) * t135;
t27 = -pkin(5) * t135 - t29;
t26 = qJ(6) * t135 + t220;
t23 = mrSges(6,1) * t63 + mrSges(6,2) * t62;
t22 = mrSges(7,1) * t63 - mrSges(7,3) * t62;
t21 = -pkin(5) * t197 - t24;
t20 = qJ(6) * t197 + t25;
t18 = (qJD(5) * t158 + qJD(6) * t134) * t137 + (-t119 + t149) * t195;
t15 = t62 * Ifges(6,4) - t63 * Ifges(6,2) + Ifges(6,6) * t177;
t14 = t62 * Ifges(7,5) + Ifges(7,6) * t177 + t63 * Ifges(7,3);
t7 = -pkin(5) * t193 - t9;
t6 = pkin(5) * t63 - qJ(6) * t62 - qJD(6) * t101 + t42;
t5 = qJ(6) * t193 + qJD(6) * t135 + t8;
t16 = [t89 * t102 + t18 * t47 + t49 * t22 + t98 * t23 + t26 * t38 + t27 * t40 + t29 * t39 + t220 * t41 - t87 * t48 + t5 * t67 + t8 * t66 + t9 * t68 + t7 * t69 + m(6) * (t12 * t9 + t13 * t8 + t220 * t3 + t29 * t4 + t42 * t98 - t52 * t87) + m(7) * (t1 * t26 + t10 * t7 + t11 * t5 + t18 * t19 + t2 * t27 + t49 * t6) + m(5) * (t74 * t89 + t75 * t94) + (t115 / 0.2e1 + t116 / 0.2e1 + t58 / 0.2e1 + t59 / 0.2e1 + t56 / 0.2e1 - t57 / 0.2e1 - t75 * mrSges(5,3) + t181 * t63 + t184 * t62 + ((m(5) / 0.2e1 + t239) * t240 + t254) * t77 + (t183 * qJD(3) + (-m(4) * t84 + m(5) * t73 - t201) * t119 + (mrSges(4,1) * t180 - t94 * mrSges(5,2) + t135 * t186) * qJD(1) + t268) * qJD(3) + t143) * t135 + (t160 * t237 + t42 * t173 + t6 * t171 + t75 * mrSges(5,2) + t14 * t226 + t15 * t227 + t76 * mrSges(4,3) - t61 * mrSges(5,1) + t164 * t236 + (t134 * t4 - t136 * t3) * mrSges(6,3) + (-t1 * t136 - t134 * t2) * mrSges(7,2) + (-m(5) * t61 / 0.2e1 + t76 * t239) * t240 + (-t52 * t172 - t19 * t170 + t161 * t235 + t165 * t234 + t34 * t228 + (t12 * t136 + t13 * t134) * mrSges(6,3) + (-t10 * t136 + t11 * t134) * mrSges(7,2) + t261 * t233 + (t252 * t134 - t136 * t264) * t230 + t251 * t227) * qJD(5) + ((t180 * mrSges(4,2) - t94 * mrSges(5,3) + (-t215 / 0.2e1 + t209 / 0.2e1 - t214 / 0.2e1 - t210 / 0.2e1 - t186) * t137 + (0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(5,3) + t182) * t135) * qJD(1) + (m(4) * t83 + m(5) * t70 + t200) * t119 + t185 * qJD(3) - t241) * qJD(3) - t246 * t62 / 0.2e1 + (t31 * qJD(5) + t263) * t229) * t137; (t22 + t23 + (t134 * t222 + t136 * t221 + t200) * qJD(3) + m(6) * (t12 * t194 + t13 * t196 + t42) + m(7) * (-t10 * t194 + t11 * t196 + t6) + m(4) * (qJD(3) * t83 + t76) + m(5) * (qJD(3) * t70 - t61)) * t135 + ((t108 - t187) * qJD(3) + m(6) * (t12 * t191 - t13 * t190 - t245) + m(7) * (-t10 * t191 - t11 * t190 - t244) + m(5) * t250 + t267) * t137 + t254 * t189 * (-t135 ^ 2 - t137 ^ 2); ((t258 + (-pkin(3) * mrSges(5,1) + t134 * t181 + t136 * t184 + t185) * qJD(3) + Ifges(5,6) * t266 + t241) * t137 + ((Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t260 - Ifges(5,2) / 0.2e1) * t197 + (Ifges(4,4) / 0.2e1 + Ifges(5,6) / 0.2e1) * t198 + (-qJ(4) * mrSges(5,1) + t183) * qJD(3) - t268) * t135) * qJD(1) + (t164 * t234 + t160 * t235 + (-m(6) * t155 + m(7) * t156 - t243) * t138 + t246 * t233 + t247 * t231 + t270) * qJD(5) + t261 * t62 / 0.2e1 + (-pkin(3) * t77 - qJ(4) * t61 - t103 * t74 + t262 * t73 - t70 * t84) * m(5) + (t82 - t28) * t47 - t103 * t102 + t104 * t22 - t64 * t48 - t25 * t66 - t20 * t67 - t24 * t68 - t21 * t69 - t76 * mrSges(4,2) - t61 * mrSges(5,3) + qJ(4) * t23 - m(6) * (t12 * t24 + t13 * t25 + t52 * t64) - m(7) * (t10 * t21 + t11 * t20 + t19 * t28) + t145 * t138 + t161 * t236 + t165 * t237 + t14 * t228 + t15 * t229 + t6 * t170 + t42 * t172 - t174 * mrSges(6,3) - t175 * mrSges(7,2) + t263 * t226 + t265 * t77 - t200 * t84 + t201 * t83 + m(6) * (t42 * qJ(4) + t52 * qJD(4) + t202 * t4 + t203 * t3) + m(7) * (t1 * t203 + t6 * t104 + t19 * t82 - t2 * t202) + t207 * qJD(4); t187 * qJD(3) + (mrSges(5,1) * t193 + (t102 - t243) * t135) * qJD(1) + (t118 * t156 + t244) * m(7) + (-t118 * t155 + t245) * m(6) + (t198 * t74 - t250) * m(5) - t267; (t10 * t100 + t101 * t11) * mrSges(7,2) + t115 + t116 + t58 + t59 + t56 - t57 - t52 * (mrSges(6,1) * t101 - mrSges(6,2) * t100) - t19 * (mrSges(7,1) * t101 + mrSges(7,3) * t100) + qJD(6) * t67 - t46 * t47 + qJ(6) * t38 - pkin(5) * t40 + t143 + (Ifges(7,3) * t101 - t213) * t235 + t34 * t232 + (t219 + t221) * t13 + (-t208 - t222) * t12 + (-t100 * t264 - t252 * t101) * t231 + (-pkin(5) * t2 + qJ(6) * t1 - t10 * t13 + t249 * t11 - t19 * t46) * m(7) + (-Ifges(6,2) * t101 + t251 - t96) * t234 + (-t269 * t100 - t218 + t31 + t95) * t233; t101 * t47 - t118 * t67 + 0.2e1 * (t2 / 0.2e1 + t19 * t232 + t11 * t231) * m(7) + t40;];
tauc  = t16(:);
