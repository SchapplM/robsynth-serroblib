% Calculate vector of inverse dynamics joint torques for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR9_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR9_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR9_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:21
% EndTime: 2019-12-31 17:09:37
% DurationCPUTime: 8.40s
% Computational Cost: add. (2324->422), mult. (5418->603), div. (0->0), fcn. (3561->10), ass. (0->208)
t151 = cos(qJ(2));
t207 = qJD(1) * t151;
t130 = qJD(4) - t207;
t237 = t130 / 0.2e1;
t145 = cos(pkin(7));
t144 = sin(pkin(7));
t148 = sin(qJ(2));
t208 = qJD(1) * t148;
t197 = t144 * t208;
t108 = qJD(2) * t145 - t197;
t195 = t145 * t208;
t109 = qJD(2) * t144 + t195;
t147 = sin(qJ(4));
t150 = cos(qJ(4));
t56 = t108 * t147 + t109 * t150;
t243 = t56 / 0.2e1;
t183 = t150 * t108 - t109 * t147;
t245 = t183 / 0.2e1;
t281 = Ifges(5,5) * t243 + Ifges(5,6) * t245 + Ifges(5,3) * t237;
t189 = -t207 / 0.2e1;
t112 = t144 * t150 + t145 * t147;
t165 = t112 * t151;
t76 = qJD(1) * t165;
t95 = t112 * qJD(4);
t280 = t76 - t95;
t170 = t144 * t147 - t145 * t150;
t164 = t170 * t151;
t77 = qJD(1) * t164;
t94 = t170 * qJD(4);
t279 = t77 - t94;
t227 = Ifges(3,4) * t148;
t176 = t151 * Ifges(3,2) + t227;
t186 = -qJ(3) * t148 - pkin(1);
t119 = -pkin(2) * t151 + t186;
t100 = t119 * qJD(1);
t136 = pkin(5) * t207;
t123 = qJD(2) * qJ(3) + t136;
t59 = t145 * t100 - t123 * t144;
t31 = -pkin(3) * t207 - pkin(6) * t109 + t59;
t60 = t144 * t100 + t145 * t123;
t32 = pkin(6) * t108 + t60;
t8 = -t147 * t32 + t150 * t31;
t9 = t147 * t31 + t150 * t32;
t278 = t8 * mrSges(5,1) - Ifges(3,6) * qJD(2) / 0.2e1 - qJD(1) * t176 / 0.2e1 - t9 * mrSges(5,2) + t109 * Ifges(4,5) / 0.2e1 + t108 * Ifges(4,6) / 0.2e1 + Ifges(4,3) * t189 + t281;
t235 = Ifges(5,4) * t56;
t20 = Ifges(5,2) * t183 + Ifges(5,6) * t130 + t235;
t249 = t20 / 0.2e1;
t52 = Ifges(5,4) * t183;
t21 = Ifges(5,1) * t56 + Ifges(5,5) * t130 + t52;
t248 = t21 / 0.2e1;
t203 = qJD(1) * qJD(2);
t117 = qJDD(1) * t148 + t151 * t203;
t149 = sin(qJ(1));
t152 = cos(qJ(1));
t277 = g(1) * t152 + g(2) * t149;
t139 = t151 * qJDD(1);
t192 = t148 * t203;
t116 = -t139 + t192;
t239 = t116 / 0.2e1;
t79 = qJDD(2) * t144 + t117 * t145;
t241 = t79 / 0.2e1;
t276 = Ifges(4,1) * t241 + Ifges(4,5) * t239;
t78 = qJDD(2) * t145 - t117 * t144;
t17 = qJD(4) * t183 + t147 * t78 + t150 * t79;
t251 = t17 / 0.2e1;
t18 = -qJD(4) * t56 - t147 * t79 + t150 * t78;
t250 = t18 / 0.2e1;
t242 = t78 / 0.2e1;
t275 = m(4) + m(5);
t110 = qJDD(4) + t116;
t240 = t110 / 0.2e1;
t274 = -t116 / 0.2e1;
t273 = t117 / 0.2e1;
t272 = qJD(2) / 0.2e1;
t212 = t145 * t151;
t169 = pkin(3) * t148 - pkin(6) * t212;
t172 = pkin(2) * t148 - qJ(3) * t151;
t114 = t172 * qJD(1);
t68 = pkin(5) * t197 + t145 * t114;
t44 = qJD(1) * t169 + t68;
t213 = t145 * t148;
t214 = t144 * t151;
t166 = -pkin(5) * t213 - pkin(6) * t214;
t96 = t144 * t114;
t57 = qJD(1) * t166 + t96;
t228 = pkin(6) + qJ(3);
t120 = t228 * t144;
t121 = t228 * t145;
t65 = -t120 * t150 - t121 * t147;
t271 = -qJD(3) * t170 + qJD(4) * t65 - t147 * t44 - t150 * t57;
t66 = -t120 * t147 + t121 * t150;
t270 = -qJD(3) * t112 - qJD(4) * t66 + t147 * t57 - t150 * t44;
t266 = -qJD(2) * mrSges(3,1) - mrSges(4,1) * t108 + mrSges(4,2) * t109 + mrSges(3,3) * t208;
t179 = -mrSges(4,1) * t145 + mrSges(4,2) * t144;
t163 = m(4) * pkin(2) - t179;
t265 = t151 * t163;
t224 = Ifges(4,4) * t145;
t175 = -Ifges(4,2) * t144 + t224;
t225 = Ifges(4,4) * t144;
t177 = Ifges(4,1) * t145 - t225;
t264 = t108 * (Ifges(4,6) * t148 + t151 * t175) + t109 * (Ifges(4,5) * t148 + t151 * t177);
t133 = pkin(5) * t139;
t106 = -pkin(5) * t192 + t133;
t107 = t117 * pkin(5);
t263 = t106 * t151 + t107 * t148;
t204 = qJD(3) * t148;
t219 = qJDD(1) * pkin(1);
t48 = pkin(2) * t116 - qJ(3) * t117 - qJD(1) * t204 - t219;
t135 = pkin(5) * t208;
t82 = qJDD(2) * qJ(3) + t133 + (qJD(3) - t135) * qJD(2);
t27 = -t144 * t82 + t145 * t48;
t28 = t144 * t48 + t145 * t82;
t262 = -t144 * t27 + t145 * t28;
t261 = -m(3) - t275;
t181 = mrSges(3,1) * t151 - mrSges(3,2) * t148;
t260 = t148 * mrSges(5,3) + mrSges(2,1) + t181;
t118 = -qJD(2) * pkin(2) + qJD(3) + t135;
t178 = t144 * mrSges(4,1) + t145 * mrSges(4,2);
t259 = -t118 * t151 * t178 - t60 * (-mrSges(4,2) * t148 - mrSges(4,3) * t214) - t59 * (mrSges(4,1) * t148 - mrSges(4,3) * t212);
t12 = pkin(3) * t116 - pkin(6) * t79 + t27;
t13 = pkin(6) * t78 + t28;
t1 = qJD(4) * t8 + t12 * t147 + t13 * t150;
t2 = -qJD(4) * t9 + t12 * t150 - t13 * t147;
t257 = -t2 * mrSges(5,1) + t1 * mrSges(5,2);
t234 = pkin(3) * t144;
t256 = -m(5) * t234 + mrSges(2,2) - mrSges(3,3) - t178;
t253 = Ifges(5,4) * t251 + Ifges(5,2) * t250 + Ifges(5,6) * t240;
t252 = Ifges(5,1) * t251 + Ifges(5,4) * t250 + Ifges(5,5) * t240;
t247 = Ifges(4,4) * t242 + t276;
t246 = -t183 / 0.2e1;
t244 = -t56 / 0.2e1;
t238 = -t130 / 0.2e1;
t233 = pkin(5) * t148;
t226 = Ifges(3,4) * t151;
t91 = -qJDD(2) * pkin(2) + qJDD(3) + t107;
t220 = t148 * t91;
t206 = qJD(2) * t148;
t200 = pkin(5) * t206;
t92 = qJD(2) * t172 - t204;
t62 = t144 * t200 + t145 * t92;
t215 = t144 * t148;
t211 = t149 * t151;
t210 = t151 * t152;
t73 = pkin(5) * t212 + t144 * t119;
t205 = qJD(2) * t151;
t201 = Ifges(5,5) * t17 + Ifges(5,6) * t18 + Ifges(5,3) * t110;
t198 = pkin(5) + t234;
t196 = t144 * t207;
t194 = m(4) * qJ(3) + mrSges(4,3);
t193 = m(5) * t228 + mrSges(5,3);
t33 = -t78 * mrSges(4,1) + t79 * mrSges(4,2);
t5 = -t18 * mrSges(5,1) + t17 * mrSges(5,2);
t185 = -t203 / 0.2e1;
t180 = mrSges(3,1) * t148 + mrSges(3,2) * t151;
t174 = Ifges(3,5) * t151 - Ifges(3,6) * t148;
t173 = Ifges(4,5) * t145 - Ifges(4,6) * t144;
t105 = t145 * t119;
t58 = -pkin(6) * t213 + t105 + (-pkin(5) * t144 - pkin(3)) * t151;
t64 = -pkin(6) * t215 + t73;
t25 = -t147 * t64 + t150 * t58;
t26 = t147 * t58 + t150 * t64;
t131 = pkin(3) * t145 + pkin(2);
t171 = t131 * t151 + t148 * t228;
t168 = pkin(1) * t180;
t167 = t148 * (Ifges(3,1) * t151 - t227);
t143 = pkin(7) + qJ(4);
t137 = sin(t143);
t138 = cos(t143);
t162 = m(5) * t131 + t138 * mrSges(5,1) - t137 * mrSges(5,2);
t156 = t151 * (Ifges(4,3) * t148 + t151 * t173);
t155 = t148 * t194 + t265;
t134 = Ifges(3,4) * t207;
t124 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t207;
t115 = t198 * t148;
t103 = t198 * t205;
t102 = pkin(3) * t196 + t136;
t99 = Ifges(3,1) * t208 + Ifges(3,5) * qJD(2) + t134;
t90 = t170 * t148;
t89 = t112 * t148;
t86 = t149 * t137 + t138 * t210;
t85 = -t137 * t210 + t149 * t138;
t84 = t137 * t152 - t138 * t211;
t83 = t137 * t211 + t138 * t152;
t80 = t144 * t92;
t75 = -mrSges(4,1) * t207 - mrSges(4,3) * t109;
t74 = mrSges(4,2) * t207 + mrSges(4,3) * t108;
t72 = -pkin(5) * t214 + t105;
t69 = -pkin(5) * t195 + t96;
t67 = -pkin(3) * t108 + t118;
t63 = -t145 * t200 + t80;
t51 = t109 * Ifges(4,1) + t108 * Ifges(4,4) - Ifges(4,5) * t207;
t50 = t109 * Ifges(4,4) + t108 * Ifges(4,2) - Ifges(4,6) * t207;
t47 = mrSges(4,1) * t116 - mrSges(4,3) * t79;
t46 = -mrSges(4,2) * t116 + mrSges(4,3) * t78;
t45 = qJD(2) * t166 + t80;
t41 = -pkin(3) * t78 + t91;
t40 = -qJD(2) * t165 + t148 * t94;
t39 = -qJD(2) * t164 - t148 * t95;
t38 = qJD(2) * t169 + t62;
t37 = mrSges(5,1) * t130 - mrSges(5,3) * t56;
t36 = -mrSges(5,2) * t130 + mrSges(5,3) * t183;
t29 = t79 * Ifges(4,4) + t78 * Ifges(4,2) + t116 * Ifges(4,6);
t24 = -mrSges(5,1) * t183 + mrSges(5,2) * t56;
t11 = -mrSges(5,2) * t110 + mrSges(5,3) * t18;
t10 = mrSges(5,1) * t110 - mrSges(5,3) * t17;
t7 = -qJD(4) * t26 - t147 * t45 + t150 * t38;
t6 = qJD(4) * t25 + t147 * t38 + t150 * t45;
t3 = [(-t213 * t27 - t215 * t28) * mrSges(4,3) + t266 * pkin(5) * t205 + t181 * t219 + (-t86 * mrSges(5,1) - t85 * mrSges(5,2) + t261 * (t152 * pkin(1) + t149 * pkin(5)) + t256 * t149 + (-m(5) * t171 - t155 - t260) * t152) * g(2) + t115 * t5 - pkin(1) * (mrSges(3,1) * t116 + mrSges(3,2) * t117) + t103 * t24 + t67 * (-mrSges(5,1) * t40 + mrSges(5,2) * t39) + t72 * t47 + t73 * t46 + t63 * t74 + t62 * t75 + t213 * t247 + t39 * t248 + t40 * t249 - t90 * t252 - t89 * t253 + t33 * t233 + t178 * t220 + t156 * t185 + (-mrSges(3,1) * t233 + Ifges(3,5) * t148 + (-mrSges(3,2) * pkin(5) + Ifges(3,6)) * t151) * qJDD(2) + (-pkin(5) * t116 * t151 + t117 * t233 + t263) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(5) * t263) + t6 * t36 + t7 * t37 + t25 * t10 + t26 * t11 + (-t84 * mrSges(5,1) - t83 * mrSges(5,2) + (m(3) * pkin(1) - m(4) * t186 + t148 * mrSges(4,3) + t265 - m(5) * (-pkin(1) - t171) + t260) * t149 + (t261 * pkin(5) + t256) * t152) * g(1) + (Ifges(5,1) * t39 + Ifges(5,4) * t40) * t243 + (t174 * t272 - t259) * qJD(2) - (Ifges(4,5) * t79 + Ifges(4,6) * t78 + Ifges(4,3) * t116 + t201) * t151 / 0.2e1 + (-t1 * t89 + t2 * t90 - t39 * t8 + t40 * t9) * mrSges(5,3) + t41 * (mrSges(5,1) * t89 - mrSges(5,2) * t90) + (-Ifges(5,1) * t90 - Ifges(5,4) * t89) * t251 + (-Ifges(5,4) * t90 - Ifges(5,2) * t89) * t250 + (-Ifges(5,5) * t90 - Ifges(5,6) * t89) * t240 + (t167 + t151 * (-Ifges(3,2) * t148 + t226)) * t203 / 0.2e1 + (t145 * t51 + t99) * t205 / 0.2e1 + t264 * t272 + t226 * t273 + t176 * t274 - t144 * t50 * t205 / 0.2e1 + (Ifges(5,4) * t39 + Ifges(5,2) * t40) * t245 + (t278 + t281) * t206 + m(5) * (t1 * t26 + t103 * t67 + t115 * t41 + t2 * t25 + t6 * t9 + t7 * t8) + (Ifges(5,5) * t39 + Ifges(5,6) * t40) * t237 + (-t27 * mrSges(4,1) + t28 * mrSges(4,2) + Ifges(3,4) * t273 - Ifges(4,5) * t241 - Ifges(5,5) * t251 + Ifges(3,2) * t274 - Ifges(4,6) * t242 - Ifges(5,6) * t250 - Ifges(4,3) * t239 - Ifges(5,3) * t240 + t257) * t151 + (Ifges(3,1) * t117 + Ifges(3,4) * t274 + t173 * t239 + t175 * t242 + t177 * t241) * t148 - t124 * t200 - t168 * t203 - t29 * t215 / 0.2e1 + m(4) * (t27 * t72 + t28 * t73 + t59 * t62 + t60 * t63 + (t118 * t205 + t220) * pkin(5)) + Ifges(2,3) * qJDD(1); t277 * (t180 + (-t193 - t194) * t151 + (t162 + t163) * t148) + (-Ifges(5,5) * t77 - Ifges(5,6) * t76) * t238 + (-Ifges(5,1) * t77 - Ifges(5,4) * t76) * t244 + (-Ifges(5,4) * t77 - Ifges(5,2) * t76) * t246 - t170 * t253 + t41 * (mrSges(5,1) * t170 + mrSges(5,2) * t112) + (Ifges(5,5) * t112 - Ifges(5,6) * t170) * t240 + (Ifges(5,4) * t112 - Ifges(5,2) * t170) * t250 + (Ifges(5,1) * t112 - Ifges(5,4) * t170) * t251 + (-t148 * t193 - t151 * t162 - t155 - t181) * g(3) + t225 * t242 + (t134 + t99) * t189 + (Ifges(5,5) * t244 - Ifges(3,2) * t189 + Ifges(5,6) * t246 + Ifges(5,3) * t238 - t278) * t208 - t131 * t5 - Ifges(3,6) * t116 + Ifges(3,5) * t117 - t102 * t24 - t106 * mrSges(3,2) - t107 * mrSges(3,1) + t66 * t11 - t69 * t74 - t68 * t75 + t112 * t252 + t174 * t185 + t224 * t241 + t279 * t248 + t280 * t249 + (-pkin(2) * t91 + qJ(3) * t262 - t118 * t136 - t59 * t68 - t60 * t69) * m(4) + t262 * mrSges(4,3) + t65 * t10 - pkin(2) * t33 + (t259 - t264 / 0.2e1 + (t168 + t156 / 0.2e1 - t167 / 0.2e1) * qJD(1)) * qJD(1) + (-Ifges(5,1) * t94 - Ifges(5,4) * t95) * t243 + (-Ifges(5,4) * t94 - Ifges(5,2) * t95) * t245 + (-Ifges(5,5) * t94 - Ifges(5,6) * t95) * t237 - t266 * t136 + t124 * t135 + t270 * t37 + (t1 * t66 - t102 * t67 - t131 * t41 + t2 * t65 + t270 * t8 + t271 * t9) * m(5) + t271 * t36 + t91 * t179 + (qJ(3) * t46 + t29 / 0.2e1 + Ifges(4,2) * t242 + Ifges(4,6) * t239 + t51 * t189 + (m(4) * t60 + t74) * qJD(3)) * t145 + (-t1 * t170 - t112 * t2 - t279 * t8 + t280 * t9) * mrSges(5,3) + (-mrSges(5,1) * t280 + mrSges(5,2) * t279) * t67 + t50 * t196 / 0.2e1 + (-qJ(3) * t47 + t247 + (-m(4) * t59 - t75) * qJD(3) + t276) * t144 + Ifges(3,3) * qJDD(2); -t108 * t74 + t109 * t75 - t183 * t36 + t56 * t37 + t33 + t5 + (t151 * g(3) - t148 * t277) * t275 + (-t183 * t9 + t56 * t8 + t41) * m(5) + (-t108 * t60 + t109 * t59 + t91) * m(4); -t67 * (mrSges(5,1) * t56 + mrSges(5,2) * t183) + (Ifges(5,1) * t183 - t235) * t244 + t20 * t243 + (Ifges(5,5) * t183 - Ifges(5,6) * t56) * t238 - t8 * t36 + t9 * t37 - g(1) * (mrSges(5,1) * t85 - mrSges(5,2) * t86) - g(2) * (-mrSges(5,1) * t83 + mrSges(5,2) * t84) - g(3) * (-mrSges(5,1) * t137 - mrSges(5,2) * t138) * t148 + (t183 * t8 + t56 * t9) * mrSges(5,3) + t201 + (-Ifges(5,2) * t56 + t21 + t52) * t246 - t257;];
tau = t3;
