% Calculate vector of inverse dynamics joint torques for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:27
% EndTime: 2019-12-05 16:02:43
% DurationCPUTime: 6.39s
% Computational Cost: add. (2052->395), mult. (4456->568), div. (0->0), fcn. (2991->10), ass. (0->196)
t221 = -m(5) - m(4);
t255 = m(5) + m(6);
t101 = sin(qJ(5));
t104 = cos(qJ(5));
t179 = qJD(4) * t104;
t105 = cos(qJ(4));
t183 = qJD(2) * t105;
t75 = -t101 * t183 + t179;
t254 = -t75 / 0.2e1;
t181 = qJD(4) * t101;
t76 = t104 * t183 + t181;
t253 = -t76 / 0.2e1;
t102 = sin(qJ(4));
t185 = qJD(2) * t102;
t95 = qJD(5) + t185;
t252 = -t95 / 0.2e1;
t138 = t101 * mrSges(6,1) + t104 * mrSges(6,2);
t215 = mrSges(3,1) - mrSges(4,2);
t251 = -t138 - mrSges(5,3) - t215;
t235 = m(4) + t255;
t250 = pkin(2) * t235 + t255 * pkin(7) - t251;
t139 = -mrSges(6,1) * t104 + mrSges(6,2) * t101;
t119 = m(6) * pkin(4) - t139;
t140 = mrSges(5,1) * t102 + mrSges(5,2) * t105;
t244 = mrSges(3,2) - mrSges(4,3);
t249 = -t119 * t102 - t140 + t244;
t172 = qJD(2) * qJD(4);
t80 = qJDD(2) * t105 - t102 * t172;
t26 = qJD(5) * t75 + qJDD(4) * t101 + t104 * t80;
t226 = t26 / 0.2e1;
t27 = -qJD(5) * t76 + qJDD(4) * t104 - t101 * t80;
t225 = t27 / 0.2e1;
t81 = -qJDD(2) * t102 - t105 * t172;
t72 = qJDD(5) - t81;
t224 = t72 / 0.2e1;
t248 = t80 / 0.2e1;
t247 = t81 / 0.2e1;
t5 = -mrSges(6,1) * t27 + mrSges(6,2) * t26;
t245 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t80 - t5;
t163 = mrSges(5,3) * t183;
t213 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t75 + mrSges(6,2) * t76 + t163;
t106 = cos(qJ(2));
t107 = -pkin(2) - pkin(7);
t177 = qJD(4) * t107;
t152 = t105 * t177;
t103 = sin(qJ(2));
t192 = t101 * t103;
t98 = sin(pkin(5));
t196 = qJD(1) * t98;
t190 = t102 * t107;
t147 = -pkin(8) * t105 + qJ(3);
t82 = pkin(4) * t102 + t147;
t51 = t101 * t82 + t104 * t190;
t143 = pkin(4) * t105 + pkin(8) * t102;
t73 = qJD(4) * t143 + qJD(3);
t243 = -qJD(5) * t51 - t101 * t152 + t104 * t73 - (-t102 * t192 + t104 * t106) * t196;
t189 = t103 * t104;
t50 = -t101 * t190 + t104 * t82;
t242 = qJD(5) * t50 + t101 * t73 + t104 * t152 - (t101 * t106 + t102 * t189) * t196;
t241 = mrSges(5,1) + t119;
t165 = m(6) * pkin(8) + mrSges(6,3);
t240 = mrSges(5,2) - t165;
t210 = Ifges(5,4) * t102;
t137 = t105 * Ifges(5,1) - t210;
t71 = Ifges(6,4) * t75;
t23 = Ifges(6,1) * t76 + Ifges(6,5) * t95 + t71;
t239 = Ifges(5,5) * qJD(4) + qJD(2) * t137 + t104 * t23;
t100 = cos(pkin(5));
t171 = qJDD(1) * t100;
t186 = qJD(1) * t100;
t161 = t106 * t196;
t141 = qJD(3) - t161;
t67 = qJD(2) * t107 + t141;
t42 = t102 * t67 + t105 * t186;
t195 = qJD(4) * t42;
t197 = t106 * t98;
t184 = qJD(2) * t103;
t160 = t98 * t184;
t89 = qJD(1) * t160;
t57 = qJDD(1) * t197 - t89;
t125 = qJDD(3) - t57;
t43 = qJDD(2) * t107 + t125;
t11 = -t102 * t171 + t105 * t43 - t195;
t29 = qJD(4) * pkin(8) + t42;
t162 = t103 * t196;
t49 = qJD(2) * t82 + t162;
t12 = -t101 * t29 + t104 * t49;
t13 = t101 * t49 + t104 * t29;
t175 = qJD(5) * t104;
t176 = qJD(5) * t101;
t238 = -t12 * t175 - t13 * t176;
t14 = mrSges(6,1) * t72 - mrSges(6,3) * t26;
t15 = -mrSges(6,2) * t72 + mrSges(6,3) * t27;
t237 = -t101 * t14 + t104 * t15;
t155 = t102 * t186;
t178 = qJD(4) * t105;
t10 = -qJD(4) * t155 + t102 * t43 + t105 * t171 + t67 * t178;
t236 = t10 * t102 + t105 * t11;
t209 = Ifges(5,4) * t105;
t134 = -Ifges(5,2) * t102 + t209;
t233 = Ifges(5,6) * qJD(4) / 0.2e1 + qJD(2) * t134 / 0.2e1 + Ifges(6,5) * t253 + Ifges(6,6) * t254 + Ifges(6,3) * t252;
t219 = Ifges(6,4) * t76;
t22 = Ifges(6,2) * t75 + Ifges(6,6) * t95 + t219;
t41 = t105 * t67 - t155;
t28 = -qJD(4) * pkin(4) - t41;
t231 = -t101 * t22 / 0.2e1 + t138 * t28;
t169 = qJDD(2) * qJ(3);
t170 = qJDD(1) * t103;
t44 = t98 * t170 + t169 + (qJD(3) + t161) * qJD(2);
t18 = -pkin(4) * t81 - pkin(8) * t80 + t44;
t8 = qJDD(4) * pkin(8) + t10;
t1 = qJD(5) * t12 + t101 * t18 + t104 * t8;
t2 = -qJD(5) * t13 - t101 * t8 + t104 * t18;
t230 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t228 = -m(6) * t147 + t105 * mrSges(6,3) + t221 * qJ(3) + t249;
t108 = qJD(2) ^ 2;
t227 = Ifges(6,1) * t226 + Ifges(6,4) * t225 + Ifges(6,5) * t224;
t222 = t76 / 0.2e1;
t9 = -qJDD(4) * pkin(4) - t11;
t218 = t105 * t9;
t212 = t98 * t103 * qJ(3) + pkin(2) * t197;
t211 = mrSges(6,3) * t101;
t208 = Ifges(6,4) * t101;
t207 = Ifges(6,4) * t104;
t204 = t102 * t98;
t203 = t104 * mrSges(6,3);
t199 = t105 * t98;
t79 = qJD(2) * qJ(3) + t162;
t198 = t106 * t79;
t194 = t100 * t103;
t193 = t100 * t106;
t191 = t101 * t105;
t188 = t104 * t105;
t182 = qJD(2) * t106;
t180 = qJD(4) * t102;
t174 = qJD(5) * t105;
t173 = qJDD(2) * mrSges(4,2);
t168 = Ifges(6,5) * t26 + Ifges(6,6) * t27 + Ifges(6,3) * t72;
t166 = t102 * t197;
t164 = mrSges(5,3) * t185;
t159 = t98 * t182;
t154 = t101 * t180;
t151 = qJD(1) * t182;
t146 = -t172 / 0.2e1;
t145 = t105 * t162;
t144 = t98 * t151;
t142 = t1 * t104 - t101 * t2;
t136 = Ifges(6,1) * t104 - t208;
t135 = Ifges(6,1) * t101 + t207;
t133 = -Ifges(6,2) * t101 + t207;
t132 = Ifges(6,2) * t104 + t208;
t131 = -Ifges(5,5) * t102 - Ifges(5,6) * t105;
t130 = Ifges(6,5) * t104 - Ifges(6,6) * t101;
t129 = Ifges(6,5) * t101 + Ifges(6,6) * t104;
t128 = t101 * t13 + t104 * t12;
t47 = -mrSges(6,2) * t95 + mrSges(6,3) * t75;
t48 = mrSges(6,1) * t95 - mrSges(6,3) * t76;
t127 = -t101 * t47 - t104 * t48;
t126 = qJ(3) * t44 + qJD(3) * t79;
t66 = t100 * t105 - t166;
t37 = -t101 * t66 + t189 * t98;
t38 = t104 * t66 + t192 * t98;
t124 = t103 * t44 + t182 * t79;
t65 = t100 * t102 + t105 * t197;
t123 = t79 * (mrSges(5,1) * t105 - mrSges(5,2) * t102);
t121 = t102 * (-Ifges(5,2) * t105 - t210);
t120 = t105 * (-Ifges(5,1) * t102 - t209);
t116 = t101 * t174 + t102 * t179;
t115 = -t104 * t174 + t154;
t112 = Ifges(6,5) * t105 - t102 * t136;
t111 = Ifges(6,6) * t105 - t102 * t133;
t110 = Ifges(6,3) * t105 - t102 * t130;
t99 = cos(pkin(9));
t97 = sin(pkin(9));
t85 = -qJD(4) * mrSges(5,2) - t164;
t78 = t143 * qJD(2);
t77 = t140 * qJD(2);
t74 = -qJD(2) * pkin(2) + t141;
t64 = t106 * t99 - t194 * t97;
t63 = t103 * t99 + t193 * t97;
t62 = t106 * t97 + t194 * t99;
t61 = t103 * t97 - t193 * t99;
t60 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t81;
t58 = (t151 + t170) * t98;
t52 = -qJDD(2) * pkin(2) + t125;
t40 = -mrSges(5,1) * t81 + mrSges(5,2) * t80;
t36 = -qJD(4) * t166 + t100 * t178 - t105 * t160;
t35 = -qJD(4) * t65 + t102 * t160;
t34 = -t61 * t102 + t199 * t99;
t32 = t102 * t63 + t199 * t97;
t20 = t101 * t78 + t104 * t41;
t19 = -t101 * t41 + t104 * t78;
t7 = qJD(5) * t37 + t101 * t159 + t104 * t35;
t6 = -qJD(5) * t38 - t101 * t35 + t104 * t159;
t3 = t26 * Ifges(6,4) + t27 * Ifges(6,2) + t72 * Ifges(6,6);
t4 = [t37 * t14 + t38 * t15 + t35 * t85 + t7 * t47 + t6 * t48 + t66 * t60 - t245 * t65 + t213 * t36 + (-m(2) - m(3) - m(6) + t221) * g(3) + m(5) * (t10 * t66 - t11 * t65 + t35 * t42 - t36 * t41) + m(6) * (t1 * t38 + t12 * t6 + t13 * t7 + t2 * t37 + t28 * t36 + t65 * t9) + (t77 * t182 + t103 * t40 + m(4) * (-t106 * t52 + t184 * t74 + t124) + m(3) * (t103 * t58 + t106 * t57) + m(5) * t124 + (-t103 * t215 - t106 * t244) * t108 + (-t103 * t244 + t106 * t215) * qJDD(2)) * t98 + (m(2) + 0.2e1 * (m(4) / 0.2e1 + m(3) / 0.2e1) * t100 ^ 2) * qJDD(1); (qJD(2) * qJD(3) - t144 + t169 + t44) * mrSges(4,3) - (t101 * t23 + t104 * t22) * t174 / 0.2e1 + (-pkin(2) * t52 + t126 - (t103 * t74 + t198) * t196) * m(4) + t28 * (-mrSges(6,1) * t115 - mrSges(6,2) * t116) + (t52 - t89) * mrSges(4,2) + t134 * t247 + t137 * t248 + t138 * t218 + (qJD(4) * t112 - t135 * t174) * t222 + t188 * t227 + t60 * t190 + qJD(4) * t123 + (t57 + t89) * mrSges(3,1) + (-t1 * t191 + t115 * t13 + t116 * t12 - t188 * t2) * mrSges(6,3) + qJDD(4) * (Ifges(5,5) * t105 - Ifges(5,6) * t102) + t51 * t15 + t50 * t14 + qJ(3) * t40 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + t121 * t146 + (-t58 + t144) * mrSges(3,2) + (Ifges(6,3) * t224 + Ifges(6,6) * t225 + Ifges(6,5) * t226 - Ifges(5,4) * t80 / 0.2e1 - Ifges(5,2) * t81 / 0.2e1 - t85 * t162 + t168 / 0.2e1 + t230) * t102 + t245 * t105 * t107 + (Ifges(5,1) * t248 + Ifges(5,4) * t247 + t130 * t224 + t133 * t225 + t136 * t226) * t105 + t44 * t140 + qJD(4) ^ 2 * t131 / 0.2e1 + (t228 * t62 + t250 * t61) * g(2) + (t228 * t64 + t250 * t63) * g(1) + t85 * t152 + t243 * t48 + (t145 * t28 + t1 * t51 + t2 * t50 + (t180 * t28 - t218) * t107 + t242 * t13 + t243 * t12) * m(6) + t242 * t47 - t239 * t180 / 0.2e1 + t141 * t77 + (-(t198 + (t102 * t42 + t105 * t41) * t103) * t196 + ((-t102 * t41 + t105 * t42) * qJD(4) + t236) * t107 + t126) * m(5) + (-t178 * t42 + t180 * t41 - t236) * mrSges(5,3) + (t12 * mrSges(6,1) - t13 * mrSges(6,2) - t233) * t178 + (-m(4) * t212 - t255 * (pkin(7) * t197 + t212) + (t251 * t106 + (t105 * t165 + t249) * t103) * t98) * g(3) + (t102 * t177 + t145) * t213 + t22 * t154 / 0.2e1 + t120 * t172 / 0.2e1 - pkin(2) * t173 + t75 * (qJD(4) * t111 - t132 * t174) / 0.2e1 + t95 * (qJD(4) * t110 - t129 * t174) / 0.2e1 - t3 * t191 / 0.2e1; t173 - t108 * mrSges(4,3) + m(4) * t52 + ((-t101 * t48 + t104 * t47 + t85) * qJD(4) + m(6) * (-t12 * t181 + t13 * t179 - t9) + m(5) * (t11 + t195) + t245) * t105 + (t60 + t127 * qJD(5) + t213 * qJD(4) + m(6) * (qJD(4) * t28 + t142 + t238) + m(5) * (-qJD(4) * t41 + t10) + t237) * t102 + (-m(6) * t128 + t221 * t79 + t127 - t77) * qJD(2) + t235 * (-g(1) * t63 - g(2) * t61 + g(3) * t197); (-t120 / 0.2e1 + t121 / 0.2e1) * t108 + (t130 * t95 + t133 * t75 + t136 * t76) * qJD(5) / 0.2e1 - (t110 * t95 + t111 * t75 + t112 * t76) * qJD(2) / 0.2e1 + (-t123 - t12 * (mrSges(6,1) * t105 + t102 * t203) - t13 * (-mrSges(6,2) * t105 + t102 * t211)) * qJD(2) + t101 * t227 + t129 * t224 + t132 * t225 + t135 * t226 + t1 * t203 + (t231 + t239 / 0.2e1) * t185 + t104 * t3 / 0.2e1 + (-t85 - t164) * t41 + Ifges(5,5) * t80 + Ifges(5,6) * t81 - t20 * t47 - t19 * t48 - pkin(4) * t5 - t10 * mrSges(5,2) + t11 * mrSges(5,1) + t131 * t146 + t231 * qJD(5) + t9 * t139 + (-m(6) * t28 + t163 - t213) * t42 + (t240 * t66 + t241 * t65) * g(3) + (-t240 * t34 - t241 * (t105 * t61 + t204 * t99)) * g(2) + (t240 * t32 - t241 * (t105 * t63 - t204 * t97)) * g(1) + (m(6) * (-qJD(5) * t128 + t142) - t48 * t175 - t47 * t176 + t237) * pkin(8) + t238 * mrSges(6,3) + t233 * t183 + (-pkin(4) * t9 - t12 * t19 - t13 * t20) * m(6) + t23 * t175 / 0.2e1 + Ifges(5,3) * qJDD(4) - t2 * t211; -t28 * (mrSges(6,1) * t76 + mrSges(6,2) * t75) + (Ifges(6,1) * t75 - t219) * t253 + t22 * t222 + (Ifges(6,5) * t75 - Ifges(6,6) * t76) * t252 - t12 * t47 + t13 * t48 - g(1) * ((-t101 * t32 + t104 * t64) * mrSges(6,1) + (-t101 * t64 - t104 * t32) * mrSges(6,2)) - g(2) * ((t101 * t34 + t104 * t62) * mrSges(6,1) + (-t101 * t62 + t104 * t34) * mrSges(6,2)) - g(3) * (mrSges(6,1) * t37 - mrSges(6,2) * t38) + (t12 * t75 + t13 * t76) * mrSges(6,3) + t168 + (-Ifges(6,2) * t76 + t23 + t71) * t254 + t230;];
tau = t4;
