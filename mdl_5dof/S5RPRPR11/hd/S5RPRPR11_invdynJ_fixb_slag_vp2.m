% Calculate vector of inverse dynamics joint torques for
% S5RPRPR11
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR11_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:11
% EndTime: 2019-12-31 18:27:28
% DurationCPUTime: 9.04s
% Computational Cost: add. (3099->414), mult. (7454->520), div. (0->0), fcn. (5223->10), ass. (0->183)
t259 = mrSges(4,1) + mrSges(5,1);
t254 = Ifges(4,1) + Ifges(5,1);
t242 = Ifges(5,4) + Ifges(4,5);
t143 = sin(pkin(8));
t144 = cos(pkin(8));
t193 = t143 ^ 2 + t144 ^ 2;
t190 = qJD(1) * qJD(2);
t124 = qJDD(1) * qJ(2) + t190;
t206 = pkin(6) + qJ(2);
t116 = t206 * t143;
t109 = qJD(1) * t116;
t117 = t206 * t144;
t110 = qJD(1) * t117;
t147 = sin(qJ(3));
t218 = cos(qJ(3));
t69 = -t218 * t109 - t147 * t110;
t232 = -t69 + qJD(4);
t258 = Ifges(4,4) - Ifges(5,5);
t241 = Ifges(5,6) - Ifges(4,6);
t108 = t143 * t218 + t147 * t144;
t100 = t108 * qJD(1);
t257 = -pkin(7) * t100 + t232;
t204 = mrSges(4,3) * t100;
t205 = mrSges(5,2) * t100;
t239 = t259 * qJD(3) - t204 - t205;
t57 = -qJD(3) * pkin(3) + t232;
t256 = -m(5) * t57 + t239;
t184 = t218 * t144;
t157 = -t147 * t143 + t184;
t101 = t157 * qJD(3);
t65 = qJD(1) * t101 + qJDD(1) * t108;
t227 = t65 / 0.2e1;
t99 = t157 * qJD(1);
t215 = Ifges(5,5) * t99;
t89 = Ifges(4,4) * t99;
t252 = t242 * qJD(3) + t254 * t100 - t215 + t89;
t80 = mrSges(5,2) * t99 + qJD(3) * mrSges(5,3);
t240 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t99 + t80;
t199 = qJDD(3) / 0.2e1;
t131 = pkin(2) * t144 + pkin(1);
t111 = -qJDD(1) * t131 + qJDD(2);
t251 = -qJ(4) * t65 - qJD(4) * t100 + t111;
t112 = -qJD(1) * t131 + qJD(2);
t250 = -qJ(4) * t100 + t112;
t138 = pkin(8) + qJ(3);
t133 = sin(t138);
t134 = cos(t138);
t249 = -t259 * t134 + (mrSges(4,2) - mrSges(5,3)) * t133;
t167 = -mrSges(3,1) * t144 + mrSges(3,2) * t143;
t245 = -m(3) * pkin(1) - mrSges(2,1) + t167 + t249;
t181 = m(3) * qJ(2) + mrSges(3,3);
t244 = mrSges(2,2) - m(6) * (-pkin(7) + t206) + mrSges(6,3) - mrSges(5,2) - t181 - mrSges(4,3);
t41 = -pkin(3) * t99 + t250;
t243 = mrSges(4,2) * t112 + t57 * mrSges(5,2) - t69 * mrSges(4,3) - mrSges(5,3) * t41;
t146 = sin(qJ(5));
t149 = cos(qJ(5));
t151 = -pkin(3) - pkin(4);
t113 = -qJ(4) * t146 + t149 * t151;
t70 = -t147 * t109 + t218 * t110;
t38 = -pkin(7) * t99 + t70;
t238 = qJD(5) * t113 - t146 * t38 + t149 * t257;
t114 = qJ(4) * t149 + t146 * t151;
t237 = -qJD(5) * t114 - t146 * t257 - t149 * t38;
t148 = sin(qJ(1));
t150 = cos(qJ(1));
t233 = g(1) * t150 + g(2) * t148;
t236 = t233 * t133;
t72 = -t147 * t116 + t218 * t117;
t196 = t134 * t150;
t197 = t133 * t150;
t235 = pkin(3) * t196 + qJ(4) * t197;
t139 = -qJD(3) + qJD(5);
t178 = qJD(3) * t218;
t192 = qJD(3) * t147;
t177 = pkin(6) * qJDD(1) + t124;
t86 = t177 * t143;
t87 = t177 * t144;
t25 = t109 * t192 - t110 * t178 - t147 * t87 - t218 * t86;
t154 = qJDD(4) - t25;
t22 = -qJDD(3) * pkin(3) + t154;
t231 = -t100 * t41 - t22;
t55 = t100 * t146 + t149 * t99;
t230 = t55 / 0.2e1;
t161 = t100 * t149 - t146 * t99;
t229 = -t161 / 0.2e1;
t228 = t161 / 0.2e1;
t226 = t99 / 0.2e1;
t225 = -t99 / 0.2e1;
t224 = m(5) + m(6);
t222 = t100 / 0.2e1;
t219 = -t139 / 0.2e1;
t216 = Ifges(6,4) * t161;
t32 = qJD(3) * t151 + t257;
t142 = qJD(3) * qJ(4);
t33 = t142 + t38;
t10 = -t146 * t33 + t149 * t32;
t212 = t10 * mrSges(6,3);
t11 = t146 * t32 + t149 * t33;
t211 = t11 * mrSges(6,3);
t128 = t134 * pkin(3);
t207 = qJD(3) / 0.2e1;
t203 = Ifges(4,4) * t100;
t201 = qJ(4) * t99;
t127 = t133 * qJ(4);
t194 = t128 + t127;
t188 = qJDD(1) * t143;
t187 = qJDD(1) * t144;
t185 = t133 * t151;
t102 = t108 * qJD(3);
t66 = qJD(1) * t102 - qJDD(1) * t184 + t147 * t188;
t180 = t66 * mrSges(4,1) + t65 * mrSges(4,2);
t179 = t66 * mrSges(5,1) - t65 * mrSges(5,3);
t47 = -qJDD(3) * mrSges(5,1) + t65 * mrSges(5,2);
t176 = -t131 - t127;
t71 = t218 * t116 + t117 * t147;
t121 = t150 * t131;
t175 = t148 * t206 + t121;
t12 = -qJD(5) * t55 + t146 * t66 + t149 * t65;
t13 = -qJD(5) * t161 - t146 * t65 + t149 * t66;
t173 = -t13 * mrSges(6,1) + t12 * mrSges(6,2);
t94 = t133 * t149 - t134 * t146;
t73 = t94 * t148;
t160 = t133 * t146 + t134 * t149;
t74 = t160 * t148;
t172 = t73 * mrSges(6,1) - t74 * mrSges(6,2);
t75 = t146 * t196 - t149 * t197;
t76 = t160 * t150;
t171 = -t75 * mrSges(6,1) - t76 * mrSges(6,2);
t170 = -mrSges(6,1) * t160 - t94 * mrSges(6,2);
t168 = -mrSges(3,1) * t187 + mrSges(3,2) * t188;
t39 = -mrSges(6,2) * t139 - mrSges(6,3) * t55;
t40 = mrSges(6,1) * t139 - mrSges(6,3) * t161;
t162 = -t146 * t40 + t149 * t39;
t44 = -pkin(7) * t108 + t71;
t45 = -pkin(7) * t157 + t72;
t17 = -t146 * t45 + t149 * t44;
t18 = t146 * t44 + t149 * t45;
t67 = -t108 * t146 - t149 * t157;
t68 = t108 * t149 - t146 * t157;
t159 = qJ(4) * t101 + qJD(4) * t108;
t158 = qJ(4) * t108 + t131;
t24 = -t109 * t178 - t110 * t192 - t147 * t86 + t218 * t87;
t155 = t134 * mrSges(5,3) + (-m(5) * pkin(3) - mrSges(5,1)) * t133;
t21 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t24;
t42 = -t116 * t178 + qJD(2) * t184 + (-qJD(2) * t143 - qJD(3) * t117) * t147;
t8 = -t65 * pkin(7) + qJDD(3) * t151 + t154;
t9 = pkin(7) * t66 + t21;
t1 = qJD(5) * t10 + t146 * t8 + t149 * t9;
t135 = -qJDD(3) + qJDD(5);
t2 = -qJD(5) * t11 - t146 * t9 + t149 * t8;
t153 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t12 + Ifges(6,6) * t13 + Ifges(6,3) * t135;
t43 = t108 * qJD(2) + qJD(3) * t72;
t132 = -qJDD(1) * pkin(1) + qJDD(2);
t120 = qJ(4) * t196;
t118 = t148 * t134 * qJ(4);
t88 = Ifges(5,5) * t100;
t64 = -pkin(3) * t157 - t158;
t60 = -mrSges(5,1) * t99 - mrSges(5,3) * t100;
t59 = pkin(3) * t100 - t201;
t58 = t142 + t70;
t50 = t99 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t203;
t49 = Ifges(5,6) * qJD(3) - t99 * Ifges(5,3) + t88;
t48 = -mrSges(5,2) * t66 + qJDD(3) * mrSges(5,3);
t46 = Ifges(6,4) * t55;
t36 = -t151 * t157 + t158;
t35 = pkin(3) * t102 - t159;
t34 = t100 * t151 + t201;
t31 = -t101 * pkin(7) + t43;
t30 = pkin(7) * t102 + t42;
t29 = -t151 * t99 - t250;
t28 = t102 * t151 + t159;
t27 = -qJD(5) * t68 - t101 * t146 + t102 * t149;
t26 = qJD(5) * t67 + t101 * t149 + t102 * t146;
t23 = mrSges(6,1) * t55 + mrSges(6,2) * t161;
t20 = Ifges(6,1) * t161 + Ifges(6,5) * t139 - t46;
t19 = -Ifges(6,2) * t55 + Ifges(6,6) * t139 + t216;
t16 = pkin(3) * t66 + t251;
t7 = -mrSges(6,2) * t135 + mrSges(6,3) * t13;
t6 = mrSges(6,1) * t135 - mrSges(6,3) * t12;
t5 = t151 * t66 - t251;
t4 = -qJD(5) * t18 - t146 * t30 + t149 * t31;
t3 = qJD(5) * t17 + t146 * t31 + t149 * t30;
t14 = [(-mrSges(4,1) * t111 - mrSges(5,1) * t16 + mrSges(5,2) * t21 + mrSges(4,3) * t24 + (-Ifges(5,3) - Ifges(4,2)) * t66 + 0.2e1 * t258 * t227 - 0.2e1 * t241 * t199) * t157 + t139 * (Ifges(6,5) * t26 + Ifges(6,6) * t27) / 0.2e1 + (-m(4) * t69 - t256) * t43 + 0.2e1 * t193 * t124 * mrSges(3,3) - t26 * t212 + (-mrSges(6,1) * t5 + mrSges(6,3) * t1 + Ifges(6,4) * t12 + Ifges(6,2) * t13 + Ifges(6,6) * t135) * t67 + (t252 / 0.2e1 + Ifges(4,4) * t226 + Ifges(5,5) * t225 + t207 * t242 + t222 * t254 + t243) * t101 + (m(4) * t70 + m(5) * t58 + t240) * t42 + t35 * t60 + t3 * t39 + t4 * t40 + t29 * (-mrSges(6,1) * t27 + mrSges(6,2) * t26) + t26 * t20 / 0.2e1 + t27 * t19 / 0.2e1 + t28 * t23 + t17 * t6 + t18 * t7 - t55 * (Ifges(6,4) * t26 + Ifges(6,2) * t27) / 0.2e1 + t36 * t173 + t132 * t167 - pkin(1) * t168 + (m(4) * t24 + m(5) * t21 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t66 + t48) * t72 + t64 * t179 + (t74 * mrSges(6,1) + t73 * mrSges(6,2) + ((-m(5) - m(4)) * t206 + t244) * t150 + (m(4) * t131 - m(5) * (t176 - t128) - m(6) * (t134 * t151 + t176) - t245) * t148) * g(1) + (-m(6) * (pkin(4) * t196 + t121 + t235) - t76 * mrSges(6,1) + t75 * mrSges(6,2) - m(5) * (t175 + t235) - m(4) * t175 + t245 * t150 + t244 * t148) * g(2) + (mrSges(6,2) * t5 - mrSges(6,3) * t2 + Ifges(6,1) * t12 + Ifges(6,4) * t13 + Ifges(6,5) * t135) * t68 + m(6) * (t1 * t18 + t10 * t4 + t11 * t3 + t17 * t2 + t28 * t29 + t36 * t5) + m(5) * (t16 * t64 + t35 * t41) + (Ifges(6,1) * t26 + Ifges(6,4) * t27) * t228 + (mrSges(4,2) * t111 + mrSges(5,2) * t22 - mrSges(4,3) * t25 - mrSges(5,3) * t16 + (Ifges(5,5) / 0.2e1 - Ifges(4,4) / 0.2e1) * t66 + t254 * t227 + t242 * t199) * t108 + (-m(4) * t25 + m(5) * t22 - qJDD(3) * mrSges(4,1) + mrSges(4,3) * t65 + t47) * t71 + (-m(4) * t111 - t180) * t131 + (t112 * mrSges(4,1) + t41 * mrSges(5,1) + t49 / 0.2e1 - t50 / 0.2e1 - t58 * mrSges(5,2) - t70 * mrSges(4,3) + Ifges(5,3) * t225 - Ifges(4,2) * t226 - t258 * t222 + t241 * t207) * t102 + (qJDD(3) * t242 + t254 * t65 - t258 * t66) * t108 / 0.2e1 + m(3) * (-pkin(1) * t132 + (t124 + t190) * qJ(2) * t193) + (Ifges(3,4) * t143 + Ifges(3,2) * t144) * t187 + (Ifges(3,1) * t143 + Ifges(3,4) * t144) * t188 + t27 * t211 + Ifges(2,3) * qJDD(1); -t240 * t99 + t239 * t100 + m(3) * t132 - t161 * t40 - t55 * t39 - t173 + t168 + t179 + t180 + (-g(1) * t148 + g(2) * t150) * (m(3) + m(4) + t224) - t181 * t193 * qJD(1) ^ 2 + (-t10 * t161 - t11 * t55 - t5) * m(6) + (-t100 * t57 - t58 * t99 + t16) * m(5) + (t100 * t69 - t70 * t99 + t111) * m(4); (t204 + t256) * t70 + t114 * t7 + t113 * t6 + (Ifges(5,2) + Ifges(4,3)) * qJDD(3) + t231 * mrSges(5,1) + qJD(4) * t80 - t59 * t60 + qJ(4) * t48 - pkin(3) * t47 - t34 * t23 - t24 * mrSges(4,2) + t21 * mrSges(5,3) - t243 * t99 - t153 + (-t100 * t112 + t25) * mrSges(4,1) + (-m(6) * (pkin(4) * t134 + t194) + t170 - m(5) * t194 + t249) * g(3) + (-t211 - t19 / 0.2e1 + t29 * mrSges(6,1) + Ifges(6,4) * t229 + Ifges(6,2) * t230 + Ifges(6,6) * t219) * t161 + (-Ifges(4,2) * t100 + t252 + t89) * t225 + (-pkin(3) * t22 - g(1) * t120 - g(2) * t118 + qJ(4) * t21 + t232 * t58 - t41 * t59) * m(5) + t233 * (mrSges(4,1) * t133 + mrSges(4,2) * t134) + (-t29 * mrSges(6,2) - t20 / 0.2e1 + Ifges(6,1) * t229 + Ifges(6,4) * t230 + t212 + Ifges(6,5) * t219) * t55 - (t254 * t99 - t203 + t49 + t88) * t100 / 0.2e1 - (t241 * t100 + t242 * t99) * qJD(3) / 0.2e1 + (Ifges(5,3) * t100 + t215) * t226 + t58 * t205 + t50 * t222 + t237 * t40 + t238 * t39 + (t1 * t114 + t113 * t2 - t29 * t34 - g(1) * (t150 * t185 + t120) - g(2) * (t148 * t185 + t118) + t238 * t11 + t237 * t10) * m(6) - t240 * t69 + t241 * t66 + t242 * t65 + (-t150 * t155 + t171) * g(1) + (-t148 * t155 + t172) * g(2); t146 * t7 + t149 * t6 + (-t23 + t60) * t100 + t162 * qJD(5) + t224 * t134 * g(3) + (-t162 - t80) * qJD(3) + t47 + (t1 * t146 - t100 * t29 - t236 + t149 * t2 + t139 * (-t10 * t146 + t11 * t149)) * m(6) + (-qJD(3) * t58 - t231 - t236) * m(5); -t29 * (mrSges(6,1) * t161 - mrSges(6,2) * t55) + (-Ifges(6,1) * t55 - t216) * t229 + t19 * t228 + (-Ifges(6,5) * t55 - Ifges(6,6) * t161) * t219 - t10 * t39 + t11 * t40 - g(1) * t171 - g(2) * t172 - g(3) * t170 + (-t10 * t55 + t11 * t161) * mrSges(6,3) + t153 + (-Ifges(6,2) * t161 + t20 - t46) * t230;];
tau = t14;
