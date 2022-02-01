% Calculate vector of inverse dynamics joint torques for
% S5RPPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 08:59:04
% EndTime: 2022-01-23 08:59:19
% DurationCPUTime: 9.15s
% Computational Cost: add. (2714->380), mult. (6797->544), div. (0->0), fcn. (5021->10), ass. (0->173)
t150 = sin(pkin(7));
t153 = cos(pkin(7));
t197 = t150 ^ 2 + t153 ^ 2;
t190 = qJD(1) * qJD(2);
t136 = qJDD(1) * qJ(2) + t190;
t148 = sin(pkin(9));
t151 = cos(pkin(9));
t124 = pkin(4) * t151 + pkin(6) * t148 + pkin(3);
t149 = sin(pkin(8));
t152 = cos(pkin(8));
t214 = qJ(4) * t152;
t168 = pkin(3) * t149 - t214;
t164 = qJ(2) + t168;
t241 = qJ(2) * (-m(3) - m(4)) - m(5) * t164 + m(6) * (-t124 * t149 - qJ(2) + t214) + mrSges(2,2) - mrSges(3,3);
t195 = qJD(1) * t150;
t181 = t149 * t195;
t194 = qJD(1) * t153;
t173 = t150 * qJ(3) + pkin(1);
t123 = pkin(2) * t153 + t173;
t101 = -qJD(1) * t123 + qJD(2);
t182 = qJ(2) * t194;
t64 = t149 * t101 + t152 * t182;
t50 = -qJ(4) * t194 + t64;
t130 = qJ(2) * t195 + qJD(3);
t75 = t168 * t195 + t130;
t25 = -t148 * t50 + t151 * t75;
t20 = -pkin(4) * t181 - t25;
t240 = m(6) * t20;
t154 = sin(qJ(5));
t156 = cos(qJ(5));
t187 = qJDD(1) * t150;
t175 = t149 * t187;
t196 = qJD(1) * t149;
t179 = t156 * t196;
t176 = t151 * t195;
t177 = t148 * t194;
t94 = t152 * t176 - t177;
t52 = t150 * t179 - t154 * t94;
t205 = t152 * t150;
t103 = -t148 * t153 + t151 * t205;
t91 = t103 * qJDD(1);
t27 = qJD(5) * t52 + t154 * t175 + t156 * t91;
t239 = Ifges(6,5) * t27;
t180 = t154 * t196;
t54 = t150 * t180 + t156 * t94;
t28 = -qJD(5) * t54 - t154 * t91 + t156 * t175;
t238 = Ifges(6,6) * t28;
t102 = t148 * t205 + t151 * t153;
t90 = t102 * qJDD(1);
t85 = qJDD(5) + t90;
t237 = Ifges(6,3) * t85;
t236 = mrSges(5,3) - mrSges(4,2);
t235 = -mrSges(6,3) + mrSges(5,2);
t218 = -mrSges(5,1) * t181 - mrSges(6,1) * t52 + mrSges(6,2) * t54 + mrSges(5,3) * t94;
t166 = -mrSges(4,1) * t153 - mrSges(4,3) * t205;
t92 = t102 * qJD(1);
t215 = -mrSges(5,1) * t92 - mrSges(5,2) * t94 + t166 * qJD(1);
t203 = t153 * t152;
t104 = t148 * t150 + t151 * t203;
t234 = t104 * t156;
t210 = t149 * t156;
t189 = qJD(1) * qJD(4);
t192 = qJD(3) * t150;
t79 = -qJD(1) * t192 - qJDD(1) * t123 + qJDD(2);
t46 = t136 * t203 + t149 * t79;
t37 = (-qJ(4) * qJDD(1) - t189) * t153 + t46;
t120 = t150 * t136 + qJDD(3);
t47 = (qJDD(1) * t168 - t152 * t189) * t150 + t120;
t13 = t148 * t47 + t151 * t37;
t11 = pkin(6) * t175 + t13;
t186 = qJDD(1) * t153;
t204 = t153 * t149;
t45 = -t136 * t204 + t152 * t79;
t38 = pkin(3) * t186 + qJDD(4) - t45;
t18 = pkin(4) * t90 - pkin(6) * t91 + t38;
t63 = t101 * t152 - t149 * t182;
t49 = pkin(3) * t194 + qJD(4) - t63;
t19 = pkin(4) * t92 - pkin(6) * t94 + t49;
t26 = t148 * t75 + t151 * t50;
t21 = pkin(6) * t181 + t26;
t5 = -t154 * t21 + t156 * t19;
t1 = qJD(5) * t5 + t11 * t156 + t154 * t18;
t6 = t154 * t19 + t156 * t21;
t2 = -qJD(5) * t6 - t11 * t154 + t156 * t18;
t231 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t169 = -mrSges(3,1) * t153 + mrSges(3,2) * t150;
t219 = t149 * qJ(4) + pkin(2);
t229 = -m(4) * t123 - m(5) * ((pkin(3) * t152 + t219) * t153 + t173) - m(6) * ((t124 * t152 + t219) * t153 + pkin(1) + (pkin(4) * t148 - pkin(6) * t151 + qJ(3)) * t150) - mrSges(2,1) - m(3) * pkin(1) + t169;
t211 = t149 * t154;
t109 = t151 * t211 + t152 * t156;
t95 = t104 * qJD(1);
t228 = -t109 * qJD(5) - t153 * t180 - t156 * t95;
t110 = t151 * t210 - t152 * t154;
t227 = -t110 * qJD(5) - t153 * t179 + t154 * t95;
t226 = qJD(1) ^ 2 * t197;
t225 = qJD(5) + t92;
t155 = sin(qJ(1));
t157 = cos(qJ(1));
t67 = -t104 * t154 + t156 * t204;
t224 = -t155 * t109 + t157 * t67;
t222 = t54 / 0.2e1;
t59 = mrSges(5,1) * t175 - mrSges(5,3) * t91;
t7 = -mrSges(6,1) * t28 + mrSges(6,2) * t27;
t221 = t59 - t7;
t220 = Ifges(6,4) * t54;
t82 = qJ(2) * t203 - t149 * t123;
t72 = -qJ(4) * t153 + t82;
t88 = t164 * t150;
t36 = t148 * t88 + t151 * t72;
t213 = t148 * t149;
t212 = t149 * t150;
t209 = t150 * t155;
t208 = t150 * t157;
t207 = t151 * t155;
t206 = t151 * t157;
t202 = t153 * t155;
t201 = t153 * t157;
t174 = t152 * t187;
t198 = mrSges(4,1) * t175 + mrSges(4,2) * t174;
t193 = qJD(2) * t153;
t191 = m(4) + m(5) + m(6);
t185 = t237 + t238 + t239;
t40 = t90 * mrSges(5,1) + t91 * mrSges(5,2);
t81 = -qJ(2) * t204 - t123 * t152;
t65 = -t103 * t154 + t150 * t210;
t66 = t103 * t156 + t150 * t211;
t172 = -t65 * mrSges(6,1) + mrSges(6,2) * t66;
t74 = t153 * pkin(3) - t81;
t170 = -mrSges(3,1) * t186 + mrSges(3,2) * t187;
t106 = -t149 * t192 + t152 * t193;
t12 = -t148 * t37 + t151 * t47;
t35 = -t148 * t72 + t151 * t88;
t29 = pkin(4) * t102 - pkin(6) * t103 + t74;
t31 = pkin(6) * t212 + t36;
t8 = -t154 * t31 + t156 * t29;
t9 = t154 * t29 + t156 * t31;
t167 = (-qJD(4) * t152 + qJD(2)) * t150;
t165 = mrSges(4,2) * t153 - mrSges(4,3) * t212;
t163 = t154 * t204 + t234;
t161 = (mrSges(4,1) * t149 + mrSges(4,2) * t152) * t150;
t111 = t149 * t202 + t152 * t157;
t113 = t149 * t201 - t152 * t155;
t159 = -g(1) * t113 - g(2) * t111 - g(3) * t212;
t141 = -qJDD(1) * pkin(1) + qJDD(2);
t115 = t165 * qJD(1);
t114 = t149 * t155 + t152 * t201;
t112 = t149 * t157 - t152 * t202;
t108 = t166 * qJDD(1);
t107 = t165 * qJDD(1);
t105 = t149 * t193 + t152 * t192;
t99 = qJD(1) * t161;
t96 = t157 * t110;
t93 = t152 * t177 - t176;
t83 = -qJD(4) * t153 + t106;
t77 = t110 * t195;
t76 = t109 * t195;
t61 = -mrSges(5,2) * t181 - mrSges(5,3) * t92;
t58 = -mrSges(5,2) * t175 - mrSges(5,3) * t90;
t57 = t66 * qJD(5);
t56 = t65 * qJD(5);
t51 = Ifges(6,4) * t52;
t43 = t148 * t167 + t151 * t83;
t33 = mrSges(6,1) * t225 - mrSges(6,3) * t54;
t32 = -mrSges(6,2) * t225 + mrSges(6,3) * t52;
t30 = -pkin(4) * t212 - t35;
t17 = t54 * Ifges(6,1) + Ifges(6,5) * t225 + t51;
t16 = Ifges(6,2) * t52 + Ifges(6,6) * t225 + t220;
t15 = -mrSges(6,2) * t85 + mrSges(6,3) * t28;
t14 = mrSges(6,1) * t85 - mrSges(6,3) * t27;
t10 = -pkin(4) * t175 - t12;
t4 = -qJD(5) * t9 + t105 * t156 - t154 * t43;
t3 = qJD(5) * t8 + t105 * t154 + t156 * t43;
t22 = [t225 * (Ifges(6,5) * t56 - Ifges(6,6) * t57) / 0.2e1 + (mrSges(5,2) * t38 - mrSges(5,3) * t12 + Ifges(5,1) * t91 - Ifges(5,4) * t90) * t103 + (-t5 * t56 - t57 * t6) * mrSges(6,3) + (-(t112 * t151 - t148 * t209) * mrSges(5,1) - t112 * mrSges(4,1) + mrSges(4,3) * t209 - t96 * mrSges(6,1) + t235 * (t112 * t148 + t150 * t207) + (t163 * mrSges(6,1) + t67 * mrSges(6,2) - t229) * t155 + t236 * t111 + (t109 * mrSges(6,2) + t241) * t157) * g(1) + (-(t114 * t151 + t148 * t208) * mrSges(5,1) - t207 * mrSges(6,1) * t210 - t224 * mrSges(6,2) - t114 * mrSges(4,1) - mrSges(4,3) * t208 + t235 * (t114 * t148 - t150 * t206) + (-mrSges(6,1) * t234 + t229) * t157 + (-t154 * mrSges(6,1) - t236) * t113 + t241 * t155) * g(2) + (-mrSges(6,3) * t2 + Ifges(6,1) * t27 + Ifges(6,4) * t28 + Ifges(6,5) * t85) * t66 + (Ifges(6,1) * t56 - Ifges(6,4) * t57) * t222 + t20 * (mrSges(6,1) * t57 + mrSges(6,2) * t56) + t52 * (Ifges(6,4) * t56 - Ifges(6,2) * t57) / 0.2e1 + (mrSges(5,1) * t12 - mrSges(5,2) * t13 + Ifges(5,5) * t91 - Ifges(5,6) * t90) * t212 + m(4) * (t106 * t64 + t45 * t81 + t46 * t82) + 0.2e1 * t197 * t136 * mrSges(3,3) + m(5) * (t12 * t35 + t13 * t36 + t26 * t43 + t38 * t74) + (-m(5) * t25 + t218 + t240) * (t148 * t83 - t151 * t167) + (-Ifges(5,4) * t91 + Ifges(5,2) * t90 + t38 * mrSges(5,1) + t237 / 0.2e1 + t239 / 0.2e1 + t238 / 0.2e1 - t13 * mrSges(5,3) + t185 / 0.2e1 + t231) * t102 + (-m(4) * t63 + m(5) * t49 - t215) * t105 + t120 * t161 + m(3) * (-pkin(1) * t141 + (t136 + t190) * qJ(2) * t197) + t106 * t115 + t82 * t107 + t81 * t108 + t74 * t40 - t57 * t16 / 0.2e1 + t36 * t58 + t35 * t59 + t43 * t61 + t56 * t17 / 0.2e1 + t30 * t7 + t3 * t32 + t4 * t33 + t8 * t14 + t9 * t15 + (Ifges(3,1) * t187 + qJD(2) * t99 + qJ(2) * t198 + (Ifges(4,1) * t152 - Ifges(4,4) * t149) * t174 + m(4) * (qJ(2) * t120 + qJD(2) * t130) + (-Ifges(4,5) * t152 + Ifges(4,6) * t149 + Ifges(3,4)) * t186) * t150 + m(6) * (t1 * t9 + t10 * t30 + t2 * t8 + t3 * t6 + t4 * t5) + (Ifges(3,4) * t187 - Ifges(4,5) * t174 + (Ifges(3,2) + Ifges(4,3)) * t186) * t153 + (mrSges(6,3) * t1 + Ifges(6,4) * t27 + Ifges(6,2) * t28 + Ifges(6,6) * t85) * t65 + t10 * t172 + t141 * t169 - pkin(1) * t170 + (Ifges(5,3) * t212 - Ifges(5,6) * t102 + Ifges(5,5) * t103 + Ifges(4,6) * t153 - (Ifges(4,4) * t152 - Ifges(4,2) * t149) * t150) * t175 + t46 * t165 + t45 * t166 + Ifges(2,3) * qJDD(1); (-t148 * t221 + t151 * t58 + t194 * t215 + t107) * t149 - mrSges(3,3) * t226 - t109 * t14 + t110 * t15 - t95 * t61 - t99 * t195 - t218 * t93 + (-t115 * t194 + t108 - t40) * t152 + t170 + t227 * t33 + t228 * t32 + (-g(1) * t155 + g(2) * t157) * (m(3) + t191) + (t1 * t110 + t10 * t213 - t109 * t2 - t20 * t93 + t227 * t5 + t228 * t6) * m(6) + (-(t130 * t150 + t203 * t64 - t204 * t63) * qJD(1) + t149 * t46 + t152 * t45) * m(4) + (-qJ(2) * t226 + t141) * m(3) + (t25 * t93 - t26 * t95 - t49 * t204 * qJD(1) - t152 * t38 + (-t12 * t148 + t13 * t151) * t149) * m(5); t77 * t32 - t76 * t33 + t221 * t151 + t191 * t153 * g(3) + (-t154 * t14 + t156 * t15 + t58 + (-t154 * t32 - t156 * t33) * qJD(5)) * t148 + m(5) * (t12 * t151 + t13 * t148) + m(4) * t120 + ((t215 * t152 + (t148 * t218 + t151 * t61 + t115) * t149 - m(5) * (-t149 * t151 * t26 + t152 * t49 + t213 * t25) - m(4) * (-t149 * t64 - t152 * t63) + t213 * t240) * qJD(1) + (-g(1) * t157 - g(2) * t155) * t191) * t150 + t198 + (-t10 * t151 + (t1 * t156 - t154 * t2 + (-t154 * t6 - t156 * t5) * qJD(5)) * t148 - t5 * t76 + t6 * t77) * m(6); t92 * t61 - t218 * t94 + (t225 * t32 + t14) * t156 + (-t225 * t33 + t15) * t154 + t40 + (t1 * t154 + t156 * t2 - t20 * t94 + t159 + t225 * (-t154 * t5 + t156 * t6)) * m(6) + (t25 * t94 + t26 * t92 + t159 + t38) * m(5); -t20 * (mrSges(6,1) * t54 + mrSges(6,2) * t52) - t54 * (Ifges(6,1) * t52 - t220) / 0.2e1 + t16 * t222 - t225 * (Ifges(6,5) * t52 - Ifges(6,6) * t54) / 0.2e1 - t5 * t32 + t6 * t33 - g(1) * (t224 * mrSges(6,1) + (-t110 * t155 - t157 * t163) * mrSges(6,2)) - g(2) * ((-(t104 * t155 - t149 * t206) * t154 + t156 * t111) * mrSges(6,1) + (-t155 * t163 + t96) * mrSges(6,2)) + g(3) * t172 + (t5 * t52 + t54 * t6) * mrSges(6,3) + t185 - (-Ifges(6,2) * t54 + t17 + t51) * t52 / 0.2e1 + t231;];
tau = t22;
