% Calculate vector of inverse dynamics joint torques for
% S5RPRPP3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:20
% EndTime: 2019-12-31 18:12:34
% DurationCPUTime: 7.32s
% Computational Cost: add. (1783->355), mult. (4305->396), div. (0->0), fcn. (2794->8), ass. (0->157)
t230 = -mrSges(5,1) - mrSges(4,3);
t229 = -mrSges(5,2) + mrSges(4,1);
t224 = Ifges(6,4) + Ifges(5,5);
t223 = Ifges(4,5) + Ifges(6,5);
t222 = Ifges(6,2) + Ifges(5,3);
t220 = Ifges(6,3) + Ifges(4,1);
t155 = qJD(1) * qJD(2);
t92 = qJDD(1) * qJ(2) + t155;
t228 = mrSges(4,2) - mrSges(6,2) - mrSges(5,3);
t109 = sin(pkin(7));
t110 = cos(pkin(7));
t157 = t109 ^ 2 + t110 ^ 2;
t184 = cos(qJ(3));
t148 = t184 * t110;
t136 = qJD(1) * t148;
t141 = qJDD(1) * t184;
t113 = sin(qJ(3));
t156 = qJD(3) * t113;
t145 = t109 * t156;
t151 = qJDD(1) * t113;
t45 = qJD(1) * t145 - qJD(3) * t136 - t109 * t141 - t110 * t151;
t227 = t45 / 0.2e1;
t76 = t109 * t184 + t113 * t110;
t72 = t76 * qJD(3);
t46 = qJD(1) * t72 + t109 * t151 - t110 * t141;
t226 = -t46 / 0.2e1;
t114 = sin(qJ(1));
t180 = g(2) * t114;
t225 = -Ifges(4,4) + Ifges(6,6);
t221 = Ifges(5,6) - Ifges(6,6);
t70 = t76 * qJD(1);
t172 = t70 * Ifges(5,6);
t62 = Ifges(6,6) * t70;
t158 = t109 * t113;
t69 = qJD(1) * t158 - t136;
t219 = qJD(3) * t224 + t222 * t69 - t172 + t62;
t182 = Ifges(6,6) * t69;
t64 = Ifges(4,4) * t69;
t218 = qJD(3) * t223 + t220 * t70 + t182 - t64;
t51 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t69;
t54 = mrSges(5,1) * t69 - qJD(3) * mrSges(5,3);
t217 = t51 - t54;
t206 = t229 * qJD(3) + t230 * t70;
t216 = -qJDD(3) / 0.2e1;
t115 = cos(qJ(1));
t202 = g(1) * t115 + t180;
t215 = -t69 * pkin(4) + qJD(5);
t167 = pkin(6) + qJ(2);
t82 = t167 * t110;
t78 = qJD(1) * t82;
t161 = t113 * t78;
t81 = t167 * t109;
t77 = qJD(1) * t81;
t47 = t184 * t77 + t161;
t124 = pkin(4) * t70 + t47;
t214 = t124 + qJD(4);
t198 = -Ifges(4,6) + t224;
t213 = -Ifges(5,4) + t223;
t106 = pkin(7) + qJ(3);
t102 = sin(t106);
t103 = cos(t106);
t212 = t102 * t228 - t103 * t229;
t144 = m(3) * qJ(2) + mrSges(3,3);
t210 = -t144 + mrSges(2,2) - m(6) * (pkin(4) + t167) - mrSges(6,1) + t230;
t131 = -mrSges(3,1) * t110 + mrSges(3,2) * t109;
t209 = -m(3) * pkin(1) - mrSges(2,1) + t131 + t212;
t100 = pkin(2) * t110 + pkin(1);
t80 = -qJD(1) * t100 + qJD(2);
t119 = -qJ(4) * t70 + t80;
t168 = pkin(3) + qJ(5);
t10 = t168 * t69 + t119;
t48 = -t113 * t77 + t184 * t78;
t36 = -qJD(3) * qJ(4) - t48;
t14 = -t36 + t215;
t20 = pkin(3) * t69 + t119;
t208 = t80 * mrSges(4,1) + t36 * mrSges(5,1) - t14 * mrSges(6,1) - t20 * mrSges(5,2) - t48 * mrSges(4,3) + t10 * mrSges(6,3);
t13 = -qJD(3) * t168 + t214;
t200 = -t47 - qJD(4);
t35 = -qJD(3) * pkin(3) - t200;
t207 = -t35 * mrSges(5,1) - t13 * mrSges(6,1) - t80 * mrSges(4,2) + t10 * mrSges(6,2) - t47 * mrSges(4,3) + t20 * mrSges(5,3);
t55 = -mrSges(6,1) * t69 + qJD(3) * mrSges(6,2);
t166 = t54 - t55;
t159 = t103 * t115;
t96 = t102 * qJ(4);
t205 = pkin(3) * t159 + t115 * t96;
t204 = t102 * t202;
t201 = -t48 - t215;
t142 = qJD(3) * t184;
t21 = t113 * (qJD(2) * t109 + qJD(3) * t82) - qJD(2) * t148 + t81 * t142;
t140 = pkin(6) * qJDD(1) + t92;
t60 = t140 * t109;
t61 = t140 * t110;
t150 = -t113 * t60 - t77 * t142 + t184 * t61;
t5 = -qJDD(3) * qJ(4) + qJD(3) * (-qJD(4) + t161) - t150;
t197 = (-g(1) * t159 - t103 * t180) * qJ(4);
t196 = -t45 / 0.2e1;
t195 = t46 / 0.2e1;
t194 = -t69 / 0.2e1;
t193 = t69 / 0.2e1;
t192 = -t70 / 0.2e1;
t191 = t70 / 0.2e1;
t186 = m(5) + m(6);
t97 = t103 * pkin(3);
t173 = t70 * Ifges(4,4);
t169 = qJD(3) / 0.2e1;
t164 = t97 + t96;
t163 = qJ(4) * t69;
t160 = qJDD(3) / 0.2e1;
t153 = qJDD(1) * t109;
t152 = qJDD(1) * t110;
t143 = t45 * mrSges(6,2) + t46 * mrSges(6,3);
t28 = -t45 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t27 = -t46 * mrSges(6,1) + qJDD(3) * mrSges(6,2);
t49 = t113 * t82 + t184 * t81;
t86 = t115 * t100;
t139 = t114 * t167 + t86;
t25 = -t45 * mrSges(6,1) - qJDD(3) * mrSges(6,3);
t138 = -t100 - t96;
t135 = -m(6) * t168 - mrSges(6,3);
t133 = -mrSges(3,1) * t152 + mrSges(3,2) * t153;
t71 = -t110 * t142 + t145;
t127 = qJ(4) * t71 - qJD(4) * t76;
t125 = -qJ(4) * t76 - t100;
t50 = -t113 * t81 + t184 * t82;
t9 = -t113 * t61 - t78 * t142 + t77 * t156 - t184 * t60;
t79 = -qJDD(1) * t100 + qJDD(2);
t120 = qJDD(4) - t9;
t117 = qJ(4) * t45 - qJD(4) * t70 + t79;
t22 = qJD(2) * t76 + qJD(3) * t50;
t101 = -qJDD(1) * pkin(1) + qJDD(2);
t75 = -t148 + t158;
t63 = Ifges(5,6) * t69;
t53 = mrSges(6,1) * t70 - qJD(3) * mrSges(6,3);
t44 = pkin(3) * t75 + t125;
t41 = t45 * mrSges(4,2);
t40 = t45 * mrSges(5,3);
t39 = -mrSges(5,2) * t69 - mrSges(5,3) * t70;
t38 = pkin(3) * t70 + t163;
t37 = -mrSges(6,2) * t70 + mrSges(6,3) * t69;
t33 = -t69 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t173;
t32 = Ifges(5,4) * qJD(3) - t70 * Ifges(5,2) + t63;
t26 = mrSges(5,1) * t46 - qJDD(3) * mrSges(5,3);
t24 = -t75 * pkin(4) + t50;
t23 = pkin(4) * t76 + t49;
t17 = t168 * t75 + t125;
t16 = pkin(3) * t72 + t127;
t15 = t168 * t70 + t163;
t12 = -t71 * pkin(4) + t22;
t11 = -pkin(4) * t72 - t21;
t8 = -t156 * t78 + t150;
t7 = -qJDD(3) * pkin(3) + t120;
t6 = qJD(5) * t75 + t168 * t72 + t127;
t4 = pkin(3) * t46 + t117;
t3 = -pkin(4) * t46 + qJDD(5) - t5;
t2 = -t45 * pkin(4) - qJD(3) * qJD(5) - qJDD(3) * t168 + t120;
t1 = qJD(5) * t69 + t168 * t46 + t117;
t18 = [-pkin(1) * t133 + (t219 / 0.2e1 - t33 / 0.2e1 + Ifges(5,6) * t192 - Ifges(4,2) * t194 + t222 * t193 + t225 * t191 + t198 * t169 + t208) * t72 + t101 * t131 + m(5) * (t16 * t20 + t4 * t44) + 0.2e1 * t157 * t92 * mrSges(3,3) + t17 * t143 + (m(4) * t8 - m(5) * t5 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t46 - t26) * t50 + t12 * t53 + t11 * t55 + t44 * (-t46 * mrSges(5,2) + t40) + t6 * t37 + t16 * t39 + t23 * t25 + t24 * t27 + (-m(4) * t9 + m(5) * t7 - qJDD(3) * mrSges(4,1) - mrSges(4,3) * t45 + t28) * t49 + (-mrSges(6,2) * t1 + Ifges(4,4) * t226 - Ifges(5,2) * t45 - mrSges(5,3) * t4 + Ifges(5,4) * t216 + mrSges(4,2) * t79 - mrSges(4,3) * t9 + mrSges(6,1) * t2 + mrSges(5,1) * t7 + Ifges(6,6) * t195 + t220 * t196 + (t226 - t195) * Ifges(5,6) + t213 * t160) * t76 + (Ifges(4,2) * t46 + Ifges(4,6) * t216 + t1 * mrSges(6,3) + Ifges(5,6) * t227 - t4 * mrSges(5,2) + t79 * mrSges(4,1) - t8 * mrSges(4,3) - t3 * mrSges(6,1) + t5 * mrSges(5,1) + Ifges(6,6) * t196 + t222 * t195 + (t227 - t196) * Ifges(4,4) + t198 * t160) * t75 + (-m(4) * t79 - t46 * mrSges(4,1) + t41) * t100 + m(3) * (-pkin(1) * t101 + (t155 + t92) * qJ(2) * t157) + (-t218 / 0.2e1 + t32 / 0.2e1 + Ifges(5,2) * t192 - Ifges(4,4) * t194 + t221 * t193 - t220 * t191 - t213 * t169 + t207) * t71 + (((-m(5) - m(4)) * t167 + t210) * t115 + (m(4) * t100 - m(5) * (t138 - t97) - m(6) * t138 - t103 * t135 - t209) * t114) * g(1) + (-m(4) * t139 - m(5) * (t139 + t205) - m(6) * (t86 + t205) + (-(m(6) * qJ(5) + mrSges(6,3)) * t103 + t209) * t115 + t210 * t114) * g(2) + (-m(4) * t48 + m(5) * t36 - t217) * t21 + (qJDD(3) * t224 + t221 * t45 + t222 * t46) * t75 / 0.2e1 + (qJDD(3) * t223 - t220 * t45 + t225 * t46) * t76 / 0.2e1 + (m(4) * t47 + m(5) * t35 - t206) * t22 + (Ifges(3,4) * t109 + Ifges(3,2) * t110) * t152 + (Ifges(3,1) * t109 + Ifges(3,4) * t110) * t153 + Ifges(2,3) * qJDD(1) + m(6) * (t1 * t17 + t10 * t6 + t11 * t14 + t12 * t13 + t2 * t23 + t24 * t3); t40 - t41 + t229 * t46 + (-t53 + t206) * t70 + (t51 - t166) * t69 + m(3) * t101 + t133 + t143 + (-g(1) * t114 + g(2) * t115) * (m(3) + m(4) + t186) - t144 * t157 * qJD(1) ^ 2 + (-t13 * t70 + t14 * t69 + t1) * m(6) + (-t35 * t70 - t36 * t69 + t4) * m(5) + (-t47 * t70 + t48 * t69 + t79) * m(4); (Ifges(4,3) + Ifges(6,1) + Ifges(5,1)) * qJDD(3) + (t228 * t103 + (m(5) * pkin(3) - t135 + t229) * t102) * t202 + (-t26 + t27) * qJ(4) + (t172 + t33) * t191 - t166 * qJD(4) + t206 * t48 + t201 * t53 + t198 * t46 + (-pkin(3) * t7 - qJ(4) * t5 - t20 * t38 + t200 * t36 - t35 * t48 + t197) * m(5) + t124 * t55 - t168 * t25 - t213 * t45 - (t198 * t70 - t213 * t69) * qJD(3) / 0.2e1 - pkin(3) * t28 - t15 * t37 - t38 * t39 - t5 * mrSges(5,3) + t7 * mrSges(5,2) - t8 * mrSges(4,2) + t9 * mrSges(4,1) - t2 * mrSges(6,3) + t3 * mrSges(6,2) + (Ifges(5,2) * t191 - t207) * t69 - t208 * t70 + t217 * t47 + (-Ifges(4,2) * t70 + t218 - t64) * t193 + (-t220 * t69 - t173 + t219 + t62) * t192 + (t222 * t70 - t182 + t32 + t63) * t194 + (-m(6) * (qJ(5) * t103 + t164) - t103 * mrSges(6,3) - m(5) * t164 + t212) * g(3) + (qJ(4) * t3 - t10 * t15 + t201 * t13 + t14 * t214 - t168 * t2 + t197) * m(6); (t37 + t39) * t70 + t166 * qJD(3) + t186 * t103 * g(3) + t25 + t28 + (-qJD(3) * t14 + t10 * t70 + t2 - t204) * m(6) + (qJD(3) * t36 + t20 * t70 - t204 + t7) * m(5); qJD(3) * t53 - t69 * t37 + (-g(3) * t102 + t13 * qJD(3) - t10 * t69 - t103 * t202 + t3) * m(6) + t27;];
tau = t18;
