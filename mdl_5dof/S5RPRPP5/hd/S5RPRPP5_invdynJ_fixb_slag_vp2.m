% Calculate vector of inverse dynamics joint torques for
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:08
% EndTime: 2019-12-31 18:16:16
% DurationCPUTime: 4.82s
% Computational Cost: add. (999->313), mult. (1879->363), div. (0->0), fcn. (747->4), ass. (0->143)
t212 = Ifges(6,4) + Ifges(5,5);
t214 = -Ifges(5,1) - Ifges(6,1);
t222 = Ifges(4,1) - t214;
t213 = Ifges(5,4) + Ifges(4,5);
t195 = -Ifges(6,5) + t213;
t78 = sin(qJ(3));
t221 = t212 * t78;
t220 = Ifges(6,2) + Ifges(5,3);
t219 = -Ifges(4,6) + Ifges(5,6);
t210 = Ifges(5,6) - Ifges(6,6);
t172 = Ifges(4,4) * t78;
t80 = cos(qJ(3));
t218 = t222 * t80 - t172 + t221;
t149 = qJD(1) * t80;
t217 = t212 * t149;
t181 = pkin(3) * t78;
t70 = t80 * qJ(4);
t216 = t70 - t181;
t153 = qJ(4) * t78;
t98 = pkin(3) * t80 + t153;
t144 = qJ(5) * qJD(1);
t83 = -pkin(1) - pkin(6);
t58 = qJD(1) * t83 + qJD(2);
t17 = (t58 + t144) * t80;
t215 = -t17 + qJD(4);
t183 = -m(3) - m(4);
t182 = -m(5) - m(6);
t46 = t78 * t58;
t33 = qJD(3) * qJ(4) + t46;
t64 = t78 * t144;
t12 = t64 + t33;
t140 = qJD(1) * qJD(5);
t148 = qJD(3) * t78;
t37 = t58 * t148;
t141 = qJD(1) * qJD(3);
t44 = qJDD(1) * t80 - t141 * t78;
t57 = qJDD(1) * t83 + qJDD(2);
t82 = -pkin(3) - pkin(4);
t2 = -qJ(5) * t44 + qJDD(4) + t37 + (-t57 - t140) * t80 + t82 * qJDD(3);
t209 = t12 * qJD(3) - t2;
t36 = t44 * mrSges(5,2);
t21 = -qJDD(3) * mrSges(5,1) + t36;
t208 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t44 - t21;
t45 = qJDD(1) * t78 + t141 * t80;
t24 = -mrSges(5,2) * t45 + qJDD(3) * mrSges(5,3);
t207 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t45 + t24;
t150 = qJD(1) * t78;
t206 = t210 * qJD(3) + t220 * t150 + t217;
t132 = mrSges(6,3) * t150;
t50 = qJD(3) * mrSges(6,2) + t132;
t133 = mrSges(5,2) * t150;
t55 = qJD(3) * mrSges(5,3) - t133;
t205 = t50 + t55;
t131 = mrSges(5,2) * t149;
t204 = -mrSges(4,3) * t149 - t131 + (mrSges(4,1) + mrSges(5,1)) * qJD(3);
t51 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t150;
t203 = t55 + t51;
t108 = t80 * mrSges(5,1) + t78 * mrSges(5,3);
t34 = t78 * t82 - qJ(2) + t70;
t11 = qJD(1) * t34 + qJD(5);
t162 = t78 * mrSges(6,2);
t49 = qJ(2) - t216;
t26 = t49 * qJD(1);
t202 = t26 * t108 + t11 * (-t80 * mrSges(6,1) - t162);
t81 = cos(qJ(1));
t178 = g(2) * t81;
t79 = sin(qJ(1));
t179 = g(1) * t79;
t196 = t178 - t179;
t201 = t196 * t80;
t72 = t80 * mrSges(6,2);
t106 = -t78 * mrSges(6,1) + t72;
t200 = -t213 * t78 + t219 * t80;
t139 = qJDD(1) * qJ(2);
t142 = qJD(1) * qJD(2);
t59 = t139 + t142;
t199 = t212 * t80;
t8 = t57 * t80 - t37;
t147 = qJD(3) * t80;
t9 = t58 * t147 + t78 * t57;
t198 = -t78 * t9 - t8 * t80;
t5 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t9;
t6 = -qJDD(3) * pkin(3) + qJDD(4) - t8;
t197 = t5 * t78 - t6 * t80;
t194 = t218 * qJD(1) + t195 * qJD(3);
t171 = Ifges(4,4) * t80;
t193 = -t80 * t171 + ((-Ifges(4,1) + t220) * t80 - t221) * t78;
t116 = qJD(4) * t80 - qJD(2);
t192 = -qJ(4) * t44 - qJD(1) * t116 + t139;
t109 = t78 * mrSges(4,1) + t80 * mrSges(4,2);
t164 = t78 * mrSges(5,1);
t191 = -mrSges(3,3) - t109 - t164 - (-m(5) * qJ(4) - mrSges(5,3)) * t80 - m(6) * (pkin(4) * t78 - t70) + t106 + mrSges(2,2);
t190 = m(6) * qJ(5) - mrSges(2,1) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t189 = qJD(1) ^ 2;
t3 = qJ(5) * t45 + t140 * t78 + t5;
t177 = t3 * t78;
t166 = t44 * mrSges(6,3);
t165 = t58 * t80;
t107 = -t80 * mrSges(5,3) + t164;
t40 = t107 * qJD(1);
t41 = t106 * qJD(1);
t158 = t40 - t41;
t156 = t98 * t79;
t154 = t81 * pkin(1) + t79 * qJ(2);
t151 = qJ(5) + t83;
t143 = qJDD(1) * mrSges(3,2);
t136 = t81 * pkin(6) + t154;
t130 = mrSges(6,3) * t149;
t124 = -t45 * mrSges(6,1) + t44 * mrSges(6,2);
t122 = (t142 + t59) * qJ(2);
t121 = -t141 / 0.2e1;
t120 = t141 / 0.2e1;
t118 = qJD(3) * t151;
t115 = m(6) * t82 - mrSges(6,1);
t110 = mrSges(4,1) * t80 - mrSges(4,2) * t78;
t102 = -t78 * Ifges(4,2) + t171;
t99 = -Ifges(6,5) * t78 + Ifges(6,6) * t80;
t25 = -qJD(3) * pkin(3) + qJD(4) - t165;
t97 = t25 * t78 + t33 * t80;
t95 = t80 * t82 - t153;
t92 = t78 * (-Ifges(4,2) * t80 - t172);
t88 = qJ(2) * t110;
t85 = qJD(3) * t97 + t197;
t71 = t81 * qJ(2);
t67 = -qJDD(1) * pkin(1) + qJDD(2);
t52 = -qJD(3) * mrSges(6,1) - t130;
t48 = t151 * t80;
t47 = t151 * t78;
t42 = t109 * qJD(1);
t39 = t98 * qJD(1);
t29 = Ifges(4,6) * qJD(3) + qJD(1) * t102;
t22 = qJDD(3) * mrSges(6,2) + mrSges(6,3) * t45;
t19 = -qJDD(3) * mrSges(6,1) - t166;
t18 = t95 * qJD(1);
t16 = t46 + t64;
t15 = qJD(3) * t98 - t116;
t14 = qJD(5) * t78 + t118 * t80;
t13 = -qJD(5) * t80 + t118 * t78;
t10 = qJD(3) * t95 + t116;
t7 = qJD(3) * t82 + t215;
t4 = pkin(3) * t45 + t192;
t1 = t45 * t82 + qJDD(5) - t192;
t20 = [((t203 * t80 - t204 * t78) * t83 + t202 + (t200 / 0.2e1 - t99 / 0.2e1) * qJD(3)) * qJD(3) + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + t220 * t78 * t45 + (t109 + 0.2e1 * mrSges(3,3)) * t59 + t88 * t141 + m(4) * (-t198 * t83 + t122) + (-m(3) * t154 - m(4) * t136 + t182 * (t79 * t181 + t136) + t190 * t81 + t191 * t79) * g(2) + (t183 * t71 + t182 * (t81 * t181 + t71) + t191 * t81 + (m(3) * pkin(1) + (-m(4) + t182) * t83 - t190) * t79) * g(1) - t45 * t102 / 0.2e1 + t4 * t107 + t1 * t106 - t194 * t148 / 0.2e1 + t198 * mrSges(4,3) + (t213 * t80 + t219 * t78) * qJDD(3) / 0.2e1 + (-t148 * t25 - t197) * mrSges(5,2) + (t206 / 0.2e1 - t33 * mrSges(5,2) - t29 / 0.2e1) * t147 + ((-Ifges(4,4) + t212) * t45 + t222 * t44 + t195 * qJDD(3)) * t80 / 0.2e1 + (t193 + (t214 * t78 + t199) * t80) * t120 + t92 * t121 - pkin(1) * t143 + m(5) * (t15 * t26 + t4 * t49 + t83 * t85) + t34 * t124 + m(3) * (-pkin(1) * t67 + t122) + t199 * t45 / 0.2e1 + (t210 * qJDD(3) + t212 * t44) * t78 / 0.2e1 + m(6) * (t1 * t34 + t10 * t11 + t12 * t14 + t13 * t7 - t2 * t48 + t3 * t47) + t218 * t44 / 0.2e1 - qJDD(3) * (Ifges(6,5) * t80 + Ifges(6,6) * t78) / 0.2e1 - t78 * (Ifges(4,4) * t44 - Ifges(4,2) * t45 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t67 * mrSges(3,2) + qJ(2) * (mrSges(4,1) * t45 + mrSges(4,2) * t44) + t47 * t22 - t48 * t19 + t49 * (mrSges(5,1) * t45 - mrSges(5,3) * t44) + t14 * t50 + t13 * t52 + t15 * t40 + t10 * t41 + qJD(2) * t42 + (t148 * t7 + t209 * t80 + t177) * mrSges(6,3) + t207 * t78 * t83 + t208 * t80 * t83; t143 + (qJ(2) * t183 - mrSges(3,3)) * t189 + (-t19 + t208) * t80 + (t22 + t207) * t78 + (-m(5) * t26 + m(6) * t11 - t158 - t42) * qJD(1) + ((t51 + t205) * t80 + (t52 - t204) * t78) * qJD(3) + m(5) * t85 - m(4) * t198 + m(6) * (-t2 * t80 + t177 + (t12 * t80 + t7 * t78) * qJD(3)) + m(3) * t67 + t196 * (-t182 - t183); (t22 + t24) * qJ(4) + (t110 - (-m(6) * qJ(4) - mrSges(6,2)) * t78 - t115 * t80 + m(5) * t98 + t108) * t178 + (-pkin(3) * t6 - g(1) * t156 + qJ(4) * t5 + qJD(4) * t33 - t26 * t39 - t58 * t97) * m(5) + (-Ifges(4,6) + t210) * t45 - t202 * qJD(1) - t203 * t165 + t204 * t46 + t205 * qJD(4) + t200 * t121 + t194 * t150 / 0.2e1 + t195 * t44 + t29 * t149 / 0.2e1 + (-t88 + t92 / 0.2e1 - t193 / 0.2e1) * t189 + (Ifges(6,3) + Ifges(5,2) + Ifges(4,3)) * qJDD(3) - t110 * t179 + t99 * t120 - t7 * t132 - t12 * t130 + (-m(6) * t156 + (-t162 - (m(6) * pkin(4) + mrSges(6,1)) * t80 - t108) * t79) * g(1) + t82 * t19 - t17 * t50 - t16 * t52 - t39 * t40 - t18 * t41 - pkin(3) * t21 - t6 * mrSges(5,1) + t8 * mrSges(4,1) - t9 * mrSges(4,2) - t2 * mrSges(6,1) + t3 * mrSges(6,2) + t5 * mrSges(5,3) - (t214 * t150 + t206 + t217) * t149 / 0.2e1 + (qJ(4) * t3 - t11 * t18 + t215 * t12 - t16 * t7 + t2 * t82) * m(6) + (-m(5) * t216 - m(6) * t70 - t115 * t78 + t107 + t109 - t72) * g(3) + t33 * t131 + t25 * t133; -t166 + t36 + (-mrSges(5,1) - mrSges(6,1)) * qJDD(3) - t205 * qJD(3) + t182 * t78 * g(3) + t158 * t149 + (-t11 * t149 - t201 - t209) * m(6) + (-qJD(3) * t33 + t149 * t26 - t201 + t6) * m(5); (-t78 * t50 + t80 * t52) * qJD(1) + (t1 + g(1) * t81 + g(2) * t79 - (t12 * t78 - t7 * t80) * qJD(1)) * m(6) + t124;];
tau = t20;
