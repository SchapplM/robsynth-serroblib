% Calculate vector of inverse dynamics joint torques for
% S5RPRPP2
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:41
% EndTime: 2019-12-31 18:10:50
% DurationCPUTime: 5.61s
% Computational Cost: add. (1087->320), mult. (2189->380), div. (0->0), fcn. (1013->8), ass. (0->140)
t213 = Ifges(5,1) + Ifges(6,1);
t211 = Ifges(6,4) + Ifges(5,5);
t212 = Ifges(5,4) + Ifges(4,5);
t195 = -Ifges(6,5) + t212;
t222 = Ifges(6,2) + Ifges(5,3);
t221 = Ifges(4,6) - Ifges(5,6);
t209 = Ifges(5,6) - Ifges(6,6);
t97 = cos(qJ(3));
t169 = Ifges(5,5) * t97;
t171 = Ifges(6,4) * t97;
t95 = sin(qJ(3));
t220 = t213 * t95 - t169 - t171;
t89 = qJ(1) + pkin(7);
t77 = sin(t89);
t78 = cos(t89);
t196 = g(1) * t78 + g(2) * t77;
t83 = t95 * qJ(4);
t130 = pkin(2) + t83;
t94 = cos(pkin(7));
t182 = pkin(1) * t94;
t183 = pkin(3) + pkin(4);
t10 = qJD(5) + (t183 * t97 + t130 + t182) * qJD(1);
t121 = mrSges(4,1) * t95 + mrSges(4,2) * t97;
t137 = pkin(2) + t182;
t164 = t97 * mrSges(5,3);
t165 = t97 * mrSges(6,2);
t87 = t97 * pkin(3);
t112 = -t130 - t87;
t19 = (t112 - t182) * qJD(1);
t219 = -t137 * qJD(1) * t121 + t19 * (mrSges(5,1) * t95 - t164) + t10 * (-t95 * mrSges(6,1) + t165);
t184 = m(5) + m(6);
t218 = -m(4) - t184;
t153 = t95 * qJD(1);
t217 = t211 * t153;
t161 = t87 + t83;
t150 = qJD(1) * qJD(3);
t51 = qJDD(1) * t95 + t150 * t97;
t53 = t137 * qJDD(1);
t216 = -qJ(4) * t51 - qJD(4) * t153 - t53;
t93 = sin(pkin(7));
t70 = pkin(1) * t93 + pkin(6);
t215 = qJD(2) * qJD(3) + t70 * qJDD(1);
t120 = t97 * mrSges(5,1) + t95 * mrSges(5,3);
t122 = mrSges(4,1) * t97 - mrSges(4,2) * t95;
t214 = -t120 - t122;
t42 = t51 * mrSges(5,2);
t28 = -qJDD(3) * mrSges(5,1) + t42;
t208 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t51 + t28;
t50 = -qJDD(1) * t97 + t150 * t95;
t29 = -mrSges(5,2) * t50 + qJDD(3) * mrSges(5,3);
t207 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t50 + t29;
t152 = t97 * qJD(1);
t206 = qJD(3) * t209 - t152 * t222 + t217;
t143 = mrSges(4,3) * t153;
t145 = mrSges(5,2) * t153;
t205 = t143 + t145 + (-mrSges(4,1) - mrSges(5,1)) * qJD(3);
t138 = mrSges(6,3) * t152;
t59 = qJD(3) * mrSges(6,2) - t138;
t144 = mrSges(5,2) * t152;
t61 = qJD(3) * mrSges(5,3) + t144;
t204 = t59 + t61;
t142 = mrSges(4,3) * t152;
t60 = -qJD(3) * mrSges(4,2) + t142;
t203 = -t61 - t60;
t84 = t95 * mrSges(6,2);
t162 = mrSges(6,1) * t97 + t84;
t201 = t212 * t97 - t221 * t95;
t151 = qJ(5) * qJD(1);
t54 = t70 * qJD(1);
t21 = qJD(2) * t97 - t54 * t95;
t12 = t151 * t95 + t21;
t200 = -t12 + qJD(4);
t199 = t211 * t95;
t156 = qJD(3) * t95;
t7 = qJDD(2) * t95 - t156 * t54 + t215 * t97;
t155 = qJD(3) * t97;
t8 = qJDD(2) * t97 - t54 * t155 - t215 * t95;
t198 = t7 * t97 - t8 * t95;
t4 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t7;
t103 = qJDD(4) - t8;
t5 = -qJDD(3) * pkin(3) + t103;
t197 = t4 * t97 + t5 * t95;
t76 = Ifges(4,4) * t152;
t194 = Ifges(4,1) * t153 + t220 * qJD(1) + qJD(3) * t195 + t76;
t193 = -mrSges(3,1) - t84 + t214;
t160 = qJ(4) * t97;
t192 = t196 * t160;
t191 = m(6) * qJ(5) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t96 = sin(qJ(1));
t181 = pkin(1) * t96;
t98 = cos(qJ(1));
t88 = t98 * pkin(1);
t174 = Ifges(4,4) * t95;
t173 = Ifges(4,4) * t97;
t168 = t51 * mrSges(6,3);
t22 = qJD(2) * t95 + t54 * t97;
t158 = qJ(5) - t70;
t154 = qJD(4) * t95;
t149 = qJD(1) * qJD(5);
t141 = mrSges(6,3) * t153;
t131 = -mrSges(6,1) * t50 + mrSges(6,2) * t51;
t39 = t158 * t97;
t129 = -t150 / 0.2e1;
t128 = t150 / 0.2e1;
t127 = -m(6) * t183 - mrSges(6,1);
t13 = -t151 * t97 + t22;
t92 = qJD(3) * qJ(4);
t11 = t13 + t92;
t9 = -qJD(3) * t183 + t200;
t123 = t11 * t97 + t9 * t95;
t117 = t97 * Ifges(4,2) + t174;
t114 = Ifges(6,5) * t97 + Ifges(6,6) * t95;
t113 = pkin(3) * t95 - t160;
t37 = -t137 - t161;
t111 = -t183 * t95 + t160;
t107 = t95 * (Ifges(4,1) * t97 - t174);
t106 = t97 * (Ifges(6,2) * t95 + t171);
t105 = t97 * (Ifges(5,3) * t95 + t169);
t86 = t97 * pkin(4);
t55 = -qJD(3) * mrSges(6,1) - t141;
t49 = t113 * qJD(1);
t48 = t162 * qJD(1);
t47 = t120 * qJD(1);
t38 = t158 * t95;
t33 = Ifges(4,6) * qJD(3) + qJD(1) * t117;
t30 = qJD(3) * t113 - t154;
t26 = -qJDD(3) * mrSges(6,1) - t168;
t24 = qJDD(3) * mrSges(6,2) + mrSges(6,3) * t50;
t23 = t111 * qJD(1);
t20 = -t37 + t86;
t18 = t92 + t22;
t17 = -qJD(3) * t39 - qJD(5) * t95;
t16 = -qJD(5) * t97 + t156 * t158;
t15 = qJD(3) * t111 + t154;
t14 = -qJD(3) * pkin(3) + qJD(4) - t21;
t6 = pkin(3) * t50 + t216;
t3 = -t183 * t50 + qJDD(5) - t216;
t2 = qJ(5) * t50 - t149 * t97 + t4;
t1 = -qJ(5) * t51 - qJDD(3) * t183 - t149 * t95 + t103;
t25 = [m(5) * (t19 * t30 + t37 * t6) + (t212 * t95 + t221 * t97) * qJDD(3) / 0.2e1 + (0.2e1 * (mrSges(3,1) * t94 - mrSges(3,2) * t93) * pkin(1) + Ifges(3,3) + Ifges(2,3) + m(3) * (t93 ^ 2 + t94 ^ 2) * pkin(1) ^ 2) * qJDD(1) + t194 * t155 / 0.2e1 + ((Ifges(4,1) + t213) * t51 + (-Ifges(4,4) + t211) * t50 + t195 * qJDD(3)) * t95 / 0.2e1 + (t201 / 0.2e1 - t114 / 0.2e1) * qJD(3) ^ 2 - t6 * t120 - t50 * t117 / 0.2e1 + (t206 / 0.2e1 - t22 * mrSges(4,3) - t18 * mrSges(5,2) - t33 / 0.2e1 + t11 * mrSges(6,3)) * t156 + t97 * (Ifges(4,4) * t51 - Ifges(4,2) * t50 + Ifges(4,6) * qJDD(3)) / 0.2e1 - qJDD(3) * (Ifges(6,5) * t95 - Ifges(6,6) * t97) / 0.2e1 + t53 * t122 + (m(4) * t53 - mrSges(4,1) * t50 - mrSges(4,2) * t51) * t137 + t219 * qJD(3) + ((t213 * t97 + t199) * t95 + t97 * (-Ifges(4,2) * t95 + t173) + t107) * t128 + (m(4) * ((-t21 * t97 - t22 * t95) * qJD(3) + t198) + m(5) * ((t14 * t97 - t18 * t95) * qJD(3) + t197) + t203 * t156 + t205 * t155 + t207 * t97 + t208 * t95) * t70 + t3 * t162 + (-t155 * t21 + t198) * mrSges(4,3) + (t14 * t155 + t197) * mrSges(5,2) + t20 * t131 + t37 * (mrSges(5,1) * t50 - mrSges(5,3) * t51) + t17 * t55 + t16 * t59 - t30 * t47 + t15 * t48 + (t106 + t105) * t129 - t38 * t26 - t39 * t24 + (t95 * Ifges(4,1) + t173 + t220) * t51 / 0.2e1 + t199 * t50 / 0.2e1 - (qJDD(3) * t209 + t211 * t51) * t97 / 0.2e1 - t222 * t97 * t50 + m(6) * (-t1 * t38 + t10 * t15 + t11 * t16 + t17 * t9 - t2 * t39 + t20 * t3) + (-t1 * t95 - t155 * t9 - t2 * t97) * mrSges(6,3) + (m(3) * t181 + mrSges(2,1) * t96 + mrSges(2,2) * t98 + t218 * (pkin(6) * t78 - t181) + t191 * t78 + (m(4) * pkin(2) - m(5) * t112 + m(6) * t130 - t127 * t97 - t193) * t77) * g(1) + (-m(3) * t88 - mrSges(2,1) * t98 + mrSges(2,2) * t96 + t191 * t77 + t218 * (pkin(2) * t78 + pkin(6) * t77 + t88) + (-t184 * t161 - (m(6) * pkin(4) + mrSges(6,1)) * t97 + t193) * t78) * g(2); m(3) * qJDD(2) + (-t26 - t208) * t97 + (t24 + t207) * t95 + ((t60 + t204) * t97 + (t55 + t205) * t95) * qJD(3) + m(4) * (t7 * t95 + t8 * t97 + (-t21 * t95 + t22 * t97) * qJD(3)) + m(5) * (t4 * t95 - t5 * t97 + (t14 * t95 + t18 * t97) * qJD(3)) + m(6) * (qJD(3) * t123 - t1 * t97 + t2 * t95) + (-m(3) + t218) * g(3); (-m(5) * t18 + t142 + t203) * t21 + t204 * qJD(4) + (-m(5) * t14 + t143 - t205) * t22 + (-Ifges(4,6) + t209) * t50 + t201 * t129 - (-Ifges(4,2) * t153 + t194 + t76) * t152 / 0.2e1 + t195 * t51 + (t29 + t24) * qJ(4) + (Ifges(6,3) + Ifges(5,2) + Ifges(4,3)) * qJDD(3) + (-pkin(3) * t5 + qJ(4) * t4 + qJD(4) * t18 - t19 * t49 - t192) * m(5) + (qJ(4) * t2 - t1 * t183 - t10 * t23 + t11 * t200 - t13 * t9 - t192) * m(6) + (t121 - t164 - t165 + (m(5) * pkin(3) + mrSges(5,1) - t127) * t95) * t196 + t18 * t145 + t9 * t138 - t183 * t26 + ((-t107 / 0.2e1 + t106 / 0.2e1 + t105 / 0.2e1) * qJD(1) - t219) * qJD(1) + t33 * t153 / 0.2e1 - t14 * t144 - t11 * t141 - t13 * t55 - t12 * t59 - t23 * t48 + t49 * t47 - pkin(3) * t28 + t2 * mrSges(6,2) + t4 * mrSges(5,3) - t5 * mrSges(5,1) - t7 * mrSges(4,2) + t8 * mrSges(4,1) - t1 * mrSges(6,1) + (-m(5) * t161 - m(6) * (t86 + t161) - t162 + t214) * g(3) - (t213 * t152 + t206 + t217) * t153 / 0.2e1 + t114 * t128; -t168 + t42 + (-mrSges(5,1) - mrSges(6,1)) * qJDD(3) - t204 * qJD(3) + t184 * t97 * g(3) + ((-t47 - t48) * qJD(1) - t184 * t196) * t95 + (-qJD(3) * t11 - t10 * t153 + t1) * m(6) + (-qJD(3) * t18 + t153 * t19 + t5) * m(5); (t55 * t95 + t59 * t97) * qJD(1) + (g(1) * t77 - g(2) * t78 + qJD(1) * t123 + t3) * m(6) + t131;];
tau = t25;
