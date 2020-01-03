% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:21
% EndTime: 2019-12-31 20:03:29
% DurationCPUTime: 3.67s
% Computational Cost: add. (1915->299), mult. (4710->384), div. (0->0), fcn. (2645->4), ass. (0->131)
t126 = sin(qJ(2));
t155 = qJD(1) * t126;
t218 = pkin(6) * t155 + qJD(3);
t167 = Ifges(5,4) + Ifges(6,4);
t217 = -pkin(7) * t155 + t218;
t125 = sin(qJ(4));
t127 = cos(qJ(4));
t128 = cos(qJ(2));
t154 = qJD(1) * t128;
t66 = -t125 * t154 + t127 * t155;
t179 = t66 / 0.2e1;
t197 = Ifges(5,1) + Ifges(6,1);
t166 = Ifges(6,5) + Ifges(5,5);
t196 = Ifges(5,2) + Ifges(6,2);
t165 = Ifges(6,6) + Ifges(5,6);
t216 = -mrSges(3,1) - mrSges(4,1);
t114 = pkin(6) * t154;
t87 = -pkin(7) * t154 + t114;
t129 = -pkin(2) - pkin(3);
t91 = -qJ(3) * t125 + t127 * t129;
t210 = t91 * qJD(4) - t125 * t87 + t217 * t127;
t92 = t127 * qJ(3) + t125 * t129;
t208 = -qJD(4) * t92 - t217 * t125 - t127 * t87;
t65 = -t125 * t155 - t127 * t154;
t215 = t167 * t65;
t214 = t167 * t66;
t134 = t125 * t128 - t126 * t127;
t189 = qJD(2) - qJD(4);
t33 = t189 * t134;
t27 = t33 * qJD(1);
t151 = qJD(4) * t127;
t152 = qJD(4) * t125;
t146 = t129 * qJD(2);
t49 = t146 + t217;
t122 = qJD(2) * qJD(3);
t153 = qJD(2) * t126;
t178 = pkin(6) - pkin(7);
t88 = t178 * t153;
t52 = -qJD(1) * t88 + t122;
t123 = qJD(2) * qJ(3);
t67 = t123 + t87;
t101 = t178 * t128;
t89 = qJD(2) * t101;
t75 = qJD(1) * t89;
t5 = t125 * t75 + t127 * t52 + t49 * t151 - t152 * t67;
t1 = -qJ(5) * t27 + qJD(5) * t65 + t5;
t182 = -t65 / 0.2e1;
t194 = -t166 * t189 + t197 * t66 + t215;
t195 = -t165 * t189 + t196 * t65 + t214;
t133 = t125 * t126 + t127 * t128;
t34 = t189 * t133;
t26 = t34 * qJD(1);
t17 = t125 * t49 + t127 * t67;
t6 = -qJD(4) * t17 - t125 * t52 + t127 * t75;
t2 = -qJ(5) * t26 - qJD(5) * t66 + t6;
t68 = -qJD(1) * pkin(1) - pkin(2) * t154 - qJ(3) * t155;
t48 = pkin(3) * t154 - t68;
t24 = -pkin(4) * t65 + qJD(5) + t48;
t213 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t1 * mrSges(6,2) - t24 * (mrSges(6,1) * t66 + mrSges(6,2) * t65) - t48 * (mrSges(5,1) * t66 + mrSges(5,2) * t65) - t165 * t27 + t166 * t26 + (-t196 * t66 + t194 + t215) * t182 + (-t197 * t65 + t195 + t214) * t179;
t212 = -t154 / 0.2e1;
t159 = qJ(5) * t66;
t211 = -t159 + t210;
t160 = qJ(5) * t65;
t209 = -t160 + t208;
t207 = -t165 * t66 + t166 * t65;
t161 = Ifges(4,5) * t128;
t198 = -qJD(2) / 0.2e1;
t199 = -qJD(1) / 0.2e1;
t200 = Ifges(3,4) * t212;
t202 = -Ifges(3,1) / 0.2e1;
t93 = -qJD(2) * pkin(2) + t218;
t204 = (m(4) * t93 + (mrSges(4,2) + mrSges(3,3)) * t155 + t216 * qJD(2)) * pkin(6) - t68 * mrSges(4,3) - (t126 * Ifges(4,1) - t161) * t199 - t155 * t202 - t200 + t93 * mrSges(4,2) - (Ifges(4,4) + Ifges(3,5)) * t198;
t110 = Ifges(4,5) * t155;
t162 = Ifges(3,4) * t126;
t201 = Ifges(4,6) / 0.2e1;
t97 = t114 + t123;
t99 = mrSges(4,2) * t154 + qJD(2) * mrSges(4,3);
t203 = -(m(4) * t97 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t154 + t99) * pkin(6) + t68 * mrSges(4,1) + qJD(2) * t201 + Ifges(4,3) * t212 + t110 / 0.2e1 + Ifges(3,6) * t198 + (t128 * Ifges(3,2) + t162) * t199 - t97 * mrSges(4,2);
t40 = mrSges(6,2) * t189 + t65 * mrSges(6,3);
t41 = mrSges(5,2) * t189 + t65 * mrSges(5,3);
t42 = -mrSges(6,1) * t189 - t66 * mrSges(6,3);
t43 = -mrSges(5,1) * t189 - mrSges(5,3) * t66;
t186 = -t125 * (t42 + t43) + t127 * (t40 + t41);
t177 = pkin(1) * mrSges(3,1);
t176 = pkin(1) * mrSges(3,2);
t175 = pkin(4) * t66;
t173 = -t189 / 0.2e1;
t168 = t27 * t92;
t100 = t178 * t126;
t38 = t125 * t100 + t127 * t101;
t108 = qJ(3) * t154;
t117 = t126 * qJD(3);
t157 = qJD(1) * t117 + qJD(2) * t108;
t156 = t128 * t123 + t117;
t94 = -t128 * pkin(2) - t126 * qJ(3) - pkin(1);
t150 = Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t149 = 0.3e1 / 0.2e1 * Ifges(4,5) - 0.3e1 / 0.2e1 * Ifges(3,4);
t148 = t201 - Ifges(3,6) / 0.2e1;
t147 = m(4) * pkin(6) + mrSges(4,2);
t145 = qJD(1) * t153;
t16 = -t125 * t67 + t127 * t49;
t37 = t127 * t100 - t101 * t125;
t76 = t128 * pkin(3) - t94;
t144 = t126 * t146;
t10 = t17 + t160;
t9 = t16 - t159;
t8 = -pkin(4) * t189 + t9;
t143 = t10 * t66 + t65 * t8;
t138 = t16 * t65 + t17 * t66;
t55 = t129 * t155 + t108;
t11 = t100 * t151 - t101 * t152 + t125 * t89 - t127 * t88;
t44 = t144 + t156;
t36 = qJD(1) * t144 + t157;
t12 = -qJD(4) * t38 + t125 * t88 + t127 * t89;
t90 = -pkin(6) * t145 + t122;
t85 = -pkin(4) + t91;
t83 = (-t128 * mrSges(4,1) - mrSges(4,3) * t126) * qJD(1);
t60 = pkin(2) * t153 - t156;
t47 = pkin(2) * t145 - t157;
t35 = pkin(4) * t133 + t76;
t30 = t55 - t175;
t29 = -mrSges(5,1) * t65 + mrSges(5,2) * t66;
t28 = -mrSges(6,1) * t65 + mrSges(6,2) * t66;
t25 = t26 * mrSges(6,2);
t23 = -qJ(5) * t133 + t38;
t22 = qJ(5) * t134 + t37;
t13 = pkin(4) * t33 + t44;
t7 = pkin(4) * t27 + t36;
t4 = -qJ(5) * t34 + qJD(5) * t134 + t12;
t3 = -qJ(5) * t33 - qJD(5) * t133 + t11;
t14 = [-t195 * t33 / 0.2e1 + t194 * t34 / 0.2e1 + (-t167 * t33 + t197 * t34) * t179 + (-t1 * t133 - t10 * t33 + t134 * t2 - t8 * t34) * mrSges(6,3) + (-t133 * t5 + t134 * t6 - t16 * t34 - t17 * t33) * mrSges(5,3) + t48 * (mrSges(5,1) * t33 + mrSges(5,2) * t34) + t24 * (mrSges(6,1) * t33 + mrSges(6,2) * t34) - (-mrSges(5,1) * t76 - mrSges(6,1) * t35 + mrSges(5,3) * t38 + mrSges(6,3) * t23 - t133 * t196 - t134 * t167) * t27 + (mrSges(5,2) * t76 - mrSges(5,3) * t37 - mrSges(6,3) * t22 - t133 * t167 - t134 * t197) * t26 + t7 * (mrSges(6,1) * t133 - mrSges(6,2) * t134) + t36 * (mrSges(5,1) * t133 - mrSges(5,2) * t134) + (-t165 * t33 + t166 * t34) * t173 + (t167 * t34 - t196 * t33) * t65 / 0.2e1 + m(4) * (t47 * t94 + t60 * t68) + m(5) * (t11 * t17 + t12 * t16 + t36 * t76 + t37 * t6 + t38 * t5 + t44 * t48) + m(6) * (t1 * t23 + t10 * t3 + t13 * t24 + t2 * t22 + t35 * t7 + t4 * t8) + (-t47 * mrSges(4,3) + (t148 * qJD(2) + (t94 * mrSges(4,1) + t126 * t149 - 0.2e1 * t177) * qJD(1) + t203) * qJD(2)) * t126 + (-t47 * mrSges(4,1) + t147 * t90 + (t150 * qJD(2) + (-0.2e1 * t176 - t94 * mrSges(4,3) - t149 * t128 + (-0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(3,1) + t147 * pkin(6)) * t126) * qJD(1) + t204) * qJD(2)) * t128 + t60 * t83 + t3 * t40 + t11 * t41 + t4 * t42 + t12 * t43 + t44 * t29 + t13 * t28 + t35 * t25; ((t200 + (t176 + t161 / 0.2e1) * qJD(1) + (-pkin(2) * mrSges(4,2) + (-m(4) * pkin(2) + t216) * pkin(6) + t150) * qJD(2) - t204) * t128 + (-t110 / 0.2e1 + (t177 + t162 / 0.2e1) * qJD(1) + (-Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + t202) * t154 + (pkin(6) * mrSges(3,2) - qJ(3) * mrSges(4,2) + t148) * qJD(2) - t203) * t126) * qJD(1) + t208 * t43 + t209 * t42 + t210 * t41 + t211 * t40 + (-t26 * t85 - t143 - t168) * mrSges(6,3) + (-t26 * t91 - t138 - t168) * mrSges(5,3) + t90 * mrSges(4,3) + qJD(3) * t99 - t55 * t29 - t30 * t28 + m(4) * (qJ(3) * t90 + qJD(3) * t97) + (-m(4) * t68 - t83) * (pkin(2) * t155 - t108) + t207 * t173 + (t1 * t92 + t10 * t211 + t2 * t85 + t209 * t8 - t24 * t30) * m(6) + (t16 * t208 + t17 * t210 - t48 * t55 + t5 * t92 + t6 * t91) * m(5) - t213; (-t28 - t29 + t83) * t155 + t186 * qJD(4) + (t147 * t154 - t186 - t99) * qJD(2) - m(4) * (qJD(2) * t97 - t155 * t68) + (mrSges(6,3) + mrSges(5,3)) * (-t125 * t27 - t127 * t26) + (t1 * t125 + t127 * t2 - t24 * t155 - t189 * (t10 * t127 - t125 * t8)) * m(6) + (t125 * t5 + t127 * t6 - t48 * t155 - t189 * (-t125 * t16 + t127 * t17)) * m(5); (-t24 * t175 - (-t8 + t9) * t10 + t2 * pkin(4)) * m(6) - t9 * t40 - t16 * t41 + t10 * t42 + t17 * t43 - t28 * t175 + (-pkin(4) * t26 + t143) * mrSges(6,3) + t138 * mrSges(5,3) + t207 * t189 / 0.2e1 + t213; t27 * mrSges(6,1) - t65 * t40 + t66 * t42 + t25 + 0.2e1 * (t7 / 0.2e1 + t10 * t182 + t8 * t179) * m(6);];
tauc = t14(:);
