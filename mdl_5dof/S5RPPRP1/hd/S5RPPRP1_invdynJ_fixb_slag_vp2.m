% Calculate vector of inverse dynamics joint torques for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:09
% EndTime: 2020-01-03 11:25:23
% DurationCPUTime: 5.10s
% Computational Cost: add. (1469->288), mult. (3110->388), div. (0->0), fcn. (1845->10), ass. (0->147)
t215 = Ifges(5,4) + Ifges(6,4);
t228 = Ifges(5,1) + Ifges(6,1);
t227 = Ifges(5,2) + Ifges(6,2);
t100 = cos(qJ(4));
t230 = t215 * t100;
t98 = sin(qJ(4));
t229 = t215 * t98;
t214 = Ifges(5,5) + Ifges(6,5);
t213 = Ifges(5,6) + Ifges(6,6);
t226 = Ifges(5,3) + Ifges(6,3);
t93 = sin(pkin(8));
t225 = (-t227 * t98 + t230) * t93;
t224 = (t228 * t100 - t229) * t93;
t90 = t93 ^ 2;
t95 = cos(pkin(8));
t221 = mrSges(4,3) * (t95 ^ 2 + t90);
t94 = sin(pkin(7));
t78 = pkin(1) * t94 + qJ(3);
t65 = qJD(1) * qJD(3) + qJDD(1) * t78;
t76 = -qJD(1) * t95 + qJD(4);
t220 = t100 * (t225 * qJD(1) + t213 * t76) + t98 * (t224 * qJD(1) + t214 * t76);
t68 = t78 * qJD(1);
t85 = t95 * qJD(2);
t44 = t68 * t93 - t85;
t186 = t44 * t93;
t158 = t98 * qJD(1);
t21 = qJD(5) - t85 + (pkin(4) * t158 + t68) * t93;
t189 = t21 * t93;
t219 = (mrSges(5,1) * t100 - mrSges(5,2) * t98) * t186 + (mrSges(6,1) * t100 - mrSges(6,2) * t98) * t189 + (-t213 * t100 - t214 * t98) * t76 * t93 / 0.2e1;
t196 = m(6) * pkin(4);
t218 = -mrSges(6,1) - mrSges(5,1);
t217 = mrSges(3,2) - mrSges(4,3);
t216 = mrSges(5,2) + mrSges(6,2);
t123 = pkin(3) * t95 + pkin(6) * t93;
t96 = cos(pkin(7));
t82 = -pkin(1) * t96 - pkin(2);
t63 = -t123 + t82;
t166 = t100 * t95;
t66 = t78 * t166;
t20 = t98 * t63 + t66;
t211 = t196 + mrSges(6,1);
t152 = qJD(1) * qJD(4);
t60 = (qJDD(1) * t100 - t152 * t98) * t93;
t61 = (-qJDD(1) * t98 - t100 * t152) * t93;
t155 = qJDD(1) * t95;
t75 = qJDD(4) - t155;
t208 = t213 * t61 + t214 * t60 + t226 * t75;
t36 = -t95 * qJDD(2) + t65 * t93;
t37 = qJDD(2) * t93 + t65 * t95;
t207 = t36 * t93 + t37 * t95;
t154 = m(4) + m(5) + m(6);
t205 = -mrSges(5,1) - t211;
t121 = mrSges(4,1) * t95 - mrSges(4,2) * t93;
t201 = -mrSges(3,1) - m(6) * (pkin(4) * t100 + pkin(3)) * t95 - m(5) * t123 - t121 + (-m(6) * (qJ(5) + pkin(6)) - mrSges(6,3) - mrSges(5,3)) * t93;
t163 = qJD(1) * t93;
t140 = t100 * t163;
t129 = mrSges(6,3) * t140;
t55 = mrSges(6,1) * t76 - t129;
t130 = mrSges(5,3) * t140;
t56 = mrSges(5,1) * t76 - t130;
t173 = t55 + t56;
t145 = t93 * t158;
t131 = mrSges(6,3) * t145;
t53 = -mrSges(6,2) * t76 - t131;
t133 = mrSges(5,3) * t145;
t54 = -mrSges(5,2) * t76 - t133;
t174 = t53 + t54;
t200 = t174 * t100 - t173 * t98;
t199 = (-(-t228 * t98 - t230) * t100 / 0.2e1 + (-t227 * t100 - t229) * t98 / 0.2e1) * t90;
t160 = qJD(5) * t93;
t136 = qJD(1) * t160;
t157 = qJD(4) * t100;
t35 = qJDD(1) * t63 + qJDD(3);
t38 = qJD(1) * t63 + qJD(3);
t150 = t100 * t37 + t38 * t157 + t98 * t35;
t45 = qJD(2) * t93 + t68 * t95;
t2 = qJ(5) * t61 + (-qJD(4) * t45 - t136) * t98 + t150;
t161 = qJD(4) * t98;
t3 = -t161 * t45 + t150;
t11 = t100 * t45 + t38 * t98;
t4 = -qJD(4) * t11 + t100 * t35 - t37 * t98;
t198 = -t4 * mrSges(5,1) + t3 * mrSges(5,2) + t2 * mrSges(6,2);
t197 = qJD(1) ^ 2;
t195 = t60 / 0.2e1;
t194 = t61 / 0.2e1;
t193 = t75 / 0.2e1;
t99 = sin(qJ(1));
t88 = t99 * pkin(1);
t101 = cos(qJ(1));
t89 = t101 * pkin(1);
t92 = qJ(1) + pkin(7);
t86 = sin(t92);
t184 = t86 * t98;
t87 = cos(t92);
t183 = t87 * t98;
t182 = t93 * t98;
t181 = t95 * t98;
t22 = mrSges(6,1) * t75 - mrSges(6,3) * t60;
t23 = mrSges(5,1) * t75 - mrSges(5,3) * t60;
t177 = t22 + t23;
t24 = -mrSges(6,2) * t75 + mrSges(6,3) * t61;
t25 = -mrSges(5,2) * t75 + mrSges(5,3) * t61;
t176 = t24 + t25;
t162 = qJD(3) * t95;
t175 = t100 * t162 + t63 * t157;
t172 = mrSges(6,2) * t100;
t169 = qJ(5) * t93;
t168 = t100 * t93;
t165 = t11 * t100;
t164 = t86 * t100;
t156 = qJDD(1) * t93;
t151 = t78 * t181;
t146 = qJ(5) * t168;
t144 = t98 * t162;
t143 = t93 * t161;
t139 = t93 * t157;
t17 = -t61 * mrSges(6,1) + t60 * mrSges(6,2);
t10 = t100 * t38 - t45 * t98;
t122 = -mrSges(4,1) * t155 + mrSges(4,2) * t156;
t120 = t45 * t95 + t186;
t119 = mrSges(5,1) * t98 + mrSges(5,2) * t100;
t41 = t181 * t87 - t164;
t39 = -t100 * t87 - t181 * t86;
t118 = (mrSges(6,1) * t98 + t172) * t93;
t8 = -qJ(5) * t140 + t10;
t67 = qJDD(1) * t82 + qJDD(3);
t64 = (pkin(4) * t157 + qJD(3)) * t93;
t62 = t119 * t93;
t59 = (pkin(4) * t98 + t78) * t93;
t58 = t119 * t163;
t57 = qJD(1) * t118;
t52 = t100 * t63;
t42 = t166 * t87 + t184;
t40 = t164 * t95 - t183;
t19 = t52 - t151;
t18 = -mrSges(5,1) * t61 + mrSges(5,2) * t60;
t16 = -t169 * t98 + t20;
t15 = -qJD(4) * t20 - t144;
t14 = -qJD(4) * t151 + t175;
t13 = -t146 + t52 + (-t78 * t98 - pkin(4)) * t95;
t12 = -pkin(4) * t61 + qJDD(5) + t36;
t9 = -qJ(5) * t145 + t11;
t7 = -t144 - t100 * t160 + (-t66 + (-t63 + t169) * t98) * qJD(4);
t6 = -t98 * t160 + (-t146 - t151) * qJD(4) + t175;
t5 = pkin(4) * t76 + t8;
t1 = pkin(4) * t75 - qJ(5) * t60 - t100 * t136 + t4;
t26 = [(t219 - t220 * t93 / 0.2e1) * qJD(4) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t96 * mrSges(3,1) - 0.2e1 * t94 * mrSges(3,2) + m(3) * (t94 ^ 2 + t96 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + m(5) * (t10 * t15 + t11 * t14 + t19 * t4 + t20 * t3) + (-t1 * t168 - t139 * t9 + t143 * t5 - t182 * t2) * mrSges(6,3) + (Ifges(4,4) * t155 + Ifges(4,1) * t156 + t78 * t18 + qJD(3) * t58 + m(5) * (qJD(3) * t44 + t36 * t78) + (t100 * t214 - t213 * t98) * t193) * t93 + m(4) * (t120 * qJD(3) + t207 * t78 + t67 * t82) + t207 * mrSges(4,3) - t199 * t152 + t12 * t118 + (t10 * t143 - t11 * t139 - t168 * t4 - t182 * t3) * mrSges(5,3) + (-t184 * t196 - m(3) * t89 - mrSges(2,1) * t101 + mrSges(2,2) * t99 + t217 * t86 + t218 * t42 + t216 * t41 - t154 * (t87 * pkin(2) + t86 * qJ(3) + t89) + t201 * t87) * g(2) + (t183 * t196 - m(3) * t88 - mrSges(2,1) * t99 - mrSges(2,2) * t101 - t217 * t87 + t218 * t40 - t216 * t39 - t154 * (t86 * pkin(2) - qJ(3) * t87 + t88) + t201 * t86) * g(3) + m(6) * (t1 * t13 + t12 * t59 + t16 * t2 + t21 * t64 + t5 * t7 + t6 * t9) + t59 * t17 + t36 * t62 + t64 * t57 + t6 * t53 + t14 * t54 + t7 * t55 + t15 * t56 + t13 * t22 + t19 * t23 + t16 * t24 + t20 * t25 + t65 * t221 - t67 * t121 + t82 * t122 + (-t1 * mrSges(6,1) + Ifges(4,4) * t156 + Ifges(4,2) * t155 - t214 * t195 - t213 * t194 - t226 * t193 + t198 - t208 / 0.2e1) * t95 - (t213 * t75 + t215 * t60 + t227 * t61) * t182 / 0.2e1 + (t214 * t75 + t215 * t61 + t228 * t60) * t168 / 0.2e1 + t224 * t195 + t225 * t194; m(3) * qJDD(2) + (-m(3) - t154) * g(1) + (-m(5) * t36 - m(6) * t12 - t17 - t18) * t95 + (-t177 * t98 + t176 * t100 + (-t100 * t173 - t174 * t98) * qJD(4) + m(5) * (-t10 * t157 + t100 * t3 - t11 * t161 - t4 * t98) + m(6) * (-t1 * t98 + t100 * t2 - t157 * t5 - t161 * t9)) * t93 + (-t36 * t95 + t37 * t93) * m(4); t176 * t98 + t177 * t100 - t197 * t221 + t200 * qJD(4) + m(6) * (t1 * t100 + t2 * t98 + (t9 * t100 - t5 * t98) * qJD(4)) + m(5) * (t100 * t4 + t3 * t98 + (-t10 * t98 + t165) * qJD(4)) + m(4) * t67 + ((-t57 - t58) * t93 - t200 * t95 - m(5) * (-t10 * t181 + t165 * t95 + t186) - m(6) * (t166 * t9 - t181 * t5 + t189) - m(4) * t120) * qJD(1) + t122 + (g(2) * t87 + g(3) * t86) * t154; -t198 + (-m(6) * (-t5 + t8) + t55 + t129) * t9 + t199 * t197 + (t56 + t130) * t11 + (-t133 - t54) * t10 + t211 * t1 + (-(-t211 * t98 - t172) * t93 + t62) * g(1) - t8 * t53 - t5 * t131 + (t205 * t39 + t216 * t40) * g(2) + (t205 * t41 - t216 * t42) * g(3) + t208 + t220 * t163 / 0.2e1 - t219 * qJD(1) + ((-m(6) * t21 - t57) * t140 + t22) * pkin(4); (g(1) * t95 + t12) * m(6) + ((-g(2) * t86 + g(3) * t87) * m(6) + (t98 * t53 + t100 * t55 - m(6) * (-t100 * t5 - t9 * t98)) * qJD(1)) * t93 + t17;];
tau = t26;
