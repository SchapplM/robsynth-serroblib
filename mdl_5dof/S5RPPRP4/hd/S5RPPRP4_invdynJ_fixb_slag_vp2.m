% Calculate vector of inverse dynamics joint torques for
% S5RPPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:03
% EndTime: 2019-12-31 17:52:09
% DurationCPUTime: 3.26s
% Computational Cost: add. (1178->269), mult. (2163->340), div. (0->0), fcn. (1051->6), ass. (0->123)
t188 = Ifges(5,1) + Ifges(6,1);
t187 = Ifges(5,2) + Ifges(6,2);
t178 = Ifges(6,5) + Ifges(5,5);
t177 = Ifges(6,6) + Ifges(5,6);
t86 = sin(qJ(4));
t152 = Ifges(6,4) * t86;
t154 = Ifges(5,4) * t86;
t87 = cos(qJ(4));
t186 = -t187 * t87 - t152 - t154;
t151 = Ifges(6,4) * t87;
t153 = Ifges(5,4) * t87;
t185 = -t188 * t86 - t151 - t153;
t163 = -t86 / 0.2e1;
t109 = mrSges(6,1) * t87 - mrSges(6,2) * t86;
t56 = mrSges(5,1) * t87 - mrSges(5,2) * t86;
t184 = t109 + t56;
t88 = -pkin(1) - pkin(2);
t57 = qJDD(1) * t88 + qJDD(2);
t128 = qJD(1) * qJD(2);
t59 = qJDD(1) * qJ(2) + t128;
t82 = sin(pkin(7));
t83 = cos(pkin(7));
t20 = t82 * t57 + t83 * t59;
t18 = -qJDD(1) * pkin(6) + t20;
t183 = qJD(3) * qJD(4) + t18;
t137 = qJD(4) * t86;
t133 = qJ(2) * qJD(1);
t58 = qJD(1) * t88 + qJD(2);
t26 = t83 * t133 + t82 * t58;
t22 = -qJD(1) * pkin(6) + t26;
t3 = t86 * qJDD(3) - t137 * t22 + t183 * t87;
t13 = qJD(3) * t86 + t22 * t87;
t74 = t87 * qJDD(3);
t4 = -qJD(4) * t13 - t18 * t86 + t74;
t111 = t3 * t87 - t4 * t86;
t12 = t87 * qJD(3) - t22 * t86;
t136 = qJD(4) * t87;
t182 = -t12 * t136 - t13 * t137 + t111;
t181 = mrSges(3,1) + mrSges(2,1);
t180 = -mrSges(3,3) + mrSges(2,2);
t179 = Ifges(5,4) + Ifges(6,4);
t176 = t186 * qJD(1) + t177 * qJD(4);
t175 = t185 * qJD(1) + t178 * qJD(4);
t134 = t87 * qJD(1);
t25 = -t82 * t133 + t58 * t83;
t21 = qJD(1) * pkin(3) - t25;
t14 = pkin(4) * t134 + qJD(5) + t21;
t174 = t21 * (-mrSges(5,1) * t86 - mrSges(5,2) * t87) + t14 * (-mrSges(6,1) * t86 - mrSges(6,2) * t87);
t173 = t177 * t86 - t178 * t87;
t171 = -m(6) - m(5) - m(4);
t78 = t87 * pkin(4);
t170 = m(6) * (t78 + pkin(3)) + m(5) * pkin(3) + mrSges(4,1) + t184;
t126 = qJD(1) * qJD(5);
t127 = qJD(1) * qJD(4);
t45 = -qJDD(1) * t86 - t127 * t87;
t1 = -t22 * t136 + qJDD(4) * pkin(4) - qJ(5) * t45 + t74 + (t126 - t183) * t86;
t132 = qJ(5) * qJD(1);
t10 = -t132 * t87 + t13;
t44 = -qJDD(1) * t87 + t127 * t86;
t2 = qJ(5) * t44 - t126 * t87 + t3;
t9 = t86 * t132 + t12;
t6 = qJD(4) * pkin(4) + t9;
t168 = t1 * t86 + t10 * t137 + t6 * t136 - t2 * t87;
t54 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t134;
t55 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t134;
t143 = t54 + t55;
t135 = t86 * qJD(1);
t52 = qJD(4) * mrSges(6,1) + mrSges(6,3) * t135;
t53 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t135;
t144 = t52 + t53;
t167 = t143 * t87 - t144 * t86;
t166 = -m(5) * pkin(6) + m(6) * (-qJ(5) - pkin(6)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t156 = cos(qJ(1));
t155 = sin(qJ(1));
t150 = t21 * t82;
t149 = t83 * t86;
t148 = t83 * t87;
t27 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t44;
t28 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t44;
t146 = t27 + t28;
t29 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t45;
t30 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t45;
t145 = t29 + t30;
t47 = t83 * qJ(2) + t82 * t88;
t142 = t156 * pkin(1) + t155 * qJ(2);
t40 = -pkin(6) + t47;
t140 = qJ(5) - t40;
t139 = qJD(2) * t82;
t138 = qJD(2) * t83;
t131 = qJDD(1) * mrSges(3,1);
t130 = qJDD(1) * mrSges(4,1);
t129 = qJDD(1) * mrSges(4,2);
t15 = -t44 * mrSges(6,1) + t45 * mrSges(6,2);
t19 = t57 * t83 - t82 * t59;
t46 = -t82 * qJ(2) + t83 * t88;
t114 = qJD(4) * t140;
t113 = -qJD(5) + t138;
t39 = pkin(3) - t46;
t17 = qJDD(1) * pkin(3) - t19;
t112 = -pkin(1) * t155 + t156 * qJ(2);
t110 = t10 * t87 - t6 * t86;
t102 = -t12 * t87 - t13 * t86;
t101 = -t12 * t86 + t13 * t87;
t100 = -t25 * t82 + t26 * t83;
t97 = t86 * (-Ifges(5,1) * t87 + t154);
t96 = t86 * (-Ifges(6,1) * t87 + t152);
t95 = t87 * (Ifges(5,2) * t86 - t153);
t94 = t87 * (Ifges(6,2) * t86 - t151);
t89 = qJD(1) ^ 2;
t70 = -qJDD(1) * pkin(1) + qJDD(2);
t49 = -pkin(4) * t137 + t139;
t42 = t56 * qJD(1);
t41 = t109 * qJD(1);
t38 = -t155 * t83 + t156 * t82;
t37 = -t155 * t82 - t156 * t83;
t31 = t39 + t78;
t24 = t140 * t87;
t23 = t140 * t86;
t16 = -mrSges(5,1) * t44 + mrSges(5,2) * t45;
t8 = -t113 * t86 + t114 * t87;
t7 = t113 * t87 + t114 * t86;
t5 = -pkin(4) * t44 + qJDD(5) + t17;
t11 = [-(t97 + t96 + t95 + t94) * t127 / 0.2e1 + m(6) * (t1 * t23 + t10 * t7 + t14 * t49 - t2 * t24 + t31 * t5 + t6 * t8) + m(4) * (qJD(2) * t100 + t19 * t46 + t20 * t47) + (t128 * t82 - t19) * mrSges(4,1) + (Ifges(4,3) + Ifges(3,2) + Ifges(2,3)) * qJDD(1) + (t128 * t83 + t20) * mrSges(4,2) + t168 * mrSges(6,3) + t5 * t109 + m(3) * (-pkin(1) * t70 + (t59 + t128) * qJ(2)) - t46 * t130 + (t179 * t44 + t188 * t45) * t163 - (t179 * t45 + t187 * t44) * t87 / 0.2e1 - t70 * mrSges(3,1) + t49 * t41 + t8 * t52 + t7 * t54 + t17 * t56 + 0.2e1 * t59 * mrSges(3,3) + t39 * t16 - t24 * t27 + t23 * t29 + t31 * t15 + (-t86 * t53 + t87 * t55) * t138 + (-m(3) * t112 + t180 * t156 + t181 * t155 + t171 * (-pkin(2) * t155 + t112) - t170 * t38 + t166 * t37) * g(1) + (-m(3) * t142 - t181 * t156 + t180 * t155 + t171 * (t156 * pkin(2) + t142) + t166 * t38 + t170 * t37) * g(2) + t186 * t44 / 0.2e1 + t185 * t45 / 0.2e1 + m(5) * (t17 * t39 + (t101 * t83 + t150) * qJD(2)) + (-t53 * t136 - t55 * t137 + m(5) * (qJD(4) * t102 + t111) + t87 * t28 - t86 * t30) * t40 + t47 * t129 + pkin(1) * t131 + t42 * t139 + t173 * qJD(4) ^ 2 / 0.2e1 + t174 * qJD(4) - t175 * t136 / 0.2e1 + t176 * t137 / 0.2e1 - t182 * mrSges(5,3) + (0.2e1 * t178 * t163 - t177 * t87) * qJDD(4); -t131 - t89 * mrSges(3,3) + (-t89 * qJ(2) + t70) * m(3) + (m(4) * t19 - m(5) * t17 - m(6) * t5 - t89 * mrSges(4,2) - t130 - t15 - t16) * t83 + (-t89 * mrSges(4,1) + t129 + t146 * t87 - t145 * t86 + (-t143 * t86 - t144 * t87) * qJD(4) + m(5) * t182 - m(6) * t168 + m(4) * t20) * t82 + ((-t41 - t42) * t82 - t167 * t83 - m(6) * (t10 * t148 + t14 * t82 - t149 * t6) - m(5) * (-t12 * t149 + t13 * t148 + t150) - m(4) * t100) * qJD(1) + (-g(1) * t155 + g(2) * t156) * (m(3) - t171); m(4) * qJDD(3) + t145 * t87 + t146 * t86 + t167 * qJD(4) + m(5) * (qJD(4) * t101 + t3 * t86 + t4 * t87) + m(6) * (qJD(4) * t110 + t1 * t87 + t2 * t86) - t171 * g(3); t4 * mrSges(5,1) + t1 * mrSges(6,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2) - t12 * t55 + t13 * t53 - t9 * t54 + t178 * t45 + t177 * t44 + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) + t184 * g(3) + (t29 + (g(3) * t87 + t1) * m(6)) * pkin(4) + (-m(6) * (-t6 + t9) + t52) * t10 + ((t96 / 0.2e1 + t94 / 0.2e1 + t97 / 0.2e1 + t95 / 0.2e1) * qJD(1) + (m(6) * t14 + t41) * t86 * pkin(4) + (-t10 * t86 - t6 * t87) * mrSges(6,3) + t102 * mrSges(5,3) + t176 * t163 + t175 * t87 / 0.2e1 - t173 * qJD(4) / 0.2e1 - t174) * qJD(1) + (g(1) * t37 + g(2) * t38) * ((-mrSges(5,2) - mrSges(6,2)) * t87 + (-m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1)) * t86); (-t86 * t52 + t87 * t54) * qJD(1) + (-g(1) * t38 + g(2) * t37 + qJD(1) * t110 + t5) * m(6) + t15;];
tau = t11;
