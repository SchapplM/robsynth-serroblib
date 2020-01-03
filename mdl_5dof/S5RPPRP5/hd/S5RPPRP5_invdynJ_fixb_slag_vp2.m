% Calculate vector of inverse dynamics joint torques for
% S5RPPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:21
% EndTime: 2019-12-31 17:53:31
% DurationCPUTime: 7.06s
% Computational Cost: add. (1319->301), mult. (3174->377), div. (0->0), fcn. (2020->6), ass. (0->121)
t125 = qJD(1) * qJD(2);
t71 = qJDD(1) * qJ(2) + t125;
t146 = -mrSges(6,1) - mrSges(5,1);
t186 = Ifges(5,1) + Ifges(6,1);
t169 = Ifges(6,4) + Ifges(5,5);
t95 = sin(pkin(7));
t93 = t95 ^ 2;
t96 = cos(pkin(7));
t94 = t96 ^ 2;
t180 = t93 + t94;
t192 = (t180 + t94) * t71;
t184 = -Ifges(5,4) + Ifges(6,5);
t145 = -mrSges(5,2) + mrSges(6,3);
t171 = m(6) * qJ(5) + t145;
t191 = m(6) * pkin(4) - t146;
t100 = cos(qJ(4));
t138 = t100 * t95;
t98 = sin(qJ(4));
t59 = -t96 * t98 + t138;
t158 = t59 / 0.2e1;
t189 = -m(6) - m(5);
t187 = mrSges(4,2) + mrSges(3,3);
t185 = Ifges(3,4) - Ifges(4,5);
t168 = Ifges(6,6) - Ifges(5,6);
t58 = t100 * t96 + t95 * t98;
t50 = t58 * qJD(4);
t25 = -qJD(1) * t50 + qJDD(1) * t59;
t12 = -qJDD(4) * mrSges(6,1) + t25 * mrSges(6,2);
t183 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t25 - t12;
t129 = qJD(4) * t100;
t120 = t95 * t129;
t132 = qJD(1) * t96;
t122 = t98 * t132;
t26 = qJD(1) * t120 - qJD(4) * t122 + qJDD(1) * t58;
t14 = -mrSges(6,2) * t26 + qJDD(4) * mrSges(6,3);
t182 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t26 + t14;
t52 = t58 * qJD(1);
t153 = Ifges(6,5) * t52;
t48 = Ifges(5,4) * t52;
t133 = qJD(1) * t95;
t53 = t100 * t133 - t122;
t181 = t169 * qJD(4) + t186 * t53 + t153 - t48;
t156 = mrSges(5,3) * t52;
t38 = -mrSges(6,2) * t52 + qJD(4) * mrSges(6,3);
t143 = -qJD(4) * mrSges(5,2) - t156 + t38;
t155 = mrSges(5,3) * t53;
t142 = mrSges(6,2) * t53 + t146 * qJD(4) + t155;
t87 = t95 * qJ(3);
t115 = pkin(1) + t87;
t179 = -t96 * pkin(2) - t115;
t150 = t95 * mrSges(4,3);
t109 = t96 * mrSges(4,1) + t150;
t110 = mrSges(3,1) * t96 - mrSges(3,2) * t95;
t177 = -t110 - mrSges(2,1) - t109;
t174 = qJ(2) * (t71 + t125);
t49 = -qJD(1) * pkin(1) - pkin(2) * t132 - qJ(3) * t133 + qJD(2);
t33 = pkin(3) * t132 - t49;
t172 = m(5) * t33 + mrSges(5,1) * t52 + mrSges(5,2) * t53 + t109 * qJD(1);
t170 = mrSges(2,2) + mrSges(6,2) + mrSges(5,3) - t187;
t144 = -pkin(6) + qJ(2);
t66 = t144 * t96;
t61 = qJD(1) * t66;
t152 = t61 * t98;
t70 = qJ(2) * t133 + qJD(3);
t55 = -pkin(6) * t133 + t70;
t27 = t100 * t55 - t152;
t167 = qJD(5) - t27;
t165 = (pkin(2) + pkin(3)) * t96 + t115;
t28 = t100 * t61 + t55 * t98;
t128 = qJDD(1) * t95;
t56 = t95 * t71 + qJDD(3);
t40 = -pkin(6) * t128 + t56;
t46 = (-pkin(6) * qJDD(1) + t71) * t96;
t5 = -qJD(4) * t28 + t100 * t40 - t46 * t98;
t162 = -t52 / 0.2e1;
t161 = t52 / 0.2e1;
t159 = t53 / 0.2e1;
t154 = Ifges(5,4) * t53;
t141 = t94 * t174;
t101 = cos(qJ(1));
t99 = sin(qJ(1));
t140 = t101 * pkin(1) + t99 * qJ(2);
t137 = t101 * t96;
t102 = qJD(1) ^ 2;
t134 = qJ(2) * t102;
t131 = qJD(3) * t95;
t130 = qJD(4) * t98;
t126 = m(4) - t189;
t123 = t100 * t46 + t55 * t129 + t98 * t40;
t54 = t96 * pkin(3) - t179;
t113 = pkin(2) * t137 + t101 * t87 + t140;
t65 = t144 * t95;
t31 = t100 * t66 + t65 * t98;
t107 = t100 * t65 - t66 * t98;
t106 = -qJD(1) * t131 + qJDD(2);
t29 = qJDD(1) * t165 - t106;
t91 = t101 * qJ(2);
t86 = -qJDD(1) * pkin(1) + qJDD(2);
t82 = t94 * t134;
t79 = mrSges(3,2) * t128;
t51 = -t130 * t96 + t120;
t47 = Ifges(6,5) * t53;
t45 = t58 * t101;
t44 = -t101 * t138 + t137 * t98;
t43 = t58 * t99;
t42 = t59 * t99;
t32 = qJDD(1) * t179 + t106;
t22 = mrSges(6,1) * t52 - mrSges(6,3) * t53;
t21 = pkin(4) * t53 + qJ(5) * t52;
t20 = qJD(4) * qJ(5) + t28;
t19 = -qJD(4) * pkin(4) + t167;
t16 = -Ifges(5,2) * t52 + Ifges(5,6) * qJD(4) + t154;
t15 = Ifges(6,6) * qJD(4) + Ifges(6,3) * t52 + t47;
t10 = pkin(4) * t58 - qJ(5) * t59 + t54;
t7 = pkin(4) * t51 + qJ(5) * t50 - qJD(5) * t59 + t131;
t6 = pkin(4) * t52 - qJ(5) * t53 + t33;
t4 = -t130 * t61 + t123;
t3 = -qJDD(4) * pkin(4) + qJDD(5) - t5;
t2 = qJDD(4) * qJ(5) + (qJD(5) - t152) * qJD(4) + t123;
t1 = pkin(4) * t26 - qJ(5) * t25 - qJD(5) * t53 + t29;
t8 = [(t169 * qJDD(4) + t186 * t25) * t158 + (t184 * t58 + t186 * t59) * t25 / 0.2e1 + (t184 * t51 - t186 * t50) * t159 + (t27 * t50 - t28 * t51 - t4 * t58 - t5 * t59) * mrSges(5,3) + (-t16 / 0.2e1 + t6 * mrSges(6,1) + t15 / 0.2e1 + t33 * mrSges(5,1) + Ifges(6,3) * t161 - Ifges(5,2) * t162) * t51 + (t71 * t93 + t192) * mrSges(3,3) + (t56 * t95 + t192) * mrSges(4,2) - t10 * t25 * mrSges(6,3) + t54 * t25 * mrSges(5,2) - pkin(1) * t79 + (-m(5) * t27 + m(6) * t19 + t142) * (-qJD(2) * t59 + qJD(4) * t31) + (m(5) * t28 + m(6) * t20 + t143) * (qJD(2) * t58 + qJD(4) * t107) - t181 * t50 / 0.2e1 + (m(5) * t4 + m(6) * t2 + t182) * t31 + (m(5) * t5 - m(6) * t3 + t183) * t107 + (t168 * t58 + t169 * t59) * qJDD(4) / 0.2e1 + (t168 * t51 - t169 * t50) * qJD(4) / 0.2e1 + (t185 * t96 + (Ifges(3,1) + Ifges(4,1)) * t95) * t128 + t58 * (Ifges(6,5) * t25 + Ifges(6,6) * qJDD(4)) / 0.2e1 + t1 * (t58 * mrSges(6,1) - t59 * mrSges(6,3)) + m(4) * (t32 * t179 + (qJ(2) * t56 + qJD(2) * t70 - qJD(3) * t49) * t95 + t141) + ((pkin(1) * mrSges(3,1) + (Ifges(3,2) + Ifges(4,3)) * t96 + t185 * t95) * t96 - t109 * t179 + Ifges(2,3)) * qJDD(1) + t7 * t22 + (-m(3) * t140 - m(4) * t113 + t189 * (pkin(3) * t137 - pkin(6) * t99 + t113) - t191 * t45 - t171 * t44 + t177 * t101 + t170 * t99) * g(2) + ((-m(4) - m(3)) * t91 + t189 * (-pkin(6) * t101 - t165 * t99 + t91) + t191 * t43 - t171 * t42 + (m(3) * pkin(1) - m(4) * t179 - t177) * t99 + t170 * t101) * g(1) + (t54 * mrSges(5,1) + t10 * mrSges(6,1) + (Ifges(5,2) + Ifges(6,3)) * t58 + 0.2e1 * t184 * t158) * t26 + (-t33 * mrSges(5,2) + t6 * mrSges(6,3) - Ifges(5,4) * t162 - Ifges(6,5) * t161) * t50 - t32 * t109 - t86 * t110 + (m(5) * t54 + t58 * mrSges(5,1) + t59 * mrSges(5,2)) * t29 + m(3) * (-pkin(1) * t86 + t93 * t174 + t141) + m(6) * (t1 * t10 + t6 * t7) + t172 * t131 - t58 * (Ifges(5,4) * t25 + Ifges(5,6) * qJDD(4)) / 0.2e1 + (-t19 * t50 - t2 * t58 - t20 * t51 + t3 * t59) * mrSges(6,2); t79 + t142 * t53 - t143 * t52 + t146 * t26 + t145 * t25 + (-t150 + (-mrSges(3,1) - mrSges(4,1)) * t96) * qJDD(1) - t187 * t102 * t180 + (-g(1) * t99 + g(2) * t101) * (m(3) + t126) + (t19 * t53 - t20 * t52 - t1) * m(6) + (-t27 * t53 - t28 * t52 - t29) * m(5) + (-t133 * t70 + t32 - t82) * m(4) + (-t134 * t93 - t82 + t86) * m(3); t182 * t98 + t183 * t100 + t126 * t96 * g(3) + (t100 * t143 + t142 * t98) * qJD(4) + m(6) * (-t100 * t3 + t2 * t98 + (t100 * t20 + t19 * t98) * qJD(4)) + m(5) * (t100 * t5 + t4 * t98 + (t100 * t28 - t27 * t98) * qJD(4)) + m(4) * t56 + (qJDD(1) * mrSges(4,2) + (m(4) * t49 - m(6) * t6 - t172 - t22) * qJD(1) + (-g(1) * t101 - g(2) * t99) * t126) * t95; (-pkin(4) * t3 + qJ(5) * t2 + t167 * t20 - t19 * t28 - t21 * t6) * m(6) + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + t169 * t25 - (t168 * t53 - t169 * t52) * qJD(4) / 0.2e1 + (t19 * t52 + t20 * t53) * mrSges(6,2) - (-t186 * t52 + t15 - t154 + t47) * t53 / 0.2e1 + (-t142 + t155) * t28 + (-t143 - t156) * t27 + t168 * t26 + (-Ifges(5,2) * t53 + t181 - t48) * t161 - t33 * (mrSges(5,1) * t53 - mrSges(5,2) * t52) - t6 * (mrSges(6,1) * t53 + mrSges(6,3) * t52) + qJD(5) * t38 - t21 * t22 - pkin(4) * t12 + qJ(5) * t14 + t2 * mrSges(6,3) - t3 * mrSges(6,1) - t4 * mrSges(5,2) + t5 * mrSges(5,1) + (-t171 * t43 - t191 * t42) * g(2) + (-t171 * t45 + t191 * t44) * g(1) + (-t171 * t59 + t191 * t58) * g(3) + t16 * t159 + (Ifges(6,3) * t53 - t153) * t162; -qJD(4) * t38 + t53 * t22 + (-g(1) * t44 + g(2) * t42 - g(3) * t58 - t20 * qJD(4) + t6 * t53 + t3) * m(6) + t12;];
tau = t8;
