% Calculate vector of inverse dynamics joint torques for
% S5PRPRP4
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:04
% EndTime: 2019-12-05 15:35:19
% DurationCPUTime: 4.98s
% Computational Cost: add. (1181->282), mult. (2448->357), div. (0->0), fcn. (1496->10), ass. (0->134)
t84 = sin(qJ(4));
t86 = cos(qJ(4));
t110 = mrSges(5,1) * t86 - mrSges(5,2) * t84;
t213 = -mrSges(4,1) - t110;
t210 = m(5) + m(6);
t212 = m(4) + t210;
t108 = t86 * mrSges(6,1) + t84 * mrSges(6,3);
t211 = -t108 + t213;
t197 = mrSges(5,1) + mrSges(6,1);
t195 = Ifges(6,4) + Ifges(5,5);
t194 = Ifges(5,6) - Ifges(6,6);
t145 = qJD(2) * t86;
t129 = mrSges(5,3) * t145;
t131 = mrSges(6,2) * t145;
t58 = qJD(4) * mrSges(6,3) + t131;
t152 = -qJD(4) * mrSges(5,2) + t129 + t58;
t146 = qJD(2) * t84;
t130 = mrSges(5,3) * t146;
t132 = mrSges(6,2) * t146;
t192 = t197 * qJD(4) - t130 - t132;
t92 = t152 * t86 - t192 * t84;
t209 = -qJD(2) * mrSges(4,2) + t92;
t107 = t84 * mrSges(6,1) - t86 * mrSges(6,3);
t109 = mrSges(5,1) * t84 + mrSges(5,2) * t86;
t87 = cos(qJ(2));
t149 = qJD(1) * t87;
t61 = qJD(2) * pkin(2) + t149;
t85 = sin(qJ(2));
t150 = qJD(1) * t85;
t80 = sin(pkin(8));
t62 = t80 * t150;
t82 = cos(pkin(8));
t20 = t61 * t82 - t62;
t102 = pkin(4) * t86 + qJ(5) * t84;
t97 = -pkin(3) - t102;
t11 = qJD(2) * t97 - t20;
t18 = -qJD(2) * pkin(3) - t20;
t204 = -t11 * t107 - t18 * t109;
t143 = qJD(4) * t86;
t141 = qJDD(2) * pkin(2);
t138 = qJD(1) * qJD(2);
t123 = t85 * t138;
t53 = t87 * qJDD(1) - t123;
t43 = t53 + t141;
t122 = t87 * t138;
t54 = qJDD(1) * t85 + t122;
t13 = t80 * t43 + t82 * t54;
t8 = qJDD(2) * pkin(6) + t13;
t135 = qJD(3) * t143 + t84 * qJDD(3) + t86 * t8;
t144 = qJD(4) * t84;
t21 = t150 * t82 + t80 * t61;
t19 = qJD(2) * pkin(6) + t21;
t3 = -t144 * t19 + t135;
t15 = qJD(3) * t84 + t19 * t86;
t4 = -t15 * qJD(4) + qJDD(3) * t86 - t8 * t84;
t114 = t3 * t86 - t4 * t84;
t166 = t19 * t84;
t14 = qJD(3) * t86 - t166;
t202 = -t14 * t143 - t15 * t144 + t114;
t1 = qJDD(4) * qJ(5) + (qJD(5) - t166) * qJD(4) + t135;
t2 = -qJDD(4) * pkin(4) + qJDD(5) - t4;
t115 = t1 * t86 + t2 * t84;
t10 = qJD(4) * qJ(5) + t15;
t140 = t10 * qJD(4);
t9 = -qJD(4) * pkin(4) + qJD(5) - t14;
t201 = -t84 * t140 + t9 * t143 + t115;
t137 = qJD(2) * qJD(4);
t52 = qJDD(2) * t84 + t137 * t86;
t29 = -qJDD(4) * mrSges(6,1) + t52 * mrSges(6,2);
t154 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t52 - t29;
t51 = -t86 * qJDD(2) + t137 * t84;
t30 = -mrSges(6,2) * t51 + qJDD(4) * mrSges(6,3);
t155 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t51 + t30;
t199 = -t154 * t84 + t155 * t86;
t79 = qJ(2) + pkin(8);
t73 = sin(t79);
t175 = g(3) * t73;
t196 = mrSges(5,2) - mrSges(6,3);
t167 = Ifges(6,5) * t86;
t106 = t84 * Ifges(6,1) - t167;
t72 = Ifges(5,4) * t145;
t193 = Ifges(5,1) * t146 + qJD(2) * t106 + t195 * qJD(4) + t72;
t190 = -t194 * t84 + t195 * t86;
t47 = t108 * qJD(2);
t186 = t213 * qJD(2) - t47;
t184 = qJD(2) ^ 2;
t182 = t84 / 0.2e1;
t181 = pkin(2) * t80;
t180 = pkin(2) * t82;
t78 = t87 * pkin(2);
t170 = Ifges(5,4) * t84;
t169 = Ifges(5,4) * t86;
t168 = Ifges(6,5) * t84;
t45 = t80 * t85 - t82 * t87;
t35 = t45 * qJD(2);
t163 = t35 * t84;
t162 = t35 * t86;
t81 = sin(pkin(7));
t159 = t81 * t84;
t158 = t81 * t86;
t83 = cos(pkin(7));
t157 = t83 * t84;
t156 = t83 * t86;
t142 = qJD(5) * t84;
t139 = qJDD(2) * mrSges(4,2);
t12 = t43 * t82 - t80 * t54;
t118 = -t137 / 0.2e1;
t116 = -g(1) * t81 + g(2) * t83;
t113 = t10 * t86 + t84 * t9;
t112 = mrSges(3,1) * t87 - mrSges(3,2) * t85;
t111 = mrSges(3,1) * t85 + mrSges(3,2) * t87;
t105 = t86 * Ifges(5,2) + t170;
t101 = pkin(4) * t84 - qJ(5) * t86;
t100 = -t14 * t84 + t15 * t86;
t99 = t80 * t87 + t82 * t85;
t7 = -qJDD(2) * pkin(3) - t12;
t94 = t84 * (Ifges(5,1) * t86 - t170);
t93 = t86 * (Ifges(6,3) * t84 + t167);
t74 = cos(t79);
t71 = Ifges(6,5) * t146;
t49 = t101 * qJD(2);
t42 = t97 - t180;
t37 = Ifges(5,6) * qJD(4) + qJD(2) * t105;
t36 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t145 + t71;
t33 = t99 * qJD(2);
t31 = qJD(4) * t101 - t142;
t26 = t156 * t74 + t159;
t25 = t157 * t74 - t158;
t24 = t158 * t74 - t157;
t23 = t159 * t74 + t156;
t17 = mrSges(5,1) * t51 + mrSges(5,2) * t52;
t16 = mrSges(6,1) * t51 - mrSges(6,3) * t52;
t5 = pkin(4) * t51 - qJ(5) * t52 - qJD(2) * t142 + t7;
t6 = [m(2) * qJDD(1) - t111 * t184 + (t16 + t17) * t45 + t186 * t33 + (-mrSges(4,1) * t45 + t112) * qJDD(2) - t209 * t35 + (-m(2) - m(3) - t212) * g(3) + m(5) * (t14 * t163 - t15 * t162 + t18 * t33 + t45 * t7) + m(6) * (-t10 * t162 + t11 * t33 - t163 * t9 + t45 * t5) + m(4) * (-t12 * t45 - t20 * t33 - t21 * t35) + m(3) * (t53 * t87 + t54 * t85) + (-t139 + (-t152 * t84 - t192 * t86) * qJD(4) + m(5) * t202 + m(6) * t201 + m(4) * t13 + t199) * t99; (-m(4) * t21 - m(5) * t100 - m(6) * t113 - t209) * (t149 * t82 - t62) + (-t37 / 0.2e1 + t36 / 0.2e1) * t144 + (t53 + t123) * mrSges(3,1) + (-t105 / 0.2e1 + t168 / 0.2e1 + (-Ifges(5,2) / 0.2e1 - Ifges(6,3)) * t86 + (Ifges(6,5) - Ifges(5,4)) * t182) * t51 + (t84 * Ifges(5,1) + t106 + t169) * t52 / 0.2e1 + (t84 * (Ifges(6,1) * t86 + t168) + t86 * (-Ifges(5,2) * t84 + t169) + t94) * t137 / 0.2e1 + (m(5) * t7 + t17) * (-pkin(3) - t180) + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) - t86 * (Ifges(6,5) * t52 + Ifges(6,6) * qJDD(4)) / 0.2e1 + (m(4) * t20 - m(5) * t18 - m(6) * t11 - t186) * t99 * qJD(1) + t86 * (Ifges(5,4) * t52 + Ifges(5,6) * qJDD(4)) / 0.2e1 + ((Ifges(5,1) + Ifges(6,1)) * t52 + t195 * qJDD(4)) * t182 + (t194 * t86 + t195 * t84) * qJDD(4) / 0.2e1 + (-t175 + t201) * mrSges(6,2) + (-t175 + t202) * mrSges(5,3) + (-t152 * t144 - t192 * t143 + m(5) * ((-t14 * t86 - t15 * t84) * qJD(4) + t114) + m(6) * ((-t10 * t84 + t86 * t9) * qJD(4) + t115) + t199) * (pkin(6) + t181) + (t190 * qJD(4) / 0.2e1 - t204) * qJD(4) - t31 * t47 + (-t54 + t122) * mrSges(3,2) + m(6) * (t11 * t31 + t42 * t5) + (t111 + (-pkin(6) * t210 + mrSges(4,2) - mrSges(6,2) - mrSges(5,3)) * t74 + (m(5) * pkin(3) - m(6) * t97 - t211) * t73 + t212 * pkin(2) * t85) * (g(1) * t83 + g(2) * t81) + m(4) * (t12 * t82 + t13 * t80) * pkin(2) + t42 * t16 - t13 * mrSges(4,2) + t193 * t143 / 0.2e1 + t93 * t118 - t5 * t108 - t7 * t110 + (t141 * t82 + t12) * mrSges(4,1) + (-m(4) * t78 + mrSges(4,2) * t73 - t112 - t210 * (t74 * pkin(3) + t73 * pkin(6) + t78) + (-m(6) * t102 + t211) * t74) * g(3) - t139 * t181; t154 * t86 + t155 * t84 + t92 * qJD(4) + (qJD(4) * t113 + t1 * t84 - t2 * t86 + t116) * m(6) + (qJD(4) * t100 + t3 * t84 + t4 * t86 + t116) * m(5) + (qJDD(3) + t116) * m(4); (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t10 - t11 * t49 - g(1) * (-pkin(4) * t25 + qJ(5) * t26) - g(2) * (-pkin(4) * t23 + qJ(5) * t24) + t101 * t175) * m(6) + t190 * t118 + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) - (Ifges(6,1) * t145 + t36 + t71) * t146 / 0.2e1 - t194 * t51 + t195 * t52 + (t196 * t26 + t197 * t25) * g(1) + (t196 * t24 + t197 * t23) * g(2) + t204 * qJD(2) + (t93 / 0.2e1 - t94 / 0.2e1) * t184 + (t107 + t109) * t175 + qJD(5) * t58 + t49 * t47 - pkin(4) * t29 + qJ(5) * t30 + t4 * mrSges(5,1) + t1 * mrSges(6,3) - t2 * mrSges(6,1) - t3 * mrSges(5,2) + (-m(6) * t10 + t129 - t152) * t14 + (-m(6) * t9 + t130 + t192) * t15 - (-Ifges(5,2) * t146 + t193 + t72) * t145 / 0.2e1 + t10 * t132 - t9 * t131 + t37 * t146 / 0.2e1; -t47 * t146 - qJD(4) * t58 + (-g(1) * t25 - g(2) * t23 + t11 * t146 - t175 * t84 - t140 + t2) * m(6) + t29;];
tau = t6;
