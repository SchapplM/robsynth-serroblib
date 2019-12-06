% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:17
% EndTime: 2019-12-05 17:59:21
% DurationCPUTime: 1.14s
% Computational Cost: add. (3320->138), mult. (5645->173), div. (0->0), fcn. (5462->4), ass. (0->83)
t102 = sin(qJ(4));
t103 = sin(qJ(3));
t104 = cos(qJ(4));
t105 = cos(qJ(3));
t93 = t102 * t103 - t104 * t105;
t148 = t93 * mrSges(6,2);
t94 = -t102 * t105 - t104 * t103;
t160 = t94 * pkin(4);
t97 = t103 * pkin(3) + qJ(2);
t70 = t97 - t160;
t163 = m(6) * t70;
t82 = t94 * mrSges(6,1);
t186 = -t82 - t148 + t163;
t106 = -pkin(1) - pkin(6);
t129 = t103 * t106;
t95 = -t103 * pkin(7) + t129;
t96 = (-pkin(7) + t106) * t105;
t137 = -t102 * t95 + t104 * t96;
t170 = t93 * qJ(5) + t137;
t61 = t102 * t96 + t104 * t95;
t47 = t94 * qJ(5) + t61;
t185 = (Ifges(5,5) + Ifges(6,5)) * t94 + (Ifges(5,6) + Ifges(6,6)) * t93 - t137 * mrSges(5,2) - t170 * mrSges(6,2) - t61 * mrSges(5,1) - t47 * mrSges(6,1);
t158 = t104 * pkin(3);
t99 = pkin(4) + t158;
t147 = t99 * t47;
t182 = m(5) * t97;
t166 = m(6) / 0.2e1;
t179 = Ifges(6,4) + Ifges(5,4);
t177 = t99 - t158;
t175 = mrSges(6,3) * t170 + t179 * t94;
t173 = pkin(3) * t102;
t144 = mrSges(5,2) + mrSges(6,2);
t171 = t144 * t93;
t169 = t93 ^ 2;
t168 = t94 ^ 2;
t161 = pkin(4) * t47;
t159 = mrSges(6,3) * t94;
t157 = t105 * pkin(3);
t146 = t99 * t93;
t145 = t99 * t94;
t83 = t94 * mrSges(5,1);
t139 = t82 + t83;
t136 = t102 * t93;
t134 = t104 * t94;
t133 = t105 * mrSges(4,2);
t121 = t168 + t169;
t14 = m(6) * (t170 * t93 + t47 * t94) + t121 * mrSges(6,3);
t128 = t14 * qJD(1);
t111 = t103 * mrSges(4,1) + t133 - t139;
t19 = mrSges(3,3) - t171 + (m(4) + m(3)) * qJ(2) + t182 + t163 + t111;
t127 = t19 * qJD(1);
t123 = t94 * t173;
t81 = t94 * mrSges(6,2);
t56 = -t93 * mrSges(6,1) + t81;
t76 = -t93 * pkin(4) + t157;
t21 = 0.2e1 * (t123 / 0.4e1 + t146 / 0.4e1 - t76 / 0.4e1) * m(6) - t56;
t126 = t21 * qJD(1);
t25 = (-t169 / 0.2e1 - t168 / 0.2e1 - 0.1e1 / 0.2e1) * m(6);
t125 = t25 * qJD(1);
t120 = m(6) * pkin(4) + mrSges(6,1);
t34 = t120 * t93 - t81;
t124 = t34 * qJD(1);
t119 = t179 * t93;
t118 = -Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2);
t116 = -t99 / 0.2e1 + t158 / 0.2e1;
t112 = -t170 * t159 + t70 * t56 + t97 * (-t93 * mrSges(5,1) + t94 * mrSges(5,2));
t1 = (qJ(2) * mrSges(4,1) - Ifges(4,4) * t105 - pkin(3) * t83) * t105 + t157 * t182 + t175 * t94 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t103 + (-Ifges(4,1) + Ifges(4,2)) * t105) * t103 + (-mrSges(5,2) * t157 + t118 * t94 - t119) * t93 + t112 + t186 * t76;
t115 = t1 * qJD(1);
t2 = (-t186 * pkin(4) - t119) * t93 + (t118 * t93 + t175) * t94 + t112;
t114 = t2 * qJD(1);
t113 = t139 + t171;
t110 = (pkin(4) / 0.2e1 + t116) * t94;
t109 = t47 * t158;
t20 = m(6) * t110;
t4 = (-t147 / 0.2e1 + t109 / 0.2e1 + t161 / 0.2e1) * m(6) + mrSges(6,3) * t110;
t43 = t144 * t158 + (t177 * m(6) + mrSges(5,1) + mrSges(6,1)) * t173;
t108 = t4 * qJD(1) - t20 * qJD(2) - t43 * qJD(3);
t77 = pkin(3) * t136;
t24 = -m(6) * t121 / 0.2e1 + t166;
t23 = (t123 + t146 + t76) * t166;
t15 = t113 + (pkin(4) + t177) * t94 * t166;
t3 = (t109 - t147 - t161) * t166 + (-pkin(4) / 0.2e1 + t116) * t159 + t185;
t5 = [t19 * qJD(2) + t1 * qJD(3) + t2 * qJD(4) + t14 * qJD(5), t24 * qJD(5) + t127, t3 * qJD(4) + t23 * qJD(5) + t115 + (-m(6) * t147 - mrSges(4,1) * t129 - mrSges(6,3) * t145 - Ifges(4,5) * t103 - Ifges(4,6) * t105 - t106 * t133 + (mrSges(6,3) * t136 + (-t134 + t136) * mrSges(5,3) + m(6) * t102 * t170 + m(5) * (t102 * t137 - t104 * t61)) * pkin(3) + t185) * qJD(3), t3 * qJD(3) + ((-m(6) * t47 - t159) * pkin(4) + t185) * qJD(4) + t114, t24 * qJD(2) + t23 * qJD(3) + t128; t25 * qJD(5) - t127, 0, t15 * qJD(4) + (t93 * mrSges(5,2) - t111 + t148 + m(5) * (pkin(3) * t134 - t77) + 0.2e1 * (-t77 + t145) * t166) * qJD(3), t15 * qJD(3) + (m(6) * t160 + t113) * qJD(4), t125; t4 * qJD(4) + t21 * qJD(5) - t115, -t20 * qJD(4), -t43 * qJD(4), (-t144 * t104 + (-mrSges(5,1) - t120) * t102) * qJD(4) * pkin(3) + t108, t126; -t4 * qJD(3) + t34 * qJD(5) - t114, t20 * qJD(3), -t108, 0, t124; -t25 * qJD(2) - t21 * qJD(3) - t34 * qJD(4) - t128, -t125, -t126, -t124, 0;];
Cq = t5;
