% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:06
% EndTime: 2019-12-31 16:36:08
% DurationCPUTime: 0.93s
% Computational Cost: add. (1686->197), mult. (4617->323), div. (0->0), fcn. (4218->8), ass. (0->116)
t114 = cos(qJ(4));
t103 = Ifges(5,5) * t114;
t111 = sin(qJ(4));
t166 = Ifges(5,6) * t111;
t186 = Ifges(4,4) - t103 / 0.2e1 + t166 / 0.2e1;
t104 = Ifges(5,4) * t114;
t185 = -Ifges(5,2) * t111 + t104;
t93 = Ifges(5,1) * t111 + t104;
t112 = sin(qJ(3));
t150 = t111 * t112;
t115 = cos(qJ(3));
t97 = pkin(3) * t112 - pkin(7) * t115;
t55 = pkin(6) * t150 + t114 * t97;
t148 = t112 * t114;
t56 = -pkin(6) * t148 + t111 * t97;
t184 = -t55 * t111 + t56 * t114;
t136 = mrSges(5,1) * t114 - mrSges(5,2) * t111;
t183 = -m(5) * pkin(3) - mrSges(4,1) - t136;
t182 = -m(5) / 0.2e1;
t181 = m(5) / 0.2e1;
t110 = sin(pkin(4));
t116 = cos(qJ(2));
t151 = t110 * t116;
t113 = sin(qJ(2));
t146 = t113 * t110;
t152 = cos(pkin(4));
t73 = t152 * t112 + t115 * t146;
t35 = -t111 * t151 + t114 * t73;
t180 = -t35 / 0.2e1;
t75 = t136 * t112;
t178 = -t75 / 0.2e1;
t141 = mrSges(5,3) * t150;
t81 = mrSges(5,2) * t115 - t141;
t177 = -t81 / 0.2e1;
t160 = t114 * mrSges(5,2);
t164 = t111 * mrSges(5,1);
t90 = t160 + t164;
t176 = t90 / 0.2e1;
t175 = -t111 / 0.2e1;
t174 = t112 / 0.2e1;
t173 = t114 / 0.2e1;
t169 = Ifges(5,4) * t111;
t168 = Ifges(5,5) * t115;
t165 = Ifges(5,6) * t115;
t140 = mrSges(5,3) * t148;
t83 = -mrSges(5,1) * t115 - t140;
t163 = t111 * t83;
t162 = t112 * mrSges(4,2);
t72 = t112 * t146 - t152 * t115;
t161 = t112 * t72;
t159 = t114 * t81;
t34 = -t73 * t111 - t114 * t151;
t158 = t34 * t111;
t157 = t35 * t114;
t144 = t115 * t116;
t51 = (-t111 * t144 + t113 * t114) * t110;
t156 = t51 * t111;
t52 = (t111 * t113 + t114 * t144) * t110;
t155 = t52 * t114;
t149 = t111 * t115;
t147 = t112 * t116;
t145 = t114 * t115;
t13 = m(5) * (-t157 + t73 + t158) * t72;
t143 = t13 * qJD(1);
t14 = 0.4e1 * (t115 * t73 / 0.4e1 + t161 / 0.4e1 - t146 / 0.4e1) * t151 * m(4) + (t161 * t151 + t34 * t51 + t35 * t52) * m(5);
t142 = t14 * qJD(1);
t139 = t110 * t147;
t94 = Ifges(5,1) * t114 - t169;
t91 = Ifges(5,2) * t114 + t169;
t134 = Ifges(5,5) * t111 + Ifges(5,6) * t114;
t117 = (t155 / 0.2e1 - t156 / 0.2e1) * mrSges(5,3) + (-pkin(3) * t139 + (t155 - t156) * pkin(7)) * t181 - t136 * t139 / 0.2e1;
t76 = t90 * t112;
t82 = -mrSges(5,2) * t112 - mrSges(5,3) * t149;
t84 = mrSges(5,1) * t112 - mrSges(5,3) * t145;
t119 = -t34 * t84 / 0.2e1 + t82 * t180 - t73 * t76 / 0.2e1;
t121 = pkin(6) * t112 * t73 + t34 * t55 + t56 * t35;
t87 = -pkin(3) * t115 - pkin(7) * t112 - pkin(2);
t53 = -pkin(6) * t149 + t114 * t87;
t54 = pkin(6) * t145 + t111 * t87;
t131 = t111 * t53 - t114 * t54;
t123 = pkin(6) * t115 + t131;
t77 = t90 * t115;
t1 = t121 * t182 + (-t77 / 0.2e1 + t123 * t182 + t159 / 0.2e1 - t163 / 0.2e1) * t72 + t117 + t119;
t60 = t112 * t185 - t165;
t61 = Ifges(5,6) * t112 + t115 * t185;
t62 = t112 * t94 - t168;
t63 = Ifges(5,5) * t112 + t94 * t115;
t3 = m(5) * (t53 * t55 + t54 * t56) + t56 * t81 + t54 * t82 + t55 * t83 + t53 * t84 + (-pkin(2) * mrSges(4,1) + pkin(6) * t77 - t186 * t112 + t63 * t173 + t61 * t175) * t112 + (-pkin(2) * mrSges(4,2) + pkin(6) * t76 + t62 * t173 + t60 * t175 + (m(5) * pkin(6) ^ 2 + Ifges(4,1) - Ifges(4,2) - Ifges(5,3)) * t112 + t186 * t115) * t115;
t133 = -t1 * qJD(1) + t3 * qJD(2);
t120 = t34 * t177 + t35 * t83 / 0.2e1 + t72 * t178;
t128 = t51 * mrSges(5,1) / 0.2e1 - t52 * mrSges(5,2) / 0.2e1;
t4 = (t157 / 0.2e1 - t158 / 0.2e1) * t112 * mrSges(5,3) + t120 + t128;
t78 = t112 * t91;
t79 = t112 * t93;
t8 = t53 * t81 - t54 * t83 + (t115 * t134 / 0.2e1 + pkin(6) * t75 - t79 * t173 - t114 * t60 / 0.2e1 + t131 * mrSges(5,3) + (t62 - t78) * t175) * t112;
t132 = -t4 * qJD(1) + t8 * qJD(2);
t130 = -t60 / 0.4e1 - t79 / 0.4e1 + pkin(7) * t177;
t129 = -pkin(7) * t83 / 0.2e1 + t62 / 0.4e1 - t78 / 0.4e1;
t127 = -t55 * mrSges(5,1) / 0.2e1 + t56 * mrSges(5,2) / 0.2e1;
t126 = pkin(3) * t178 - t115 * t103 / 0.4e1;
t124 = t160 / 0.2e1 + t164 / 0.2e1;
t10 = (-t90 / 0.2e1 + t124) * t72;
t18 = -pkin(3) * t90 + (t93 / 0.2e1 + t185 / 0.2e1) * t114 + (t94 / 0.2e1 - t91 / 0.2e1) * t111;
t106 = t111 ^ 2;
t108 = t114 ^ 2;
t118 = pkin(6) * t176 + (t94 / 0.4e1 - t91 / 0.4e1) * t114 + (-t93 / 0.4e1 - t185 / 0.4e1) * t111 + (-t108 / 0.2e1 - t106 / 0.2e1) * pkin(7) * mrSges(5,3);
t7 = (-t168 / 0.2e1 + t129) * t114 + (0.3e1 / 0.4e1 * t165 + t130) * t111 + (-Ifges(5,3) / 0.2e1 + t118) * t112 + t126 + t127;
t122 = t10 * qJD(1) - t7 * qJD(2) - t18 * qJD(3);
t109 = t115 ^ 2;
t107 = t112 ^ 2;
t88 = t107 * pkin(6) * t151;
t11 = (t124 + t176) * t72;
t6 = Ifges(5,5) * t145 / 0.2e1 - Ifges(5,6) * t149 / 0.2e1 + Ifges(5,3) * t174 + (t165 / 0.4e1 + t130) * t111 + t129 * t114 + t118 * t112 + t126 - t127;
t5 = t140 * t180 + t34 * t141 / 0.2e1 - t120 + t128;
t2 = (t123 * t72 + t121) * t181 - t72 * t159 / 0.2e1 + t117 - t119 + (t77 + t163) * t72 / 0.2e1 - (mrSges(4,1) * t112 + mrSges(4,2) * t115) * t151;
t9 = [t14 * qJD(2) + t13 * qJD(3), t2 * qJD(3) + t5 * qJD(4) + t142 + (t51 * t83 + t52 * t81 + m(4) * t88 + 0.2e1 * (t51 * t53 + t52 * t54 + t88) * t181 + (t76 * t147 + (-mrSges(3,2) + m(4) * pkin(6) * t109 + (t107 + t109) * mrSges(4,3)) * t116 + (-m(4) * pkin(2) - t115 * mrSges(4,1) - mrSges(3,1) + t162) * t113) * t110) * qJD(2), t143 + t2 * qJD(2) + (t183 * t73 + (mrSges(4,2) + (m(5) * pkin(7) + mrSges(5,3)) * (-t106 - t108)) * t72) * qJD(3) + t11 * qJD(4), t5 * qJD(2) + t11 * qJD(3) + (-mrSges(5,1) * t35 - mrSges(5,2) * t34) * qJD(4); -qJD(3) * t1 - qJD(4) * t4 - t142, qJD(3) * t3 + qJD(4) * t8, t6 * qJD(4) + t133 + (pkin(6) * t162 - pkin(3) * t77 + t134 * t174 + t111 * t63 / 0.2e1 + t61 * t173 - Ifges(4,6) * t112 + (m(5) * t184 - t111 * t84 + t114 * t82) * pkin(7) + (t183 * pkin(6) + t93 * t173 + t91 * t175 + Ifges(4,5)) * t115 + t184 * mrSges(5,3)) * qJD(3), t6 * qJD(3) + (-mrSges(5,1) * t54 - mrSges(5,2) * t53 - t112 * t134) * qJD(4) + t132; qJD(2) * t1 - qJD(4) * t10 - t143, qJD(4) * t7 - t133, t18 * qJD(4), (-t136 * pkin(7) + t103 - t166) * qJD(4) - t122; t4 * qJD(2) + t10 * qJD(3), -qJD(3) * t7 - t132, t122, 0;];
Cq = t9;
