% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:59
% EndTime: 2019-12-31 17:06:02
% DurationCPUTime: 1.32s
% Computational Cost: add. (2624->182), mult. (5399->271), div. (0->0), fcn. (5503->6), ass. (0->109)
t161 = sin(qJ(2));
t130 = t161 * pkin(2);
t181 = m(4) * t130;
t133 = sin(pkin(7));
t134 = cos(pkin(7));
t162 = cos(qJ(2));
t83 = t133 * t161 - t134 * t162;
t84 = -t133 * t162 - t134 * t161;
t180 = -t84 * mrSges(4,1) - t83 * mrSges(4,2);
t109 = cos(qJ(4));
t102 = Ifges(5,5) * t109;
t108 = sin(qJ(4));
t149 = Ifges(5,6) * t108;
t179 = Ifges(4,4) - t102 / 0.2e1 + t149 / 0.2e1;
t103 = Ifges(5,4) * t109;
t153 = Ifges(5,1) * t108;
t92 = t103 + t153;
t106 = t108 ^ 2;
t107 = t109 ^ 2;
t178 = t106 + t107;
t56 = -t84 * pkin(3) + t83 * pkin(6) + t130;
t105 = t162 * pkin(5);
t132 = t162 * qJ(3) + t105;
t129 = t161 * pkin(5);
t88 = -t161 * qJ(3) - t129;
t66 = t133 * t132 - t134 * t88;
t29 = t108 * t66 + t109 * t56;
t30 = t108 * t56 - t109 * t66;
t177 = -t29 * t108 + t30 * t109;
t118 = Ifges(5,2) * t108 - t103;
t176 = t134 * t132 + t133 * t88;
t175 = t84 ^ 2;
t174 = m(5) / 0.2e1;
t173 = m(4) * pkin(2);
t172 = -mrSges(5,1) / 0.2e1;
t171 = mrSges(5,2) / 0.2e1;
t170 = -t83 / 0.2e1;
t169 = t83 / 0.4e1;
t168 = -t84 / 0.2e1;
t122 = t133 * pkin(2);
t98 = t122 + pkin(6);
t167 = -t98 / 0.2e1;
t166 = -t108 / 0.2e1;
t165 = t108 / 0.2e1;
t164 = -t109 / 0.2e1;
t163 = t109 / 0.2e1;
t160 = mrSges(5,3) * t84;
t159 = Ifges(5,5) * t83;
t158 = Ifges(5,6) * t83;
t156 = t66 * t84;
t155 = t83 * mrSges(5,1);
t154 = t83 * mrSges(5,2);
t152 = Ifges(5,4) * t108;
t150 = Ifges(5,2) * t109;
t136 = t83 * t108;
t142 = t109 * t83;
t100 = -t162 * pkin(2) - pkin(1);
t55 = t83 * pkin(3) + t84 * pkin(6) + t100;
t26 = t108 * t55 + t109 * t176;
t143 = t109 * t26;
t25 = -t108 * t176 + t109 * t55;
t131 = t108 * t160;
t58 = t131 - t154;
t60 = t109 * t160 + t155;
t144 = t109 * mrSges(5,2);
t146 = t108 * mrSges(5,1);
t89 = t144 + t146;
t5 = -t60 * t136 + t58 * t142 - t175 * t89 - m(5) * (-t156 + (t108 * t25 - t143) * t83) - m(4) * (-t176 * t83 - t156) + (-t83 ^ 2 - t175) * mrSges(4,3);
t148 = qJD(1) * t5;
t38 = t118 * t84 + t158;
t145 = t108 * t38;
t37 = -Ifges(5,6) * t84 + t118 * t83;
t93 = Ifges(5,1) * t109 - t152;
t39 = -Ifges(5,5) * t84 - t93 * t83;
t40 = -t84 * t93 + t159;
t54 = t89 * t83;
t57 = t84 * mrSges(5,2) + mrSges(5,3) * t136;
t59 = -t84 * mrSges(5,1) + mrSges(5,3) * t142;
t1 = -t66 * t54 + t26 * t57 + t30 * t58 + t25 * t59 + t29 * t60 - pkin(1) * (t161 * mrSges(3,1) + t162 * mrSges(3,2)) + m(5) * (t66 * t176 + t25 * t29 + t26 * t30) + (mrSges(4,1) * t130 + t145 / 0.2e1 + t40 * t164 + t179 * t83) * t83 + (-t176 * t89 - mrSges(4,2) * t130 + t37 * t165 + t39 * t164 - t179 * t84 + (-Ifges(5,3) - Ifges(4,2) + Ifges(4,1)) * t83) * t84 + (-Ifges(3,2) + Ifges(3,1)) * t162 * t161 + (-t161 ^ 2 + t162 ^ 2) * Ifges(3,4) + (t180 + t181) * t100;
t147 = t1 * qJD(1);
t141 = t109 * t98;
t117 = Ifges(5,5) * t108 + Ifges(5,6) * t109;
t87 = -t109 * mrSges(5,1) + t108 * mrSges(5,2);
t90 = t150 + t152;
t4 = -t87 * t156 + t26 * t60 + (t40 * t166 + t38 * t164 - mrSges(5,3) * t143 + t117 * t170 + (t92 * t163 + t90 * t166) * t84) * t84 + (-t58 + t131) * t25;
t138 = t4 * qJD(1);
t123 = t134 * pkin(2);
t99 = -t123 - pkin(3);
t110 = (-t178 * t98 * t83 - t99 * t84) * t174 + t87 * t168 + (-t133 * t83 + t134 * t84) * t173 / 0.2e1 + t178 * mrSges(5,3) * t170;
t112 = (t108 * t30 + t109 * t29) * t174 + t57 * t165 + t59 * t163 + t181 / 0.2e1;
t6 = t110 - t112 - t180;
t137 = t6 * qJD(1);
t9 = (t154 / 0.2e1 - t58 / 0.2e1) * t109 + (t155 / 0.2e1 + t60 / 0.2e1) * t108;
t135 = t9 * qJD(1);
t128 = t66 * t89 / 0.2e1;
t126 = t58 * t167;
t125 = -t142 / 0.2e1;
t124 = t136 / 0.2e1;
t119 = t60 * t167 + t40 / 0.4e1;
t24 = t99 * t89 + (t92 / 0.2e1 - t118 / 0.2e1) * t109 + (t93 / 0.2e1 - t90 / 0.2e1) * t108;
t111 = (t107 / 0.2e1 + t106 / 0.2e1) * t98 * mrSges(5,3) + (-t93 / 0.4e1 + t90 / 0.4e1 + t99 * t172 + t150 / 0.4e1) * t109 + (t92 / 0.4e1 - t118 / 0.4e1 + t103 / 0.2e1 + t99 * t171 + t153 / 0.4e1) * t108;
t115 = t30 * t171 + t29 * t172;
t3 = t128 + t102 * t169 + (t159 / 0.2e1 + t119) * t109 + (t126 - 0.3e1 / 0.4e1 * t158 - t38 / 0.4e1) * t108 + (Ifges(5,3) / 0.2e1 + t111) * t84 + t115;
t116 = -t3 * qJD(1) - t24 * qJD(2);
t113 = t84 * t117;
t10 = t58 * t163 + t60 * t166 + (t144 / 0.2e1 + t146 / 0.2e1) * t83;
t7 = t110 + t112;
t2 = t108 * t126 + t128 + (t102 - t149) * t169 - t145 / 0.4e1 + Ifges(5,5) * t125 + Ifges(5,6) * t124 + Ifges(5,3) * t168 + t119 * t109 + t111 * t84 - t115;
t8 = [qJD(2) * t1 - qJD(3) * t5 - qJD(4) * t4, t147 + (-Ifges(4,5) * t83 + Ifges(4,6) * t84 + Ifges(3,5) * t162 - Ifges(3,6) * t161 + t92 * t125 + t90 * t124 + mrSges(3,2) * t129 - mrSges(3,1) * t105 + t57 * t141 - t113 / 0.2e1 + t39 * t165 + t37 * t163 - t99 * t54 + (m(5) * t177 - t108 * t59) * t98 - (t133 * t173 - mrSges(4,2)) * t66 + (m(5) * t99 - t134 * t173 - mrSges(4,1) + t87) * t176 + t177 * mrSges(5,3) + (t84 * t122 + t83 * t123) * mrSges(4,3)) * qJD(2) + t7 * qJD(3) + t2 * qJD(4), qJD(2) * t7 + qJD(4) * t10 - t148, -t138 + t2 * qJD(2) + t10 * qJD(3) + (-t26 * mrSges(5,1) - t25 * mrSges(5,2) + t113) * qJD(4); qJD(3) * t6 + qJD(4) * t3 - t147, t24 * qJD(4), t137, (-mrSges(5,1) * t141 + t102 + (mrSges(5,2) * t98 - Ifges(5,6)) * t108) * qJD(4) - t116; -qJD(2) * t6 - qJD(4) * t9 + t148, -t137, 0, -t89 * qJD(4) - t135; -qJD(2) * t3 + qJD(3) * t9 + t138, t116, t135, 0;];
Cq = t8;
