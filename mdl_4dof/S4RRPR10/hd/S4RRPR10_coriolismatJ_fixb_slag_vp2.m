% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:11
% EndTime: 2019-12-31 17:11:13
% DurationCPUTime: 0.81s
% Computational Cost: add. (1346->180), mult. (2743->248), div. (0->0), fcn. (2102->4), ass. (0->101)
t85 = cos(qJ(2));
t159 = -0.2e1 * t85;
t158 = -Ifges(3,4) - Ifges(4,6);
t83 = sin(qJ(2));
t104 = t83 * pkin(2) - qJ(3) * t85;
t56 = pkin(6) * t83 + t104;
t144 = pkin(3) + pkin(5);
t70 = t144 * t85;
t82 = sin(qJ(4));
t84 = cos(qJ(4));
t25 = -t56 * t82 + t84 * t70;
t121 = t84 * t25;
t26 = t56 * t84 + t82 * t70;
t129 = t82 * t26;
t157 = -t129 - t121;
t141 = mrSges(5,3) * t85;
t154 = t84 ^ 2;
t156 = (t82 ^ 2 + t154) * t141 / 0.2e1;
t152 = -Ifges(5,2) / 0.2e1;
t151 = Ifges(5,3) / 0.2e1;
t150 = -t82 / 0.2e1;
t148 = t82 / 0.2e1;
t147 = -t84 / 0.2e1;
t145 = t84 / 0.2e1;
t86 = -pkin(2) - pkin(6);
t143 = qJ(3) / 0.2e1;
t142 = mrSges(5,3) * t83;
t140 = Ifges(5,1) * t82;
t139 = Ifges(5,4) * t82;
t138 = Ifges(5,4) * t84;
t137 = Ifges(5,5) * t82;
t136 = Ifges(5,5) * t83;
t135 = Ifges(5,5) * t84;
t134 = Ifges(5,6) * t82;
t133 = Ifges(5,6) * t83;
t132 = Ifges(5,6) * t84;
t131 = t82 * mrSges(5,1);
t130 = t82 * mrSges(5,2);
t127 = t82 * t85;
t110 = mrSges(5,3) * t127;
t125 = t83 * mrSges(5,1);
t58 = t110 + t125;
t128 = t82 * t58;
t126 = t82 * t86;
t124 = t83 * mrSges(5,2);
t123 = t84 * mrSges(5,1);
t122 = t84 * mrSges(5,2);
t119 = t84 * t85;
t109 = mrSges(5,3) * t119;
t60 = -t109 - t124;
t120 = t84 * t60;
t118 = t84 * t86;
t117 = t85 * mrSges(5,1);
t116 = t85 * mrSges(5,2);
t115 = qJ(3) * t83;
t98 = -pkin(2) * t85 - t115;
t63 = -pkin(1) + t98;
t64 = t85 * mrSges(4,2) - t83 * mrSges(4,3);
t105 = m(4) * t63 + t64;
t53 = t85 * t86 - pkin(1) - t115;
t69 = t144 * t83;
t23 = -t53 * t82 + t69 * t84;
t24 = t53 * t84 + t69 * t82;
t6 = (m(5) * (t23 * t82 - t24 * t84) - t120 + t128 - t105) * t83;
t114 = qJD(1) * t6;
t102 = t123 - t130;
t100 = Ifges(5,2) * t84 + t139;
t40 = Ifges(5,6) * t85 + t100 * t83;
t101 = t138 + t140;
t41 = Ifges(5,5) * t85 + t101 * t83;
t52 = t102 * t83;
t57 = -t142 * t82 + t117;
t59 = t142 * t84 - t116;
t94 = t137 / 0.2e1 + t132 / 0.2e1;
t99 = t132 + t137;
t1 = -t70 * t52 + t23 * t57 + t25 * t58 + t24 * t59 + t26 * t60 + m(5) * (t23 * t25 + t24 * t26 - t69 * t70) + t105 * t104 + (-pkin(1) * mrSges(3,2) - t63 * mrSges(4,3) - t69 * t102 + t40 * t147 + t41 * t150 + (-t94 - t158) * t85) * t85 + (-t63 * mrSges(4,2) - pkin(1) * mrSges(3,1) + (t99 + t158) * t83 + (t154 * t152 + Ifges(5,3) + Ifges(3,1) - Ifges(3,2) - Ifges(4,3) + Ifges(4,2) + (-t138 - t140 / 0.2e1) * t82) * t85) * t83;
t113 = t1 * qJD(1);
t65 = t122 + t131;
t89 = t65 * t85;
t4 = t70 * t89 + ((-Ifges(5,4) * t119 + t136) * t84 + (Ifges(5,4) * t127 - t133 + (-Ifges(5,1) + Ifges(5,2)) * t119) * t82) * t85 + (-t110 + t58) * t24 + (-t60 - t109) * t23;
t112 = t4 * qJD(1);
t108 = -t141 / 0.2e1;
t7 = (t124 / 0.2e1 - t60 / 0.2e1 + t84 * t108) * t84 + (t125 / 0.2e1 + t58 / 0.2e1 + t82 * t108) * t82;
t111 = t7 * qJD(1);
t106 = -t86 * mrSges(5,3) / 0.2e1;
t68 = Ifges(5,1) * t84 - t139;
t67 = -Ifges(5,2) * t82 + t138;
t92 = t145 * t67 + t148 * t68;
t14 = qJ(3) * t102 + t100 * t148 + t101 * t147 - t92;
t90 = t70 * t102;
t91 = -t120 / 0.2e1 + t128 / 0.2e1;
t93 = t25 * mrSges(5,1) / 0.2e1 - t26 * mrSges(5,2) / 0.2e1;
t2 = -t90 / 0.2e1 + t91 * t86 + t99 * t83 + (t151 + (mrSges(5,2) * t143 + t68 / 0.4e1 + (Ifges(5,1) / 0.4e1 + t152 + t106) * t84) * t84 + (-0.3e1 / 0.2e1 * t138 + mrSges(5,1) * t143 - t67 / 0.4e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.4e1 + t106) * t82) * t82) * t85 + t93;
t96 = t2 * qJD(1) - t14 * qJD(2);
t51 = mrSges(4,3) + (m(4) + m(5)) * qJ(3) + t65;
t9 = (t117 / 0.2e1 - t57 / 0.2e1) * t84 + (-t116 / 0.2e1 - t59 / 0.2e1) * t82 + 0.2e1 * (t70 / 0.4e1 - t121 / 0.4e1 - t129 / 0.4e1) * m(5);
t95 = qJD(1) * t9 + qJD(2) * t51;
t8 = (t122 / 0.2e1 + t131 / 0.2e1) * t83 - t91 + t156;
t5 = t59 * t148 + t57 * t145 + (-t130 / 0.2e1 + t123 / 0.2e1 + m(4) * pkin(5) + mrSges(4,1)) * t85 + (t70 - t157) * m(5) / 0.2e1;
t3 = -t68 * t119 / 0.2e1 + t67 * t127 / 0.2e1 - t89 * t143 + t90 / 0.2e1 + t60 * t118 / 0.2e1 - t58 * t126 / 0.2e1 + t85 * t151 + t93 - (t101 * t159 + t136) * t82 / 0.4e1 - (t100 * t159 + t133) * t84 / 0.4e1 + (-t99 / 0.4e1 + t94) * t83 + t86 * t156;
t10 = [qJD(2) * t1 + qJD(3) * t6 - qJD(4) * t4, t5 * qJD(3) + t3 * qJD(4) + t113 + (t41 * t145 + t40 * t150 - qJ(3) * t52 + m(5) * (-qJ(3) * t69 - t157 * t86) - t69 * t65 + t59 * t126 + t57 * t118 + (-Ifges(4,4) + Ifges(3,5) + t135 / 0.2e1 - t134 / 0.2e1 - pkin(2) * mrSges(4,1)) * t85 + (-qJ(3) * mrSges(4,1) + Ifges(4,5) - Ifges(3,6) + t92) * t83 + (m(4) * t98 - t85 * mrSges(3,1) + t83 * mrSges(3,2) + t64) * pkin(5) + t157 * mrSges(5,3)) * qJD(2), qJD(2) * t5 + qJD(4) * t8 + t114, -t112 + t3 * qJD(2) + t8 * qJD(3) + (-mrSges(5,1) * t24 - mrSges(5,2) * t23 + (t134 - t135) * t85) * qJD(4); qJD(3) * t9 - qJD(4) * t2 - t113, qJD(3) * t51 + qJD(4) * t14, t95, ((-mrSges(5,2) * t86 - Ifges(5,6)) * t84 + (-mrSges(5,1) * t86 - Ifges(5,5)) * t82) * qJD(4) - t96; -qJD(2) * t9 - qJD(4) * t7 - t114, -t95, 0, -qJD(4) * t65 - t111; qJD(2) * t2 + qJD(3) * t7 + t112, t96, t111, 0;];
Cq = t10;
