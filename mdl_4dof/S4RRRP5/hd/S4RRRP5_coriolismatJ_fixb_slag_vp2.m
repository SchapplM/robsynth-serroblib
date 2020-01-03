% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:39
% EndTime: 2019-12-31 17:16:41
% DurationCPUTime: 0.73s
% Computational Cost: add. (1628->107), mult. (3388->138), div. (0->0), fcn. (2983->4), ass. (0->67)
t125 = cos(qJ(2));
t76 = -t125 * pkin(2) - pkin(1);
t122 = sin(qJ(3));
t123 = sin(qJ(2));
t124 = cos(qJ(3));
t65 = -t122 * t125 - t124 * t123;
t106 = t65 * qJ(4);
t64 = t122 * t123 - t124 * t125;
t92 = t64 * pkin(3) + t106;
t23 = t92 + t76;
t136 = m(5) * t23 + mrSges(5,1) * t64 + mrSges(5,3) * t65;
t81 = t125 * pkin(5);
t110 = t125 * pkin(6) + t81;
t102 = t123 * pkin(5);
t87 = -t123 * pkin(6) - t102;
t130 = t124 * t110 + t122 * t87;
t137 = t130 * mrSges(5,1);
t138 = t130 * mrSges(4,1);
t43 = t122 * t110 - t124 * t87;
t144 = t43 * mrSges(5,3);
t145 = t43 * mrSges(4,2);
t97 = (Ifges(4,6) - Ifges(5,6)) * t65 + (-Ifges(5,4) - Ifges(4,5)) * t64;
t148 = t97 - t137 - t138 + t145 - t144;
t147 = t137 / 0.2e1 + t138 / 0.2e1 + t144 / 0.2e1 - t145 / 0.2e1;
t141 = m(5) * (-pkin(3) * t130 - qJ(4) * t43);
t139 = -Ifges(5,5) + Ifges(4,4);
t78 = m(5) * qJ(4) + mrSges(5,3);
t135 = qJD(3) * t78;
t134 = t78 * qJD(4);
t129 = m(5) / 0.2e1;
t101 = t122 * pkin(2);
t104 = t124 * pkin(2);
t131 = t65 * t101 + t64 * t104;
t128 = m(5) * pkin(2);
t127 = -mrSges(5,2) / 0.2e1;
t126 = m(5) * t130;
t113 = t64 * mrSges(5,2);
t6 = t136 * t65;
t109 = qJD(1) * t6;
t100 = t139 * t65;
t103 = t123 * pkin(2);
t29 = -t65 * pkin(3) + t64 * qJ(4);
t86 = (Ifges(4,1) + Ifges(5,1) - Ifges(4,2) - Ifges(5,3)) * t65 + t139 * t64;
t91 = t23 * (-mrSges(5,1) * t65 + mrSges(5,3) * t64) + t76 * (-mrSges(4,1) * t65 - mrSges(4,2) * t64);
t1 = -pkin(1) * (t123 * mrSges(3,1) + t125 * mrSges(3,2)) + m(4) * t76 * t103 + (-mrSges(4,2) * t103 - t100) * t65 + (mrSges(4,1) * t103 + t86) * t64 + t91 + (-Ifges(3,2) + Ifges(3,1)) * t125 * t123 + (-t123 ^ 2 + t125 ^ 2) * Ifges(3,4) + t136 * (t103 + t29);
t108 = t1 * qJD(1);
t2 = -t100 * t65 + t136 * t29 + t86 * t64 + t91;
t107 = t2 * qJD(1);
t71 = t101 + qJ(4);
t105 = t71 * t65 * mrSges(5,2);
t57 = -t113 / 0.2e1;
t94 = 0.2e1 * t57;
t75 = -t104 - pkin(3);
t84 = (-mrSges(4,1) - mrSges(5,1)) * t101 + (-mrSges(4,2) + mrSges(5,3)) * t104;
t14 = -(t122 * t75 + t124 * t71) * t128 - t84;
t82 = ((t101 - t71) * t43 + (t75 + t104) * t130) * t129 + t105 / 0.2e1 + t75 * t57 + t131 * t127 - t147;
t83 = -t141 / 0.2e1 + pkin(3) * t57 + t106 * t127 + t147;
t4 = t82 + t83;
t90 = t4 * qJD(1) - t14 * qJD(2);
t68 = m(5) * t71 + mrSges(5,3);
t89 = qJD(2) * t68;
t88 = -qJD(2) * t78 - t135;
t55 = 0.2e1 * t101 * t129 + t78;
t10 = t94 + t126;
t7 = t126 / 0.2e1 + t130 * t129 + t94;
t3 = t82 - t83 + t97;
t5 = [qJD(2) * t1 + qJD(3) * t2 + qJD(4) * t6, t108 + (-t75 * t113 + t105 + mrSges(3,2) * t102 - mrSges(3,1) * t81 + m(5) * (t130 * t75 - t71 * t43) + m(4) * (-t122 * t43 - t124 * t130) * pkin(2) + Ifges(3,5) * t125 - Ifges(3,6) * t123 + t131 * mrSges(4,3) + t148) * qJD(2) + t3 * qJD(3) + t7 * qJD(4), t107 + t3 * qJD(2) + (t92 * mrSges(5,2) + t141 + t148) * qJD(3) + t10 * qJD(4), qJD(2) * t7 + qJD(3) * t10 + t109; qJD(3) * t4 - t108, -qJD(3) * t14 + qJD(4) * t68, ((-t122 * pkin(3) + t124 * qJ(4)) * t128 + t84) * qJD(3) + t55 * qJD(4) + t90, qJD(3) * t55 + t89; -qJD(2) * t4 - t107, -t90 + t134, t134, -t88; -t109, -t89 - t135, t88, 0;];
Cq = t5;
