% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:37
% EndTime: 2019-12-05 15:37:40
% DurationCPUTime: 0.74s
% Computational Cost: add. (1389->122), mult. (3475->170), div. (0->0), fcn. (3382->6), ass. (0->71)
t125 = m(6) / 0.2e1;
t85 = cos(qJ(2));
t120 = sin(qJ(4));
t121 = cos(qJ(4));
t82 = sin(pkin(8));
t83 = cos(pkin(8));
t88 = t120 * t83 + t121 * t82;
t53 = t88 * t85;
t66 = t120 * t82 - t121 * t83;
t75 = -pkin(3) * t83 - pkin(2);
t27 = pkin(4) * t66 - qJ(5) * t88 + t75;
t43 = mrSges(6,1) * t66 - mrSges(6,3) * t88;
t135 = m(6) * t27 + t43;
t55 = t66 * t85;
t134 = (mrSges(5,2) / 0.2e1 - mrSges(6,3) / 0.2e1) * t55 + (-pkin(4) * t53 - qJ(5) * t55) * t125;
t133 = t88 ^ 2;
t132 = m(5) + m(6);
t130 = mrSges(5,1) + mrSges(6,1);
t129 = -Ifges(6,5) + Ifges(5,4);
t79 = m(6) * qJ(5) + mrSges(6,3);
t127 = 0.2e1 * m(6);
t126 = m(5) / 0.2e1;
t124 = -t53 / 0.2e1;
t84 = sin(qJ(2));
t123 = t84 / 0.2e1;
t122 = t88 * pkin(4);
t119 = t53 * t88;
t118 = t66 * mrSges(6,2);
t117 = t84 * t85;
t115 = pkin(6) + qJ(3);
t80 = t82 ^ 2;
t81 = t83 ^ 2;
t114 = t80 + t81;
t52 = t88 * t84;
t54 = t66 * t84;
t6 = m(4) * (-0.1e1 + t114) * t117 + t132 * (t52 * t53 + t54 * t55 - t117);
t113 = t6 * qJD(1);
t112 = t66 * qJ(5);
t111 = qJD(4) * t79;
t41 = mrSges(6,1) * t88 + t66 * mrSges(6,3);
t42 = mrSges(5,1) * t88 - t66 * mrSges(5,2);
t95 = -t112 - t122;
t11 = (-t95 / 0.4e1 + t122 / 0.4e1 + t112 / 0.4e1) * t127 + t42 + t41;
t109 = t11 * qJD(2);
t29 = m(6) * t88;
t108 = t29 * qJD(2);
t107 = m(4) * t123;
t106 = t53 / 0.2e1;
t103 = t115 * t82;
t102 = t114 * mrSges(4,3);
t70 = t115 * t83;
t44 = t121 * t103 + t120 * t70;
t45 = -t120 * t103 + t121 * t70;
t100 = t44 * t53 - t45 * t55;
t96 = t114 * qJ(3);
t1 = t27 * t41 + t75 * t42 - t129 * t133 + (-(Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t88 + t129 * t66) * t66 - t135 * t95;
t86 = (-t42 / 0.2e1 - t41 / 0.2e1 + t95 * t125) * t85;
t2 = (mrSges(5,1) / 0.2e1 + mrSges(6,1) / 0.2e1) * t53 + t86 - t134;
t94 = t2 * qJD(1) + t1 * qJD(2);
t4 = m(4) * t96 + t102 + t132 * (t44 * t88 - t45 * t66) + (t66 ^ 2 + t133) * (mrSges(5,3) + mrSges(6,2));
t87 = (t126 + t125) * (t52 * t88 + t66 * t54);
t9 = (-m(5) / 0.2e1 - m(6) / 0.2e1 + (t80 / 0.2e1 + t81 / 0.2e1 - 0.1e1 / 0.2e1) * m(4)) * t84 + t87;
t93 = qJD(1) * t9 + qJD(2) * t4;
t12 = t135 * t88;
t23 = (t106 + t124) * m(6);
t92 = -qJD(1) * t23 + qJD(2) * t12;
t22 = m(6) * t106 + t53 * t125;
t16 = m(6) * t45 - t118;
t8 = t114 * t107 + t132 * t123 + t107 + t87;
t3 = t130 * t124 + t134 + t86;
t5 = [qJD(2) * t6, t8 * qJD(3) + t3 * qJD(4) + t22 * qJD(5) + t113 + (mrSges(6,2) * t119 + t55 * t118 + (-mrSges(3,2) + t102) * t85 + (-t83 * mrSges(4,1) + t66 * mrSges(5,1) + t82 * mrSges(4,2) + mrSges(5,2) * t88 - mrSges(3,1) + t43) * t84 + m(4) * (-t84 * pkin(2) + t85 * t96) + 0.2e1 * (t75 * t84 + t100) * t126 + 0.2e1 * (t27 * t84 + t100) * t125 + (t55 * t66 + t119) * mrSges(5,3)) * qJD(2), t8 * qJD(2), t3 * qJD(2) + (t130 * t54 + (mrSges(5,2) - mrSges(6,3)) * t52) * qJD(4) + ((pkin(4) * t54 - qJ(5) * t52) * qJD(4) / 0.2e1 - t54 * qJD(5) / 0.2e1) * t127, -m(6) * t54 * qJD(4) + t22 * qJD(2); qJD(3) * t9 + qJD(4) * t2 + qJD(5) * t23 - t113, qJD(3) * t4 + qJD(4) * t1 - qJD(5) * t12, t93, t16 * qJD(5) + t94 + ((-m(6) * pkin(4) - t130) * t45 + (mrSges(5,2) - t79) * t44 - (qJ(5) * mrSges(6,2) + Ifges(5,6) - Ifges(6,6)) * t88 + (pkin(4) * mrSges(6,2) - Ifges(6,4) - Ifges(5,5)) * t66) * qJD(4), qJD(4) * t16 - t92; -t9 * qJD(2), qJD(4) * t11 - qJD(5) * t29 - t93, 0, t109, -t108; -qJD(2) * t2, -qJD(3) * t11 - t94, -t109, t79 * qJD(5), t111; -t23 * qJD(2), qJD(3) * t29 + t92, t108, -t111, 0;];
Cq = t5;
