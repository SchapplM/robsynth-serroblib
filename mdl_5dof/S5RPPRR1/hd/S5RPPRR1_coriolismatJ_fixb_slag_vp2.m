% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:00
% EndTime: 2019-12-05 17:38:02
% DurationCPUTime: 0.52s
% Computational Cost: add. (1560->87), mult. (2765->116), div. (0->0), fcn. (2460->4), ass. (0->58)
t71 = qJ(2) - pkin(6);
t76 = sin(qJ(4));
t108 = t76 * t71;
t64 = -t76 * pkin(7) + t108;
t78 = cos(qJ(4));
t65 = (-pkin(7) + t71) * t78;
t75 = sin(qJ(5));
t77 = cos(qJ(5));
t38 = t77 * t64 + t75 * t65;
t59 = t75 * t76 - t77 * t78;
t60 = -t75 * t78 - t77 * t76;
t88 = -t75 * t64 + t77 * t65;
t3 = -t38 * mrSges(6,1) - t88 * mrSges(6,2) + Ifges(6,5) * t60 + Ifges(6,6) * t59;
t136 = t3 * qJD(5);
t89 = t60 * mrSges(6,1) + t59 * mrSges(6,2);
t135 = qJD(5) * t89;
t134 = t60 ^ 2;
t122 = m(6) * pkin(4);
t72 = pkin(1) + qJ(3);
t66 = t76 * pkin(4) + t72;
t133 = m(6) * t66;
t130 = -t75 * t59 + t77 * t60;
t112 = t60 * mrSges(6,2);
t114 = t59 * mrSges(6,1);
t2 = (-Ifges(6,4) * t59 + (-Ifges(6,1) + Ifges(6,2)) * t60) * t59 + Ifges(6,4) * t134 + t66 * (t112 - t114);
t128 = t59 ^ 2 + t134;
t85 = t114 / 0.2e1 - t112 / 0.2e1;
t124 = t78 ^ 2;
t123 = -m(6) / 0.2e1;
t119 = t78 * pkin(4);
t113 = t59 * t77;
t111 = t60 * t75;
t106 = t78 * mrSges(5,2);
t101 = t2 * qJD(1);
t91 = t76 ^ 2 + t124;
t87 = m(5) * t91;
t11 = -mrSges(4,2) - mrSges(3,3) + t91 * mrSges(5,3) - m(6) * (-t60 * t38 - t59 * t88) - t71 * t87 + (-m(4) - m(3)) * qJ(2) + t128 * mrSges(6,3);
t99 = t11 * qJD(1);
t90 = t78 * mrSges(5,1) - t76 * mrSges(5,2);
t12 = (-t113 / 0.2e1 - t111 / 0.2e1 + t78 / 0.2e1) * t122 - 0.2e1 * t85 + t90;
t98 = t12 * qJD(1);
t83 = -t76 * mrSges(5,1) - t106 + t89;
t17 = mrSges(4,3) + t133 + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t72 - t83;
t97 = t17 * qJD(1);
t18 = 0.2e1 * t85;
t96 = t18 * qJD(1);
t81 = t87 / 0.2e1 + m(6) * t128 / 0.2e1;
t94 = -m(5) / 0.2e1 + t123;
t23 = -m(4) - t81 + t94;
t95 = t23 * qJD(1);
t1 = t72 * t90 - t124 * Ifges(5,4) + (Ifges(5,4) * t76 + (-Ifges(5,1) + Ifges(5,2)) * t78) * t76 + t2 + (-t89 + t133) * t119;
t86 = t1 * qJD(1);
t63 = (mrSges(6,1) * t75 + mrSges(6,2) * t77) * pkin(4);
t82 = t63 * qJD(4);
t58 = t63 * qJD(5);
t24 = t81 + t94;
t13 = t119 * t123 + (-t111 - t113) * t122 / 0.2e1;
t4 = [-t11 * qJD(2) + t17 * qJD(3) + t1 * qJD(4) + t2 * qJD(5), t24 * qJD(3) + t13 * qJD(4) - t99, t24 * qJD(2) + t97, t13 * qJD(2) + t136 + t86 + (-mrSges(5,1) * t108 - Ifges(5,5) * t76 - Ifges(5,6) * t78 - t71 * t106 + (m(6) * (-t38 * t77 + t75 * t88) - t130 * mrSges(6,3)) * pkin(4) + t3) * qJD(4), t3 * qJD(4) + t101 + t136; t23 * qJD(3) - t12 * qJD(4) + t18 * qJD(5) + t99, 0, t95, -t98, t96; -t23 * qJD(2) - t97, -t95, 0, (t130 * t122 + t83) * qJD(4) + t135, qJD(4) * t89 + t135; t12 * qJD(2) - t86, t98, 0, -t58, -t58 - t82; -t18 * qJD(2) - t101, -t96, 0, t82, 0;];
Cq = t4;
