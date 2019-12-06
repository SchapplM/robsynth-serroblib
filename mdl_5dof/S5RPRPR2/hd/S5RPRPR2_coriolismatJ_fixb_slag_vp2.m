% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:24
% EndTime: 2019-12-05 17:49:25
% DurationCPUTime: 0.63s
% Computational Cost: add. (2660->94), mult. (5222->129), div. (0->0), fcn. (4838->8), ass. (0->70)
t76 = sin(pkin(9));
t78 = cos(pkin(9));
t101 = t76 ^ 2 + t78 ^ 2;
t119 = t101 * mrSges(5,3);
t123 = -mrSges(4,2) + t119;
t111 = pkin(1) * sin(pkin(8));
t69 = cos(pkin(8)) * pkin(1) + pkin(2);
t80 = sin(qJ(3));
t82 = cos(qJ(3));
t55 = -t80 * t111 + t69 * t82;
t79 = sin(qJ(5));
t81 = cos(qJ(5));
t64 = t76 * t81 + t78 * t79;
t32 = t64 * t55;
t63 = -t76 * t79 + t78 * t81;
t33 = t63 * t55;
t16 = -t32 * t63 + t33 * t64;
t122 = m(6) * t16;
t56 = t82 * t111 + t80 * t69;
t54 = qJ(4) + t56;
t96 = t101 * t54;
t121 = m(5) * t96;
t45 = t64 * mrSges(6,1) + mrSges(6,2) * t63;
t120 = (qJD(1) + qJD(3)) * t45;
t49 = (-pkin(7) - t54) * t76;
t73 = t78 * pkin(7);
t50 = t54 * t78 + t73;
t22 = t49 * t81 - t50 * t79;
t23 = t49 * t79 + t50 * t81;
t105 = -t22 * t64 + t23 * t63;
t65 = (-pkin(7) - qJ(4)) * t76;
t67 = qJ(4) * t78 + t73;
t47 = t65 * t81 - t67 * t79;
t48 = t65 * t79 + t67 * t81;
t103 = -t47 * t64 + t48 * t63;
t118 = -mrSges(5,1) * t78 - mrSges(6,1) * t63 + mrSges(5,2) * t76 + mrSges(6,2) * t64 - mrSges(4,1);
t117 = (t32 * t64 + t33 * t63) * mrSges(6,3);
t115 = t64 ^ 2;
t114 = -m(6) / 0.2e1;
t110 = t78 * pkin(4);
t102 = Ifges(6,5) * t63 - Ifges(6,6) * t64;
t92 = t119 + (t63 ^ 2 + t115) * mrSges(6,3);
t11 = m(6) * t105 + t121 + t92;
t100 = qJD(1) * t11;
t43 = t45 * qJD(5);
t98 = qJD(1) * t122;
t95 = t101 * qJ(4);
t94 = -pkin(3) - t55;
t93 = 0.2e1 * (m(6) / 0.4e1 + m(5) / 0.4e1) * t56;
t15 = m(5) * t95 + m(6) * t103 + t92;
t83 = -m(5) * (t96 + t95) / 0.2e1 + (t103 + t105) * t114 - t92;
t7 = t93 + t83;
t91 = -qJD(1) * t7 + qJD(3) * t15;
t90 = -t32 * mrSges(6,1) / 0.2e1 - t33 * mrSges(6,2) / 0.2e1;
t51 = t94 - t110;
t3 = m(6) * (-t22 * t32 + t23 * t33) + t117 + (m(5) * t94 + m(6) * t51 + t118) * t56 + (t121 + t123) * t55;
t88 = qJD(2) * t114 * t16 - t3 * qJD(1);
t85 = (Ifges(6,4) * t63 + (Ifges(6,1) - Ifges(6,2)) * t64) * t63 - t115 * Ifges(6,4);
t4 = t51 * t45 + t85;
t87 = t4 * qJD(1);
t70 = -pkin(3) - t110;
t84 = (t51 / 0.2e1 + t70 / 0.2e1) * t45 + t85;
t1 = t84 - t90;
t6 = t70 * t45 + t85;
t86 = -t1 * qJD(1) - t6 * qJD(3);
t44 = t45 * qJD(4);
t8 = t93 - t83;
t5 = qJD(3) * t122 / 0.2e1;
t2 = t84 + t90;
t9 = [qJD(3) * t3 + qJD(4) * t11 + qJD(5) * t4, t5, t8 * qJD(4) + t2 * qJD(5) - t88 + (t117 + t118 * t56 + t123 * t55 + (-t32 * t47 + t33 * t48 + t56 * t70) * m(6) + m(5) * (-pkin(3) * t56 + t55 * t95)) * qJD(3), qJD(3) * t8 + t100, t2 * qJD(3) + (-mrSges(6,1) * t23 - mrSges(6,2) * t22 + t102) * qJD(5) + t87; t5, 0, t98 / 0.2e1, 0, -t43; -t7 * qJD(4) + t1 * qJD(5) + t88, -t98 / 0.2e1, qJD(4) * t15 + qJD(5) * t6, t91, (-mrSges(6,1) * t48 - mrSges(6,2) * t47 + t102) * qJD(5) - t86; qJD(3) * t7 - t100 + t43, 0, -t91 + t43, 0, t120; -t1 * qJD(3) - t44 - t87, 0, -t44 + t86, -t120, 0;];
Cq = t9;
