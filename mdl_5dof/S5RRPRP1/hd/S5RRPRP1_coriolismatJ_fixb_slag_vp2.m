% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:57
% EndTime: 2020-01-03 11:58:58
% DurationCPUTime: 0.70s
% Computational Cost: add. (1465->123), mult. (2959->149), div. (0->0), fcn. (2151->6), ass. (0->88)
t108 = Ifges(5,4) + Ifges(6,4);
t70 = sin(qJ(4));
t134 = t108 * t70;
t144 = Ifges(5,2) + Ifges(6,2);
t72 = cos(qJ(4));
t147 = t144 * t72 + t134;
t146 = t70 / 0.2e1;
t145 = Ifges(5,1) + Ifges(6,1);
t67 = t72 ^ 2;
t143 = t70 ^ 2 + t67;
t142 = (mrSges(5,2) + mrSges(6,2)) * t72;
t130 = m(6) * pkin(4);
t96 = mrSges(6,1) + t130;
t140 = -mrSges(5,1) - t96;
t106 = t143 * mrSges(6,3);
t109 = t72 * mrSges(6,2);
t47 = t96 * t70 + t109;
t138 = (qJD(1) + qJD(2)) * t47;
t69 = cos(pkin(8));
t71 = sin(qJ(2));
t110 = t69 * t71;
t73 = cos(qJ(2));
t125 = pkin(1) * t73;
t61 = pkin(2) + t125;
t68 = sin(pkin(8));
t90 = pkin(1) * t110 + t68 * t61;
t37 = pkin(7) + t90;
t104 = qJ(5) + t37;
t26 = t104 * t72;
t22 = t26 * t72;
t25 = t104 * t70;
t137 = t25 * t70 + t22;
t60 = pkin(2) * t68 + pkin(7);
t103 = qJ(5) + t60;
t46 = t103 * t72;
t38 = t46 * t72;
t45 = t103 * t70;
t136 = t45 * t70 + t38;
t133 = t145 * t72;
t126 = pkin(1) * t71;
t59 = t68 * t126;
t44 = t69 * t125 - t59;
t92 = t143 * t44;
t131 = m(6) / 0.2e1;
t129 = -t70 / 0.2e1;
t124 = pkin(4) * t70;
t123 = t72 * pkin(4);
t122 = mrSges(5,1) * t72;
t107 = -Ifges(5,6) - Ifges(6,6);
t43 = (t68 * t73 + t110) * pkin(1);
t86 = -t72 * mrSges(6,1) + t70 * mrSges(6,2);
t75 = -mrSges(3,1) * t126 - mrSges(3,2) * t125 + mrSges(5,3) * t92 + (-mrSges(4,2) + t106) * t44 + (mrSges(5,2) * t70 - mrSges(4,1) - t122 + t86) * t43;
t91 = t61 * t69 - t59;
t89 = -pkin(3) - t91;
t81 = t89 - t123;
t3 = m(6) * (t137 * t44 + t81 * t43) + m(5) * (t37 * t92 + t89 * t43) + m(4) * (-t91 * t43 + t90 * t44) + t75;
t105 = t3 * qJD(1);
t12 = m(6) * t137 + t106;
t102 = qJD(1) * t12;
t101 = qJD(4) * t70;
t99 = t43 * t131;
t95 = -pkin(2) * t69 - pkin(3);
t56 = t86 * t124;
t88 = -t108 * t67 - t56;
t87 = t70 * mrSges(5,1) + t72 * mrSges(5,2);
t85 = t70 * mrSges(6,1) + t109;
t84 = -mrSges(6,3) * t123 + (Ifges(5,5) + Ifges(6,5)) * t72;
t15 = m(6) * t136 + t106;
t77 = -(t22 + t38 + (t25 + t45) * t70) * m(6) / 0.2e1 - t106;
t8 = t99 + t77;
t83 = -qJD(1) * t8 + qJD(2) * t15;
t82 = t95 - t123;
t21 = t81 * t85;
t24 = t89 * t87;
t76 = -t133 + t147;
t4 = -t21 - t24 + (-t81 * t130 + t76) * t70 + t88;
t79 = t4 * qJD(1);
t35 = t82 * t85;
t42 = t95 * t87;
t74 = t56 + t21 / 0.2e1 + t24 / 0.2e1 + t35 / 0.2e1 + t42 / 0.2e1 + (-0.2e1 * t123 - 0.2e1 * pkin(3) + t59 + (-pkin(2) - t61) * t69) * t124 * t131 + t147 * t129 + (t133 - t134) * t146 + ((-t144 + t145) * t146 + t108 * t72) * t72;
t1 = -t74 + (-t142 / 0.2e1 - t140 * t129) * t44;
t5 = -t35 - t42 + (-t82 * t130 + t76) * t70 + t88;
t78 = t1 * qJD(1) + t5 * qJD(2);
t41 = t47 * qJD(5);
t40 = t47 * qJD(4);
t9 = t99 - t77;
t2 = t74 + ((-mrSges(5,2) / 0.2e1 - mrSges(6,2) / 0.2e1) * t72 + (-mrSges(5,1) / 0.2e1 - mrSges(6,1) / 0.2e1 - t130 / 0.2e1) * t70) * t44;
t6 = [qJD(2) * t3 - qJD(4) * t4 + qJD(5) * t12, t105 + (m(6) * (t136 * t44 + t82 * t43) + m(5) * (t95 * t43 + t60 * t92) + m(4) * (-t43 * t69 + t44 * t68) * pkin(2) + t75) * qJD(2) + t2 * qJD(4) + t9 * qJD(5), 0, t2 * qJD(2) + (mrSges(6,2) * t25 - t37 * t122 - t96 * t26 + t84) * qJD(4) + (mrSges(5,2) * t37 + t107) * t101 - t79, qJD(2) * t9 + t102; -qJD(4) * t1 - qJD(5) * t8 - t105, -qJD(4) * t5 + qJD(5) * t15, 0, (mrSges(6,2) * t45 - t60 * t122 - t96 * t46 + t84) * qJD(4) + (mrSges(5,2) * t60 + t107) * t101 - t78, t83; 0, 0, 0, (t140 * t70 - t142) * qJD(4), 0; t1 * qJD(2) - t41 + t79, -t41 + t78, 0, 0, -t138; qJD(2) * t8 - t102 + t40, t40 - t83, 0, t138, 0;];
Cq = t6;
