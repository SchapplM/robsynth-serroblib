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
% m [6x1]
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:19:35
% EndTime: 2022-01-20 10:19:38
% DurationCPUTime: 0.73s
% Computational Cost: add. (1465->123), mult. (2959->150), div. (0->0), fcn. (2151->6), ass. (0->88)
t106 = Ifges(6,4) + Ifges(5,4);
t68 = sin(qJ(4));
t133 = t106 * t68;
t143 = Ifges(5,2) + Ifges(6,2);
t70 = cos(qJ(4));
t146 = t143 * t70 + t133;
t145 = t68 / 0.2e1;
t144 = Ifges(5,1) + Ifges(6,1);
t65 = t70 ^ 2;
t104 = t68 ^ 2 + t65;
t142 = (mrSges(5,2) + mrSges(6,2)) * t70;
t107 = t70 * mrSges(6,2);
t130 = m(6) * pkin(4);
t95 = mrSges(6,1) + t130;
t47 = t95 * t68 + t107;
t140 = (qJD(1) + qJD(2)) * t47;
t139 = -mrSges(5,1) - t95;
t67 = cos(pkin(8));
t69 = sin(qJ(2));
t109 = t67 * t69;
t71 = cos(qJ(2));
t122 = t71 * pkin(1);
t61 = pkin(2) + t122;
t66 = sin(pkin(8));
t88 = pkin(1) * t109 + t66 * t61;
t37 = pkin(7) + t88;
t102 = qJ(5) + t37;
t26 = t102 * t70;
t22 = t26 * t70;
t25 = t102 * t68;
t137 = t25 * t68 + t22;
t60 = t66 * pkin(2) + pkin(7);
t101 = qJ(5) + t60;
t46 = t101 * t70;
t38 = t46 * t70;
t45 = t101 * t68;
t136 = t45 * t68 + t38;
t132 = t144 * t70;
t131 = m(6) / 0.2e1;
t129 = -t68 / 0.2e1;
t43 = (t66 * t71 + t109) * pkin(1);
t126 = m(6) * t43;
t125 = pkin(1) * t69;
t124 = t68 * pkin(4);
t123 = t70 * pkin(4);
t108 = t70 * mrSges(5,1);
t105 = -Ifges(5,6) - Ifges(6,6);
t59 = t66 * t125;
t44 = t67 * t122 - t59;
t84 = -t70 * mrSges(6,1) + t68 * mrSges(6,2);
t73 = -mrSges(3,1) * t125 - mrSges(3,2) * t122 + (t68 * mrSges(5,2) - mrSges(4,1) - t108 + t84) * t43 + (-mrSges(4,2) + (mrSges(5,3) + mrSges(6,3)) * t104) * t44;
t89 = t67 * t61 - t59;
t87 = -pkin(3) - t89;
t79 = t87 - t123;
t90 = t104 * t44;
t3 = m(6) * (t137 * t44 + t79 * t43) + m(5) * (t37 * t90 + t87 * t43) + m(4) * (-t89 * t43 + t88 * t44) + t73;
t103 = t3 * qJD(1);
t100 = qJD(4) * t68;
t91 = t104 * mrSges(6,3);
t12 = m(6) * t137 + t91;
t99 = t12 * qJD(1);
t94 = -t67 * pkin(2) - pkin(3);
t56 = t84 * t124;
t86 = -t106 * t65 - t56;
t85 = t68 * mrSges(5,1) + t70 * mrSges(5,2);
t83 = t68 * mrSges(6,1) + t107;
t82 = -mrSges(6,3) * t123 + (Ifges(5,5) + Ifges(6,5)) * t70;
t15 = m(6) * t136 + t91;
t75 = t91 + (t22 + t38 + (t25 + t45) * t68) * t131;
t9 = -t126 / 0.2e1 + t75;
t81 = t9 * qJD(1) + t15 * qJD(2);
t80 = t94 - t123;
t21 = t79 * t83;
t24 = t87 * t85;
t74 = -t132 + t146;
t4 = -t21 - t24 + (-t79 * t130 + t74) * t68 + t86;
t77 = t4 * qJD(1);
t35 = t80 * t83;
t42 = t94 * t85;
t72 = t56 + t21 / 0.2e1 + t24 / 0.2e1 + t35 / 0.2e1 + t42 / 0.2e1 + (-0.2e1 * t123 - 0.2e1 * pkin(3) + t59 + (-pkin(2) - t61) * t67) * t124 * t131 + t146 * t129 + (t132 - t133) * t145 + ((-t143 + t144) * t145 + t106 * t70) * t70;
t1 = -t72 + (-t142 / 0.2e1 - t139 * t129) * t44;
t5 = -t35 - t42 + (-t80 * t130 + t74) * t68 + t86;
t76 = t1 * qJD(1) + t5 * qJD(2);
t41 = t47 * qJD(5);
t40 = t47 * qJD(4);
t8 = t126 / 0.2e1 + t75;
t2 = t72 + ((-mrSges(5,2) / 0.2e1 - mrSges(6,2) / 0.2e1) * t70 + (-mrSges(5,1) / 0.2e1 - t130 / 0.2e1 - mrSges(6,1) / 0.2e1) * t68) * t44;
t6 = [t3 * qJD(2) - t4 * qJD(4) + t12 * qJD(5), t103 + (m(6) * (t136 * t44 + t80 * t43) + m(5) * (t94 * t43 + t60 * t90) + m(4) * (-t67 * t43 + t44 * t66) * pkin(2) + t73) * qJD(2) + t2 * qJD(4) + t8 * qJD(5), 0, t2 * qJD(2) + (t25 * mrSges(6,2) - t37 * t108 - t95 * t26 + t82) * qJD(4) + (mrSges(5,2) * t37 + t105) * t100 - t77, t8 * qJD(2) + t99; -t1 * qJD(4) + t9 * qJD(5) - t103, -t5 * qJD(4) + t15 * qJD(5), 0, (t45 * mrSges(6,2) - t60 * t108 - t95 * t46 + t82) * qJD(4) + (mrSges(5,2) * t60 + t105) * t100 - t76, t81; 0, 0, 0, (t139 * t68 - t142) * qJD(4), 0; t1 * qJD(2) - t41 + t77, -t41 + t76, 0, 0, -t140; -t9 * qJD(2) + t40 - t99, t40 - t81, 0, t140, 0;];
Cq = t6;
