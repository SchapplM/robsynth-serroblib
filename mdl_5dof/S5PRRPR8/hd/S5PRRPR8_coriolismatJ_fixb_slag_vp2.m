% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:22
% EndTime: 2019-12-31 17:42:23
% DurationCPUTime: 0.57s
% Computational Cost: add. (1410->100), mult. (3223->151), div. (0->0), fcn. (3355->8), ass. (0->75)
t115 = sin(qJ(2));
t116 = cos(qJ(2));
t77 = sin(qJ(3));
t79 = cos(qJ(3));
t56 = -t77 * t115 + t79 * t116;
t57 = -t79 * t115 - t77 * t116;
t74 = sin(pkin(9));
t75 = cos(pkin(9));
t36 = t75 * t56 + t57 * t74;
t76 = sin(qJ(5));
t78 = cos(qJ(5));
t99 = t76 ^ 2 + t78 ^ 2;
t132 = t36 * t99;
t131 = t78 / 0.2e1;
t129 = -Ifges(6,2) * t78 / 0.2e1 - Ifges(6,4) * t76 + Ifges(6,1) * t131;
t58 = -mrSges(6,1) * t78 + mrSges(6,2) * t76;
t88 = t56 * t74 - t75 * t57;
t128 = t58 * t88;
t104 = t74 * t77;
t119 = pkin(2) * t79;
t69 = pkin(3) + t119;
t44 = -pkin(2) * t104 + t69 * t75;
t41 = -pkin(4) - t44;
t101 = t78 * mrSges(6,2);
t102 = t76 * mrSges(6,1);
t59 = t101 + t102;
t107 = t41 * t59;
t71 = Ifges(6,4) * t78;
t61 = -Ifges(6,2) * t76 + t71;
t62 = Ifges(6,1) * t76 + t71;
t89 = t129 * t76 + (t61 + t62) * t131;
t18 = t89 + t107;
t110 = t36 * t59;
t85 = -t101 / 0.2e1 - t102 / 0.2e1;
t83 = t85 * t36;
t16 = t110 / 0.2e1 + t83;
t97 = t16 * qJD(1);
t127 = t18 * qJD(2) - t97;
t66 = -pkin(3) * t75 - pkin(4);
t105 = t66 * t59;
t20 = t89 + t105;
t50 = (t75 * t79 - t104) * pkin(2);
t121 = -t50 / 0.2e1;
t9 = (-t41 / 0.2e1 - t66 / 0.2e1) * t59 + (mrSges(6,2) * t121 - t62 / 0.2e1 - t61 / 0.2e1) * t78 + (mrSges(6,1) * t121 - t129) * t76;
t126 = t9 * qJD(2) - t20 * qJD(3) + t97;
t103 = t75 * t77;
t49 = (t74 * t79 + t103) * pkin(2);
t125 = (t99 * mrSges(6,3) - mrSges(5,2)) * t50 + (-mrSges(5,1) + t58) * t49 + (-mrSges(4,1) * t77 - mrSges(4,2) * t79) * pkin(2);
t124 = 0.2e1 * qJD(3);
t123 = m(5) / 0.2e1;
t122 = m(6) / 0.2e1;
t5 = m(6) * (-t36 + t132) * t88;
t98 = t5 * qJD(1);
t118 = -t16 * qJD(5) - t98;
t17 = -t110 / 0.2e1 + t83;
t117 = t17 * qJD(5) + t98;
t106 = t49 * t36;
t45 = pkin(2) * t103 + t74 * t69;
t100 = t45 * t36 - t44 * t88;
t95 = qJD(2) + qJD(3);
t94 = pkin(3) * t123;
t42 = pkin(7) + t45;
t93 = t132 * t42 + t41 * t88;
t91 = t99 * t50;
t90 = Ifges(6,5) * t78 - Ifges(6,6) * t76;
t80 = (t88 * t91 - t106 + t93) * t122 + (t50 * t88 + t100 - t106) * t123;
t65 = pkin(3) * t74 + pkin(7);
t81 = (t132 * t65 + t66 * t88) * t122 + (t36 * t74 - t75 * t88) * t94;
t2 = t80 - t81;
t7 = m(6) * (t41 * t49 + t42 * t91) + m(5) * (-t44 * t49 + t45 * t50) + t125;
t87 = t2 * qJD(1) + t7 * qJD(2);
t82 = t57 * mrSges(4,1) - mrSges(5,1) * t88 - t56 * mrSges(4,2) - t36 * mrSges(5,2) + mrSges(6,3) * t132 + t128;
t10 = t107 / 0.2e1 + t105 / 0.2e1 + t85 * t50 + t89;
t1 = t80 + t81 + t82;
t3 = [t95 * t5, t1 * qJD(3) + t117 + (-t115 * mrSges(3,1) - t116 * mrSges(3,2) + t82 + 0.2e1 * t93 * t122 + m(4) * (pkin(2) * t56 * t77 + t57 * t119) + 0.2e1 * t100 * t123) * qJD(2), t1 * qJD(2) + t82 * qJD(3) + t81 * t124 + t117, 0, qJD(5) * t128 + t95 * t17; qJD(3) * t2 + t118, qJD(3) * t7 + qJD(5) * t18, t10 * qJD(5) + ((t66 * t49 + t65 * t91) * t122 + (-t49 * t75 + t50 * t74) * t94) * t124 + t87 + t125 * qJD(3), 0, t10 * qJD(3) + (t58 * t42 + t90) * qJD(5) + t127; -qJD(2) * t2 + t118, -qJD(5) * t9 - t87, t20 * qJD(5), 0, (t58 * t65 + t90) * qJD(5) - t126; 0, 0, 0, 0, -t59 * qJD(5); t95 * t16, qJD(3) * t9 - t127, t126, 0, 0;];
Cq = t3;
