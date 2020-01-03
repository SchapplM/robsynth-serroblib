% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:44
% EndTime: 2019-12-31 18:40:45
% DurationCPUTime: 0.47s
% Computational Cost: add. (1041->115), mult. (1995->148), div. (0->0), fcn. (1432->6), ass. (0->75)
t69 = sin(qJ(4));
t71 = cos(qJ(4));
t115 = t69 ^ 2 + t71 ^ 2;
t104 = pkin(1) * sin(pkin(8));
t57 = cos(pkin(8)) * pkin(1) + pkin(2);
t70 = sin(qJ(3));
t72 = cos(qJ(3));
t30 = -t70 * t104 + t72 * t57;
t114 = t115 * t30;
t88 = t69 * qJ(5);
t80 = -t71 * pkin(4) - t88;
t41 = -pkin(3) + t80;
t18 = -t30 + t41;
t113 = t41 + t18;
t43 = -t71 * mrSges(6,1) - t69 * mrSges(6,3);
t112 = -t71 * mrSges(5,1) + t69 * mrSges(5,2) + t43;
t31 = t72 * t104 + t70 * t57;
t111 = (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t115) * t30 + (-mrSges(4,1) + t112) * t31;
t110 = t30 / 0.2e1;
t45 = -t69 * pkin(4) + t71 * qJ(5);
t106 = m(6) * t45;
t105 = m(6) * t71;
t95 = t71 * mrSges(5,2);
t97 = t69 * mrSges(5,1);
t47 = t95 + t97;
t103 = pkin(3) * t47;
t102 = Ifges(5,4) * t69;
t101 = Ifges(6,5) * t71;
t28 = -pkin(3) - t30;
t100 = t28 * t47;
t96 = t69 * mrSges(6,1);
t59 = t71 * mrSges(6,2);
t94 = t71 * mrSges(6,3);
t29 = pkin(7) + t31;
t93 = t114 * t29;
t92 = t114 * pkin(7);
t1 = m(6) * (t18 * t31 + t93) + m(5) * (t28 * t31 + t93) + t111;
t90 = t1 * qJD(1);
t27 = t45 * t43;
t61 = Ifges(6,5) * t69;
t48 = -Ifges(6,3) * t71 + t61;
t49 = Ifges(6,3) * t69 + t101;
t50 = Ifges(5,2) * t71 + t102;
t64 = Ifges(5,4) * t71;
t51 = -Ifges(5,2) * t69 + t64;
t52 = Ifges(6,1) * t69 - t101;
t53 = Ifges(6,1) * t71 + t61;
t54 = Ifges(5,1) * t69 + t64;
t55 = Ifges(5,1) * t71 - t102;
t77 = -t69 * t50 / 0.2e1 - t27 - t71 * t49 / 0.2e1 + (t48 + t53 + t55) * t69 / 0.2e1 + (t51 + t52 + t54) * t71 / 0.2e1;
t46 = -t94 + t96;
t83 = t46 - t106;
t4 = t83 * t18 + t100 + t77;
t89 = t4 * qJD(1);
t11 = (m(6) * t18 + t43) * t69;
t87 = t11 * qJD(1);
t58 = m(6) * qJ(5) + mrSges(6,3);
t86 = t58 * qJD(4);
t85 = -t18 / 0.2e1 - t41 / 0.2e1;
t84 = m(6) * t113;
t81 = -t84 / 0.2e1;
t73 = (-t96 / 0.2e1 + t94 / 0.2e1 - t95 / 0.2e1 - t97 / 0.2e1 + t106 / 0.2e1) * t30;
t2 = t27 + (-t28 / 0.2e1 + pkin(3) / 0.2e1) * t47 + t85 * t46 + t45 * t84 / 0.2e1 + (-t54 / 0.2e1 - t51 / 0.2e1 - t52 / 0.2e1 + t49 / 0.2e1) * t71 + (-t55 / 0.2e1 + t50 / 0.2e1 - t53 / 0.2e1 - t48 / 0.2e1) * t69 + t73;
t5 = t83 * t41 - t103 + t77;
t79 = t2 * qJD(1) - t5 * qJD(3);
t16 = (m(6) * t41 + t43) * t69;
t6 = (t43 + (t110 - t85) * m(6)) * t69;
t78 = t6 * qJD(1) + t16 * qJD(3);
t75 = (-mrSges(6,2) * t88 - pkin(4) * t59 + (Ifges(6,4) + Ifges(5,5)) * t71 + (-Ifges(5,6) + Ifges(6,6)) * t69) * qJD(4);
t74 = qJD(4) * (m(6) * t80 + t112);
t42 = pkin(7) * t105 + t59;
t17 = t29 * t105 + t59;
t7 = (m(6) * t110 - t43 + t81) * t69;
t3 = t77 + t73 + t45 * t81 - t103 / 0.2e1 + t100 / 0.2e1 + t113 * t46 / 0.2e1;
t8 = [t1 * qJD(3) + t4 * qJD(4) - t11 * qJD(5), 0, t3 * qJD(4) + t7 * qJD(5) + t90 + (m(6) * (t41 * t31 + t92) + m(5) * (-pkin(3) * t31 + t92) + t111) * qJD(3), t3 * qJD(3) + t17 * qJD(5) + t29 * t74 + t75 + t89, t7 * qJD(3) + t17 * qJD(4) - t87; 0, 0, 0, (t106 + (-mrSges(5,2) + mrSges(6,3)) * t71) * qJD(4) + ((-mrSges(5,1) - mrSges(6,1)) * qJD(4) + m(6) * qJD(5)) * t69, m(6) * qJD(4) * t69; -t2 * qJD(4) - t6 * qJD(5) - t90, 0, t5 * qJD(4) - t16 * qJD(5), pkin(7) * t74 + t42 * qJD(5) + t75 - t79, t42 * qJD(4) - t78; t2 * qJD(3) - t89, 0, t79, t58 * qJD(5), t86; t6 * qJD(3) + t87, 0, t78, -t86, 0;];
Cq = t8;
