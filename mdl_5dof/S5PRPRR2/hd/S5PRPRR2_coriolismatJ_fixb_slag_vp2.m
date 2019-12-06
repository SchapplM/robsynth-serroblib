% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:50
% EndTime: 2019-12-05 15:44:52
% DurationCPUTime: 0.65s
% Computational Cost: add. (1543->100), mult. (3354->150), div. (0->0), fcn. (3556->8), ass. (0->70)
t102 = cos(qJ(4));
t101 = sin(qJ(2));
t103 = cos(qJ(2));
t85 = sin(pkin(9));
t86 = cos(pkin(9));
t49 = -t86 * t101 - t85 * t103;
t65 = sin(qJ(4));
t69 = -t85 * t101 + t86 * t103;
t34 = t102 * t69 + t49 * t65;
t64 = sin(qJ(5));
t62 = t64 ^ 2;
t66 = cos(qJ(5));
t63 = t66 ^ 2;
t89 = t62 + t63;
t117 = t89 * t34;
t108 = mrSges(6,3) * t117;
t51 = -mrSges(6,1) * t66 + t64 * mrSges(6,2);
t67 = -t102 * t49 + t65 * t69;
t111 = t67 * t51;
t95 = t34 * mrSges(5,2);
t99 = t67 * mrSges(5,1);
t122 = t108 + t111 - t99 - t95;
t90 = t66 * mrSges(6,2);
t91 = t64 * mrSges(6,1);
t73 = t90 / 0.2e1 + t91 / 0.2e1;
t72 = t73 * t34;
t52 = t90 + t91;
t96 = t34 * t52;
t119 = -t96 / 0.2e1 - t72;
t121 = t119 * qJD(5);
t120 = t66 / 0.2e1;
t70 = m(6) * (-pkin(4) * t67 + pkin(7) * t117);
t98 = t34 * t67;
t59 = t86 * pkin(2) + pkin(3);
t81 = pkin(2) * t85;
t41 = t102 * t59 - t65 * t81;
t39 = -pkin(4) - t41;
t78 = t89 * t41;
t116 = t39 + t78;
t42 = t102 * t81 + t65 * t59;
t40 = pkin(7) + t42;
t114 = -t89 * t40 + t42;
t113 = -Ifges(6,2) * t66 / 0.2e1 - Ifges(6,4) * t64 + Ifges(6,1) * t120;
t107 = m(6) / 0.2e1;
t106 = -t41 / 0.2e1;
t104 = pkin(4) * t52;
t92 = t39 * t52;
t3 = m(6) * (t117 * t67 - t98);
t88 = t3 * qJD(1);
t4 = m(6) * (0.1e1 - t89) * t98;
t87 = t4 * qJD(1);
t77 = Ifges(6,5) * t66 - Ifges(6,6) * t64;
t61 = Ifges(6,4) * t66;
t54 = -Ifges(6,2) * t64 + t61;
t55 = Ifges(6,1) * t64 + t61;
t76 = t113 * t64 + (t54 + t55) * t120;
t2 = (-t114 * t34 + t116 * t67) * t107 - t70 / 0.2e1;
t68 = -t41 * mrSges(5,2) + mrSges(6,3) * t78 + (t51 - mrSges(5,1)) * t42;
t6 = m(6) * (t39 * t42 + t40 * t78) + t68;
t75 = t2 * qJD(1) + t6 * qJD(2);
t11 = t96 / 0.2e1 - t72;
t15 = t76 + t92;
t74 = -t11 * qJD(1) + t15 * qJD(2);
t10 = (-t52 / 0.2e1 + t73) * t34;
t16 = t76 - t104;
t8 = (pkin(4) / 0.2e1 - t39 / 0.2e1) * t52 + (mrSges(6,2) * t106 - t55 / 0.2e1 - t54 / 0.2e1) * t66 + (mrSges(6,1) * t106 - t113) * t64;
t71 = -t10 * qJD(1) + t8 * qJD(2) - t16 * qJD(4);
t9 = -t104 / 0.2e1 + t92 / 0.2e1 - t73 * t41 + t76;
t1 = t70 / 0.2e1 + t111 / 0.2e1 - t95 / 0.2e1 - t99 / 0.2e1 + (t116 * t107 + t51 / 0.2e1 - mrSges(5,1) / 0.2e1) * t67 - (mrSges(5,2) / 0.2e1 + t114 * t107 + (-t63 / 0.2e1 - t62 / 0.2e1) * mrSges(6,3)) * t34 + t108 / 0.2e1;
t5 = [qJD(2) * t3 - qJD(4) * t4, t88 + (m(6) * (t117 * t40 + t39 * t67) + m(5) * (t34 * t42 - t41 * t67) + m(4) * (t86 * t49 + t85 * t69) * pkin(2) - t103 * mrSges(3,2) - t101 * mrSges(3,1) - t69 * mrSges(4,2) + t49 * mrSges(4,1) + t122) * qJD(2) + t1 * qJD(4) + t121, 0, -t87 + t1 * qJD(2) + (t70 + t122) * qJD(4) + t121, qJD(5) * t111 + (qJD(2) + qJD(4)) * t119; qJD(4) * t2 - qJD(5) * t11 - t88, qJD(4) * t6 + qJD(5) * t15, 0, (m(6) * (-pkin(4) * t42 + pkin(7) * t78) + t68) * qJD(4) + t9 * qJD(5) + t75, t9 * qJD(4) + (t51 * t40 + t77) * qJD(5) + t74; 0, 0, 0, 0, -t52 * qJD(5); -qJD(2) * t2 + qJD(5) * t10 + t87, -qJD(5) * t8 - t75, 0, t16 * qJD(5), (t51 * pkin(7) + t77) * qJD(5) - t71; qJD(2) * t11 - qJD(4) * t10, qJD(4) * t8 - t74, 0, t71, 0;];
Cq = t5;
