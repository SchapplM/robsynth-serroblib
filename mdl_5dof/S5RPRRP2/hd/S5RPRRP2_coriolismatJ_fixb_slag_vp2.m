% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP2
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:10
% EndTime: 2020-01-03 11:45:12
% DurationCPUTime: 0.71s
% Computational Cost: add. (1429->115), mult. (2805->134), div. (0->0), fcn. (2011->6), ass. (0->76)
t102 = Ifges(5,4) + Ifges(6,4);
t69 = sin(qJ(4));
t125 = t102 * t69;
t133 = Ifges(5,2) + Ifges(6,2);
t71 = cos(qJ(4));
t136 = t133 * t71 + t125;
t135 = t69 / 0.2e1;
t134 = Ifges(5,1) + Ifges(6,1);
t115 = t71 * pkin(4);
t89 = -pkin(3) - t115;
t132 = t102 * t71;
t116 = t69 * pkin(4);
t59 = m(6) * t116;
t98 = t69 * mrSges(6,1) + t71 * mrSges(6,2);
t42 = -t59 - t98;
t130 = (qJD(1) + qJD(3)) * t42;
t121 = m(6) * pkin(4);
t129 = -mrSges(6,1) - t121;
t117 = pkin(1) * sin(pkin(8));
t58 = cos(pkin(8)) * pkin(1) + pkin(2);
t70 = sin(qJ(3));
t72 = cos(qJ(3));
t35 = t72 * t117 + t70 * t58;
t33 = pkin(7) + t35;
t96 = qJ(5) + t33;
t25 = t96 * t71;
t21 = t25 * t71;
t24 = t96 * t69;
t128 = t24 * t69 + t21;
t100 = -qJ(5) - pkin(7);
t56 = t100 * t71;
t52 = t56 * t71;
t54 = t100 * t69;
t127 = -t54 * t69 - t52;
t124 = t134 * t71;
t67 = t71 ^ 2;
t120 = -t69 / 0.2e1;
t62 = t67 * mrSges(6,3);
t103 = t71 * mrSges(5,1);
t101 = -Ifges(5,6) - Ifges(6,6);
t66 = t69 ^ 2;
t99 = t66 * mrSges(6,3) + t62;
t34 = -t70 * t117 + t72 * t58;
t83 = -t71 * mrSges(6,1) + t69 * mrSges(6,2);
t74 = (t69 * mrSges(5,2) - mrSges(4,1) - t103 + t83) * t35 + ((mrSges(5,3) + mrSges(6,3)) * t66 + t67 * mrSges(5,3) - mrSges(4,2) + t62) * t34;
t85 = -pkin(3) - t34;
t80 = t85 - t115;
t86 = (t66 + t67) * t34;
t3 = m(6) * (t128 * t34 + t80 * t35) + m(5) * (t33 * t86 + t85 * t35) + t74;
t97 = t3 * qJD(1);
t95 = qJD(4) * t69;
t11 = m(6) * t128 + t99;
t94 = t11 * qJD(1);
t92 = m(6) * t35 / 0.2e1;
t84 = t69 * mrSges(5,1) + t71 * mrSges(5,2);
t82 = -mrSges(6,3) * t115 + (Ifges(5,5) + Ifges(6,5)) * t71;
t15 = m(6) * t127 + t99;
t76 = -(t21 - t52 + (t24 - t54) * t69) * m(6) / 0.2e1 - t99;
t8 = t92 + t76;
t81 = -t8 * qJD(1) + t15 * qJD(3);
t20 = t80 * t98;
t23 = t85 * t84;
t51 = t83 * t116;
t75 = -t124 + t136;
t4 = -t20 - t23 - t51 - t102 * t67 + (-t80 * t121 + t75) * t69;
t78 = t4 * qJD(1);
t39 = t89 * t98;
t73 = t51 + t20 / 0.2e1 + t23 / 0.2e1 + t39 / 0.2e1 + t136 * t120 + (t124 - t125) * t135 - pkin(3) * t84 / 0.2e1 + ((-t133 + t134) * t135 + t132) * t71 + (-t34 / 0.2e1 + t89) * t59;
t1 = -t73 + (-(mrSges(5,2) + mrSges(6,2)) * t71 / 0.2e1 + (mrSges(5,1) - t129) * t120) * t34;
t5 = -t39 - t51 + (pkin(3) * mrSges(5,2) - t132) * t71 + (pkin(3) * mrSges(5,1) - t89 * t121 + t75) * t69;
t77 = t1 * qJD(1) + t5 * qJD(3);
t38 = t42 * qJD(4);
t37 = t42 * qJD(5);
t9 = t92 - t76;
t2 = ((-mrSges(5,2) / 0.2e1 - mrSges(6,2) / 0.2e1) * t71 + (-mrSges(5,1) / 0.2e1 - t121 / 0.2e1 - mrSges(6,1) / 0.2e1) * t69) * t34 + t73;
t6 = [t3 * qJD(3) - t4 * qJD(4) + t11 * qJD(5), 0, t97 + (m(6) * (t127 * t34 + t89 * t35) + m(5) * (-pkin(3) * t35 + pkin(7) * t86) + t74) * qJD(3) + t2 * qJD(4) + t9 * qJD(5), t2 * qJD(3) + (t24 * mrSges(6,2) - t33 * t103 + t129 * t25 + t82) * qJD(4) + (mrSges(5,2) * t33 + t101) * t95 - t78, t9 * qJD(3) + t94; 0, 0, 0, (-t84 + t42) * qJD(4), 0; -t1 * qJD(4) - t8 * qJD(5) - t97, 0, -t5 * qJD(4) + t15 * qJD(5), (-t54 * mrSges(6,2) - pkin(7) * t103 - t129 * t56 + t82) * qJD(4) + (mrSges(5,2) * pkin(7) + t101) * t95 - t77, t81; t1 * qJD(3) + t37 + t78, 0, t37 + t77, 0, t130; t8 * qJD(3) - t38 - t94, 0, -t38 - t81, -t130, 0;];
Cq = t6;
