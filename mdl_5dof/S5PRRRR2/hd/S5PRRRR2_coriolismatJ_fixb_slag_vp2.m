% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:45
% EndTime: 2019-12-05 17:04:47
% DurationCPUTime: 0.70s
% Computational Cost: add. (1096->111), mult. (2793->153), div. (0->0), fcn. (1846->6), ass. (0->82)
t118 = m(6) / 0.2e1;
t57 = sin(qJ(4));
t108 = cos(qJ(3));
t78 = t108 * pkin(2) + pkin(3);
t58 = sin(qJ(3));
t60 = cos(qJ(4));
t96 = t58 * t60;
t35 = pkin(2) * t96 + t57 * t78;
t31 = pkin(6) + t35;
t56 = sin(qJ(5));
t54 = t56 ^ 2;
t59 = cos(qJ(5));
t55 = t59 ^ 2;
t93 = t55 + t54;
t128 = t31 * t93;
t71 = -t35 + t128;
t132 = -mrSges(5,2) / 0.2e1 + t71 * t118;
t131 = -Ifges(6,1) + Ifges(6,2);
t120 = t93 * mrSges(6,3);
t76 = -mrSges(5,2) + t120;
t85 = m(6) * t93;
t79 = pkin(6) * t85;
t130 = t79 + t76;
t44 = -mrSges(6,1) * t59 + mrSges(6,2) * t56;
t94 = t44 - mrSges(5,1);
t98 = t57 * t58;
t34 = pkin(2) * t98 - t60 * t78;
t87 = t44 / 0.2e1 - mrSges(5,1) / 0.2e1;
t126 = (t34 * t118 + t87) * t57;
t125 = (-t55 / 0.2e1 - t54 / 0.2e1) * mrSges(6,3);
t107 = Ifges(6,4) * t56;
t121 = t131 * t59;
t62 = -Ifges(6,4) * t55 + (t107 + t121) * t56;
t123 = qJD(4) * t62;
t122 = qJD(5) * t62;
t39 = (t108 * t60 - t98) * pkin(2);
t119 = t76 * t39 + (-t58 * mrSges(4,1) - t108 * mrSges(4,2)) * pkin(2);
t115 = mrSges(5,2) / 0.2e1;
t114 = -t39 / 0.2e1;
t111 = t59 / 0.2e1;
t110 = pkin(3) * t60;
t109 = t57 * pkin(3);
t38 = (t108 * t57 + t96) * pkin(2);
t104 = t34 * t38;
t103 = t38 * mrSges(5,1);
t102 = t38 * t44;
t101 = t38 * t60;
t100 = t56 * mrSges(6,1);
t95 = t59 * mrSges(6,2);
t83 = t94 * t35;
t3 = t83 + (-m(6) * t71 - t76) * t34;
t92 = t3 * qJD(2);
t4 = t94 * t38 + m(6) * (t39 * t128 + t104) + m(5) * (t35 * t39 + t104) + t119;
t91 = t4 * qJD(2);
t15 = 0.2e1 * Ifges(6,4) * t59 * t111 + (-t107 + t131 * (-t59 / 0.2e1 - t111)) * t56;
t74 = t95 + t100;
t68 = t34 * t74;
t12 = t68 + t15;
t90 = t12 * qJD(2);
t82 = t94 * t57;
t52 = pkin(6) + t109;
t81 = t93 * t52;
t80 = Ifges(6,5) * t59 - Ifges(6,6) * t56;
t65 = t74 * t110;
t77 = -t65 / 0.2e1 + t15;
t13 = pkin(3) * t82 + (t76 + (-t109 + t81) * m(6)) * t110;
t66 = t115 - t52 * t85 / 0.2e1;
t2 = (t115 - t79 / 0.2e1) * t39 + t66 * t34 + (t132 * t60 + t126) * pkin(3) + (t110 / 0.2e1 - t34 / 0.2e1 + t114) * t120 + t94 * (t35 / 0.2e1 - t38 / 0.2e1);
t73 = t2 * qJD(2) + t13 * qJD(3);
t61 = (-t55 + t54) * Ifges(6,4) + t56 * t121;
t14 = t65 + t61;
t69 = t95 / 0.2e1 + t100 / 0.2e1;
t64 = t69 * t39;
t6 = t65 / 0.2e1 - t68 / 0.2e1 - t64 + t61;
t72 = -t6 * qJD(2) - t14 * qJD(3);
t63 = -t15 * qJD(4) + (qJD(2) + qJD(3)) * t62;
t18 = t68 / 0.2e1;
t11 = -t69 * t110 + t77;
t9 = t69 * t34 + t15 + t18;
t7 = t18 - t64 + t77;
t1 = t102 / 0.2e1 + mrSges(5,2) * t114 - t103 / 0.2e1 + t87 * t35 + (t66 + t125) * t34 + (t126 + (-t125 + t132) * t60) * pkin(3) + (t79 + t120) * t39 / 0.2e1;
t5 = [0, 0, 0, 0, -t74 * qJD(5); 0, qJD(3) * t4 + qJD(4) * t3 + qJD(5) * t12, t1 * qJD(4) + t7 * qJD(5) + t91 + (t102 - t103 + 0.2e1 * (-pkin(3) * t101 + t39 * t81) * t118 + m(5) * (t39 * t57 - t101) * pkin(3) + t119) * qJD(3), t1 * qJD(3) + t9 * qJD(5) + t92 + (-t130 * t34 + t83) * qJD(4), t90 + t7 * qJD(3) + t9 * qJD(4) + (t44 * t31 + t80) * qJD(5); 0, qJD(4) * t2 - qJD(5) * t6 - t91, qJD(4) * t13 - qJD(5) * t14, t11 * qJD(5) + (t130 * t60 + t82) * qJD(4) * pkin(3) + t73, t11 * qJD(4) + (t44 * t52 + t80) * qJD(5) + t72; 0, -qJD(3) * t2 - t122 - t92, -t73 - t122, t15 * qJD(5), (t44 * pkin(6) + t80) * qJD(5) - t63; 0, qJD(3) * t6 + t123 - t90, -t72 + t123, t63, 0;];
Cq = t5;
