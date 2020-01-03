% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:54
% EndTime: 2019-12-31 17:12:55
% DurationCPUTime: 0.50s
% Computational Cost: add. (888->112), mult. (1898->135), div. (0->0), fcn. (1196->4), ass. (0->73)
t55 = sin(qJ(3));
t91 = Ifges(4,4) + Ifges(5,4);
t124 = t91 * t55;
t123 = Ifges(4,1) + Ifges(5,1);
t122 = Ifges(4,2) + Ifges(5,2);
t57 = cos(qJ(3));
t121 = t91 * t57;
t54 = t57 ^ 2;
t87 = t55 ^ 2 + t54;
t110 = m(5) * pkin(3);
t115 = mrSges(5,1) + t110;
t92 = t57 * mrSges(5,2);
t26 = t115 * t55 + t92;
t119 = (qJD(1) + qJD(2)) * t26;
t117 = t122 * t57;
t114 = t123 * t57;
t113 = m(5) / 0.2e1;
t112 = m(4) * pkin(1);
t111 = m(5) * pkin(1);
t56 = sin(qJ(2));
t106 = pkin(2) * t56;
t105 = t55 * pkin(3);
t104 = t56 * pkin(1);
t103 = t57 * pkin(3);
t58 = cos(qJ(2));
t102 = t58 * pkin(1);
t93 = t57 * mrSges(4,1);
t90 = -Ifges(4,6) - Ifges(5,6);
t89 = -qJ(4) - pkin(6);
t88 = t87 * mrSges(5,3);
t76 = -pkin(2) - t103;
t62 = t76 - t102;
t68 = t55 * mrSges(5,1) + t92;
t17 = t62 * t68;
t70 = t55 * mrSges(4,1) + t57 * mrSges(4,2);
t22 = (-pkin(2) - t102) * t70;
t69 = -t57 * mrSges(5,1) + t55 * mrSges(5,2);
t37 = t69 * t105;
t61 = -t114 + t117 + t124;
t3 = -t17 - t22 - t37 - t91 * t54 + (-t62 * t110 + t61) * t55;
t86 = t3 * qJD(1);
t48 = pkin(6) + t104;
t60 = (t55 * mrSges(4,2) - mrSges(3,1) + t69 - t93) * t104 + (-mrSges(3,2) + t87 * (mrSges(4,3) + mrSges(5,3))) * t102;
t84 = qJ(4) + t48;
t27 = t84 * t55;
t28 = t84 * t57;
t67 = t27 * t55 + t28 * t57;
t71 = t76 * t56;
t5 = (t71 + (t67 - t104) * t58) * t111 + (-t106 + (t87 * t48 - t104) * t58) * t112 + t60;
t85 = t5 * qJD(1);
t83 = qJD(3) * t55;
t12 = m(5) * t67 + t88;
t82 = t12 * qJD(1);
t78 = t104 / 0.2e1;
t40 = t89 * t55;
t75 = (t27 - t40) * t55;
t42 = t89 * t57;
t74 = (t28 - t42) * t57;
t66 = -t40 * t55 - t42 * t57;
t65 = -mrSges(5,3) * t103 + (Ifges(4,5) + Ifges(5,5)) * t57;
t23 = t76 * t68;
t59 = t37 + t17 / 0.2e1 + t22 / 0.2e1 + t23 / 0.2e1 - pkin(2) * t70 / 0.2e1 + (-0.2e1 * pkin(2) - t102 - 0.2e1 * t103) * t105 * t113 + t121 * t57 + (-t124 - t117 / 0.2e1 + t114 / 0.2e1 + (-t122 + t123) * t57 / 0.2e1) * t55;
t1 = -t59 - ((mrSges(4,2) + mrSges(5,2)) * t57 + (mrSges(4,1) + t115) * t55) * t102 / 0.2e1;
t4 = -t23 - t37 + (pkin(2) * mrSges(4,2) - t121) * t57 + (pkin(2) * mrSges(4,1) - t76 * t110 + t61) * t55;
t64 = -t1 * qJD(1) - t4 * qJD(2);
t13 = m(5) * t66 + t88;
t7 = (t78 - t74 / 0.2e1 - t75 / 0.2e1) * m(5) - t88;
t63 = -t7 * qJD(1) + t13 * qJD(2);
t21 = t26 * qJD(4);
t20 = t26 * qJD(3);
t8 = (t74 + t75) * t113 + m(5) * t78 + t88;
t2 = ((-mrSges(4,2) / 0.2e1 - mrSges(5,2) / 0.2e1) * t57 + (-mrSges(4,1) / 0.2e1 - t110 / 0.2e1 - mrSges(5,1) / 0.2e1) * t55) * t102 + t59;
t6 = [t5 * qJD(2) - t3 * qJD(3) + t12 * qJD(4), t85 + ((t66 * t58 + t71) * t111 + (t87 * t58 * pkin(6) - t106) * t112 + t60) * qJD(2) + t2 * qJD(3) + t8 * qJD(4), -t86 + t2 * qJD(2) + (t27 * mrSges(5,2) - t115 * t28 - t48 * t93 + t65) * qJD(3) + (mrSges(4,2) * t48 + t90) * t83, t8 * qJD(2) + t82; -t1 * qJD(3) - t7 * qJD(4) - t85, -t4 * qJD(3) + t13 * qJD(4), (-t40 * mrSges(5,2) - pkin(6) * t93 + t115 * t42 + t65) * qJD(3) + (mrSges(4,2) * pkin(6) + t90) * t83 + t64, t63; t1 * qJD(2) - t21 + t86, -t21 - t64, 0, -t119; t7 * qJD(2) + t20 - t82, t20 - t63, t119, 0;];
Cq = t6;
