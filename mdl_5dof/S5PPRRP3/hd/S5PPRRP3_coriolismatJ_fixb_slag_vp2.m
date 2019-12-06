% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:36
% EndTime: 2019-12-05 15:10:37
% DurationCPUTime: 0.49s
% Computational Cost: add. (634->81), mult. (2026->128), div. (0->0), fcn. (1616->6), ass. (0->55)
t65 = sin(qJ(4));
t67 = cos(qJ(4));
t116 = t65 ^ 2 + t67 ^ 2;
t119 = pkin(6) * t116;
t87 = m(5) / 0.4e1 + m(6) / 0.4e1;
t101 = t67 * mrSges(6,3);
t104 = t65 * mrSges(6,1);
t79 = t65 * pkin(4) - t67 * qJ(5);
t74 = m(6) * t79;
t81 = -t101 + t104;
t118 = t74 + t65 * mrSges(5,1) + t104 / 0.2e1 + t67 * mrSges(5,2) - t101 / 0.2e1 + t81 / 0.2e1;
t52 = -t67 * mrSges(5,1) + t65 * mrSges(5,2);
t117 = t52 - mrSges(4,1);
t96 = m(6) * qJD(5);
t114 = (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t116) * qJD(3) + t65 * t96;
t113 = 0.2e1 * qJD(3);
t112 = m(5) / 0.2e1;
t110 = m(6) / 0.2e1;
t64 = sin(pkin(8));
t68 = cos(qJ(3));
t106 = t64 * t68;
t66 = sin(qJ(3));
t103 = t66 * t64;
t100 = t68 * t66;
t99 = Ifges(5,4) - Ifges(6,5);
t98 = t68 * t119;
t94 = cos(pkin(8));
t93 = qJD(3) * t66;
t92 = qJD(3) * t68;
t91 = qJD(4) * t67;
t80 = -t67 * pkin(4) - t65 * qJ(5);
t49 = -pkin(3) + t80;
t51 = -t67 * mrSges(6,1) - t65 * mrSges(6,3);
t84 = m(6) * t49 + t51;
t19 = t84 * t65;
t90 = t19 * qJD(3);
t61 = m(6) * qJ(5) + mrSges(6,3);
t89 = t61 * qJD(4);
t88 = -0.1e1 + t116;
t32 = t65 * t106 + t94 * t67;
t33 = t67 * t106 - t94 * t65;
t78 = t32 * t65 + t33 * t67;
t5 = 0.2e1 * (-t88 * t64 * t66 ^ 2 + (-t106 + t78) * t68) * t87;
t7 = (m(5) + m(6)) * (t64 ^ 2 * t100 - t78 * t103);
t77 = -t7 * qJD(1) - t5 * qJD(2);
t4 = -t79 * t51 + (pkin(3) * mrSges(5,2) - t99 * t67) * t67 + (pkin(3) * mrSges(5,1) + t99 * t65 + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,3)) * t67) * t65 + (-t74 - t81) * t49;
t76 = t4 * qJD(3);
t14 = 0.4e1 * t87 * t88 * t100;
t75 = t5 * qJD(1) + t14 * qJD(2);
t69 = (m(6) * t80 + t51 + t52) * qJD(4);
t50 = (m(6) * pkin(6) + mrSges(6,2)) * t67;
t10 = t118 * t68;
t3 = t5 * qJD(3);
t2 = t118 * t103;
t1 = [t7 * qJD(3), t3, t2 * qJD(4) + ((-m(5) * pkin(3) + t117 + t84) * t92 + (-(t110 + t112) * t113 * t119 - t114) * t66) * t64 - t77, t2 * qJD(3) + ((-mrSges(5,1) - mrSges(6,1)) * t33 + (mrSges(5,2) - mrSges(6,3)) * t32) * qJD(4) + 0.2e1 * ((-pkin(4) * t33 - t32 * qJ(5)) * qJD(4) / 0.2e1 + t33 * qJD(5) / 0.2e1) * m(6), (-t64 * t65 * t93 + t33 * qJD(4)) * m(6); t3, t14 * qJD(3), -t10 * qJD(4) + (t51 + t117) * t93 + ((t49 * t66 + t98) * t110 + (-pkin(3) * t66 + t98) * t112) * t113 + t114 * t68 + t75, -t10 * qJD(3) + (t67 * t96 + t69) * t66, (t65 * t92 + t66 * t91) * m(6); t77, -t75, -t4 * qJD(4) - t19 * qJD(5), t50 * qJD(5) + (-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t91 + (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * qJD(4) * t65 + pkin(6) * t69 - t76, t50 * qJD(4) - t90; 0, 0, t76, t61 * qJD(5), t89; 0, 0, t90, -t89, 0;];
Cq = t1;
