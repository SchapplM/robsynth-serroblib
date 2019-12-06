% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:52
% EndTime: 2019-12-05 16:39:54
% DurationCPUTime: 0.60s
% Computational Cost: add. (947->112), mult. (2049->132), div. (0->0), fcn. (1285->4), ass. (0->74)
t100 = Ifges(5,4) + Ifges(6,4);
t64 = sin(qJ(4));
t137 = t100 * t64;
t115 = t64 * pkin(4);
t66 = cos(qJ(4));
t77 = -t66 * mrSges(6,1) + t64 * mrSges(6,2);
t42 = t77 * t115;
t63 = t66 ^ 2;
t136 = -t100 * t63 - t42;
t67 = cos(qJ(3));
t112 = t67 * pkin(2);
t135 = -t112 / 0.2e1;
t134 = Ifges(5,1) + Ifges(6,1);
t133 = Ifges(5,2) + Ifges(6,2);
t113 = t66 * pkin(4);
t85 = -pkin(3) - t113;
t95 = t64 ^ 2 + t63;
t55 = m(6) * t115;
t96 = t64 * mrSges(6,1) + t66 * mrSges(6,2);
t31 = -t55 - t96;
t131 = (qJD(2) + qJD(3)) * t31;
t120 = m(6) * pkin(4);
t129 = -mrSges(6,1) - t120;
t65 = sin(qJ(3));
t114 = t65 * pkin(2);
t54 = pkin(7) + t114;
t93 = qJ(5) + t54;
t33 = t93 * t66;
t22 = t33 * t66;
t32 = t93 * t64;
t128 = t32 * t64 + t22;
t98 = -qJ(5) - pkin(7);
t48 = t98 * t66;
t43 = t48 * t66;
t46 = t98 * t64;
t127 = -t46 * t64 - t43;
t126 = t133 * t66;
t124 = t134 * t66;
t122 = m(5) * pkin(2);
t121 = m(6) * pkin(2);
t116 = pkin(3) * t65;
t101 = t66 * mrSges(5,1);
t99 = -Ifges(5,6) - Ifges(6,6);
t97 = t95 * mrSges(6,3);
t69 = (t64 * mrSges(5,2) - mrSges(4,1) - t101 + t77) * t114 + (-mrSges(4,2) + t95 * (mrSges(5,3) + mrSges(6,3))) * t112;
t79 = t85 * t65;
t5 = (t79 + (-t114 + t128) * t67) * t121 + (-t116 + (t95 * t54 - t114) * t67) * t122 + t69;
t94 = t5 * qJD(2);
t92 = qJD(4) * t64;
t12 = m(6) * t128 + t97;
t91 = t12 * qJD(2);
t87 = t114 / 0.2e1;
t82 = (t32 - t46) * t64;
t78 = t64 * mrSges(5,1) + t66 * mrSges(5,2);
t76 = -mrSges(6,3) * t113 + (Ifges(5,5) + Ifges(6,5)) * t66;
t14 = m(6) * t127 + t97;
t7 = (t87 - t22 / 0.2e1 + t43 / 0.2e1 - t82 / 0.2e1) * m(6) - t97;
t75 = -t7 * qJD(2) + t14 * qJD(3);
t74 = t85 - t112;
t21 = t74 * t96;
t27 = (-pkin(3) - t112) * t78;
t70 = -t124 + t126 + t137;
t3 = -t21 - t27 + (-t74 * t120 + t70) * t64 + t136;
t72 = t3 * qJD(2);
t28 = t85 * t96;
t68 = t21 / 0.2e1 + t27 / 0.2e1 + t28 / 0.2e1 - pkin(3) * t78 / 0.2e1 + (t135 + t85) * t55 + (-t137 - t126 / 0.2e1 + t124 / 0.2e1 + (-t133 + t134) * t66 / 0.2e1) * t64 - t136;
t1 = -t68 + ((mrSges(5,2) + mrSges(6,2)) * t66 + (mrSges(5,1) - t129) * t64) * t135;
t4 = -t28 - t42 + (pkin(3) * mrSges(5,2) - t100 * t66) * t66 + (pkin(3) * mrSges(5,1) - t85 * t120 + t70) * t64;
t71 = t1 * qJD(2) + t4 * qJD(3);
t26 = t31 * qJD(4);
t25 = t31 * qJD(5);
t8 = m(6) * (t22 - t43 + t82) / 0.2e1 + m(6) * t87 + t97;
t2 = ((-mrSges(5,2) / 0.2e1 - mrSges(6,2) / 0.2e1) * t66 + (-mrSges(5,1) / 0.2e1 - t120 / 0.2e1 - mrSges(6,1) / 0.2e1) * t64) * t112 + t68;
t6 = [0, 0, 0, (-t78 + t31) * qJD(4), 0; 0, t5 * qJD(3) - t3 * qJD(4) + t12 * qJD(5), t94 + ((t127 * t67 + t79) * t121 + (t95 * t67 * pkin(7) - t116) * t122 + t69) * qJD(3) + t2 * qJD(4) + t8 * qJD(5), t2 * qJD(3) + (t32 * mrSges(6,2) - t54 * t101 + t129 * t33 + t76) * qJD(4) + (mrSges(5,2) * t54 + t99) * t92 - t72, t8 * qJD(3) + t91; 0, -t1 * qJD(4) - t7 * qJD(5) - t94, -t4 * qJD(4) + t14 * qJD(5), (-t46 * mrSges(6,2) - pkin(7) * t101 - t129 * t48 + t76) * qJD(4) + (mrSges(5,2) * pkin(7) + t99) * t92 - t71, t75; 0, t1 * qJD(3) + t25 + t72, t25 + t71, 0, t131; 0, t7 * qJD(3) - t26 - t91, -t26 - t75, -t131, 0;];
Cq = t6;
