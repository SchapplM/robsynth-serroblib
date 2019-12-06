% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRP2
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:40
% EndTime: 2019-12-05 16:41:41
% DurationCPUTime: 0.52s
% Computational Cost: add. (691->111), mult. (1449->138), div. (0->0), fcn. (904->4), ass. (0->66)
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t34 = -t56 * mrSges(6,1) - t54 * mrSges(6,3);
t36 = -t54 * pkin(4) + t56 * qJ(5);
t17 = t36 * t34;
t53 = t56 ^ 2;
t81 = Ifges(6,5) - Ifges(5,4);
t108 = -t81 * t53 - t17;
t107 = t81 * t54;
t106 = t54 ^ 2 + t53;
t105 = Ifges(5,1) + Ifges(6,1);
t92 = cos(qJ(3)) * pkin(2);
t104 = t106 * t92;
t103 = -Ifges(6,3) + t105;
t102 = -t56 * mrSges(5,1) + t54 * mrSges(5,2) + t34;
t93 = sin(qJ(3)) * pkin(2);
t101 = (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t106) * t92 + (-mrSges(4,1) + t102) * t93;
t100 = -mrSges(6,1) / 0.2e1;
t75 = t54 * qJ(5);
t66 = -t56 * pkin(4) - t75;
t32 = -pkin(3) + t66;
t20 = t32 - t92;
t99 = t20 / 0.2e1;
t96 = m(6) * t36;
t95 = m(6) * t56;
t83 = t56 * mrSges(5,2);
t85 = t54 * mrSges(5,1);
t37 = t83 + t85;
t94 = pkin(3) * t37;
t91 = m(6) * qJ(5);
t46 = -pkin(3) - t92;
t86 = t46 * t37;
t84 = t54 * mrSges(6,1);
t48 = t56 * mrSges(6,2);
t82 = t56 * mrSges(6,3);
t45 = pkin(7) + t93;
t80 = t104 * t45;
t79 = t104 * pkin(7);
t62 = ((-Ifges(5,2) - Ifges(6,3) / 0.2e1 + t103 / 0.2e1 + t105 / 0.2e1) * t56 + t107) * t54 + t108;
t67 = -t82 + t84;
t3 = t86 + (t67 - t96) * t20 + t62;
t77 = t3 * qJD(2);
t5 = m(6) * (t20 * t93 + t80) + m(5) * (t46 * t93 + t80) + t101;
t76 = t5 * qJD(2);
t10 = (m(6) * t20 + t34) * t54;
t74 = t10 * qJD(2);
t47 = mrSges(6,3) + t91;
t73 = t47 * qJD(4);
t70 = t92 / 0.2e1;
t69 = m(6) * (-t20 - t32);
t68 = t69 / 0.2e1;
t16 = t32 * t67;
t61 = -t107 + (Ifges(5,2) - t103) * t56;
t1 = t17 - t16 / 0.2e1 + (-t46 / 0.2e1 + pkin(3) / 0.2e1) * t37 - t36 * t69 / 0.2e1 + (mrSges(6,3) * t99 + t81 * t56 + (-mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1 + t91 / 0.2e1) * t92) * t56 + (t20 * t100 + (-mrSges(5,1) / 0.2e1 + t100 - m(6) * pkin(4) / 0.2e1) * t92 + t61) * t54;
t4 = t32 * t96 + t61 * t54 - t108 - t16 + t94;
t64 = t1 * qJD(2) + t4 * qJD(3);
t11 = (m(6) * t32 + t34) * t54;
t6 = (t34 + (t70 + t99 + t32 / 0.2e1) * m(6)) * t54;
t63 = t6 * qJD(2) + t11 * qJD(3);
t59 = (-mrSges(6,2) * t75 - pkin(4) * t48 + (Ifges(6,4) + Ifges(5,5)) * t56 + (-Ifges(5,6) + Ifges(6,6)) * t54) * qJD(4);
t58 = qJD(4) * (m(6) * t66 + t102);
t33 = pkin(7) * t95 + t48;
t19 = t45 * t95 + t48;
t7 = (m(6) * t70 - t34 + t68) * t54;
t2 = t62 + t36 * t68 + (-t83 / 0.2e1 - t85 / 0.2e1 - t84 / 0.2e1 + t82 / 0.2e1 + t96 / 0.2e1) * t92 + t16 / 0.2e1 + t67 * t99 + t86 / 0.2e1 - t94 / 0.2e1;
t8 = [0, 0, 0, (t96 + (-mrSges(5,2) + mrSges(6,3)) * t56) * qJD(4) + ((-mrSges(5,1) - mrSges(6,1)) * qJD(4) + m(6) * qJD(5)) * t54, m(6) * qJD(4) * t54; 0, t5 * qJD(3) + t3 * qJD(4) - t10 * qJD(5), t2 * qJD(4) + t7 * qJD(5) + t76 + (m(6) * (t32 * t93 + t79) + m(5) * (-pkin(3) * t93 + t79) + t101) * qJD(3), t2 * qJD(3) + t19 * qJD(5) + t45 * t58 + t59 + t77, t7 * qJD(3) + t19 * qJD(4) - t74; 0, -t1 * qJD(4) - t6 * qJD(5) - t76, -t4 * qJD(4) - t11 * qJD(5), pkin(7) * t58 + t33 * qJD(5) + t59 - t64, t33 * qJD(4) - t63; 0, t1 * qJD(3) - t77, t64, t47 * qJD(5), t73; 0, t6 * qJD(3) + t74, t63, -t73, 0;];
Cq = t8;
