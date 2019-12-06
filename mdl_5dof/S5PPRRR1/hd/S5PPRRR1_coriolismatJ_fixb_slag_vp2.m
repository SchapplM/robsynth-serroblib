% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:40
% EndTime: 2019-12-05 15:12:41
% DurationCPUTime: 0.63s
% Computational Cost: add. (1201->99), mult. (2779->148), div. (0->0), fcn. (2977->8), ass. (0->68)
t81 = sin(pkin(9));
t82 = cos(pkin(9));
t97 = sin(qJ(3));
t98 = cos(qJ(3));
t43 = t98 * t81 + t97 * t82;
t60 = sin(qJ(4));
t62 = cos(qJ(4));
t65 = -t97 * t81 + t98 * t82;
t33 = t60 * t43 - t62 * t65;
t61 = cos(qJ(5));
t86 = t61 * mrSges(6,2);
t59 = sin(qJ(5));
t87 = t59 * mrSges(6,1);
t69 = t86 / 0.2e1 + t87 / 0.2e1;
t68 = t69 * t33;
t46 = t86 + t87;
t92 = t33 * t46;
t117 = t92 / 0.2e1 + t68;
t120 = t117 * qJD(5);
t119 = t61 / 0.2e1;
t57 = t59 ^ 2;
t58 = t61 ^ 2;
t85 = t57 + t58;
t118 = t85 * t62;
t114 = t85 * t33;
t100 = t60 * pkin(3);
t53 = pkin(7) + t100;
t99 = t62 * pkin(3);
t54 = -pkin(4) - t99;
t63 = t62 * t43 + t60 * t65;
t116 = -t53 * t114 + t54 * t63;
t66 = m(6) * (-pkin(4) * t63 - pkin(7) * t114);
t115 = t33 * t60;
t93 = t33 * t63;
t45 = -t61 * mrSges(6,1) + t59 * mrSges(6,2);
t108 = t63 * t45;
t95 = t63 * mrSges(5,1);
t113 = t108 - t95;
t112 = -Ifges(6,2) * t61 / 0.2e1 - Ifges(6,4) * t59 + Ifges(6,1) * t119;
t106 = t85 * mrSges(6,3);
t104 = m(6) / 0.2e1;
t103 = -t33 / 0.2e1;
t101 = pkin(4) * t46;
t90 = t54 * t46;
t3 = m(6) * (-t114 * t63 + t93);
t84 = t3 * qJD(1);
t4 = m(6) * (0.1e1 - t85) * t93;
t83 = t4 * qJD(1);
t80 = -t99 / 0.2e1;
t73 = Ifges(6,5) * t61 - Ifges(6,6) * t59;
t56 = Ifges(6,4) * t61;
t48 = -Ifges(6,2) * t59 + t56;
t49 = Ifges(6,1) * t59 + t56;
t72 = t112 * t59 + (t48 + t49) * t119;
t64 = -t62 * mrSges(5,2) + mrSges(6,3) * t118 + (t45 - mrSges(5,1)) * t60;
t14 = (m(6) * (t118 * t53 + t54 * t60) + t64) * pkin(3);
t2 = (t33 / 0.2e1 + t103) * mrSges(5,2) + ((t118 * t63 + t115) * pkin(3) + t116) * t104 - t66 / 0.2e1;
t71 = t2 * qJD(1) + t14 * qJD(3);
t15 = t72 + t90;
t7 = -t92 / 0.2e1 + t68;
t70 = -t7 * qJD(1) + t15 * qJD(3);
t11 = (-t54 / 0.2e1 + pkin(4) / 0.2e1) * t46 + (mrSges(6,2) * t80 - t49 / 0.2e1 - t48 / 0.2e1) * t61 + (mrSges(6,1) * t80 - t112) * t59;
t16 = t72 - t101;
t6 = (-t46 / 0.2e1 + t69) * t33;
t67 = t6 * qJD(1) + t11 * qJD(3) - t16 * qJD(4);
t12 = t90 / 0.2e1 - t101 / 0.2e1 - t69 * t99 + t72;
t1 = t66 / 0.2e1 + t108 / 0.2e1 - t95 / 0.2e1 + ((pkin(3) * t118 + t54) * t104 + t45 / 0.2e1 - mrSges(5,1) / 0.2e1) * t63 + t103 * t106 + (mrSges(5,2) + (-t85 * t53 + t100) * t104 + (-t58 / 0.2e1 - t57 / 0.2e1) * mrSges(6,3)) * t33;
t5 = [qJD(3) * t3 + qJD(4) * t4, 0, t1 * qJD(4) + t120 + t84 + (-t43 * mrSges(4,1) - t65 * mrSges(4,2) - (-mrSges(5,2) + t106) * t33 + 0.2e1 * t116 * t104 + m(5) * (-t62 * t63 - t115) * pkin(3) + t113) * qJD(3), t83 + t1 * qJD(3) + (t33 * mrSges(5,2) - mrSges(6,3) * t114 + t113 + t66) * qJD(4) + t120, qJD(5) * t108 + (qJD(3) + qJD(4)) * t117; 0, 0, 0, 0, -t46 * qJD(5); qJD(4) * t2 - qJD(5) * t7 - t84, 0, qJD(4) * t14 + qJD(5) * t15, t12 * qJD(5) + (m(6) * (-pkin(4) * t60 + pkin(7) * t118) + t64) * qJD(4) * pkin(3) + t71, t12 * qJD(4) + (t45 * t53 + t73) * qJD(5) + t70; -qJD(3) * t2 - qJD(5) * t6 - t83, 0, -qJD(5) * t11 - t71, t16 * qJD(5), (t45 * pkin(7) + t73) * qJD(5) - t67; qJD(3) * t7 + qJD(4) * t6, 0, qJD(4) * t11 - t70, t67, 0;];
Cq = t5;
