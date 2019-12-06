% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:47
% EndTime: 2019-12-05 15:59:49
% DurationCPUTime: 0.77s
% Computational Cost: add. (2182->112), mult. (4534->159), div. (0->0), fcn. (4339->6), ass. (0->70)
t95 = sin(qJ(5));
t96 = sin(qJ(4));
t98 = cos(qJ(5));
t99 = cos(qJ(4));
t112 = t95 * t99 + t98 * t96;
t101 = -pkin(2) - pkin(6);
t124 = t96 * t101;
t82 = -t96 * pkin(7) + t124;
t83 = (-pkin(7) + t101) * t99;
t115 = -t95 * t82 + t98 * t83;
t56 = t98 * t82 + t95 * t83;
t79 = t95 * t96 - t98 * t99;
t7 = -t56 * mrSges(6,1) - t115 * mrSges(6,2) - Ifges(6,5) * t112 + Ifges(6,6) * t79;
t164 = t7 * qJD(5);
t135 = t99 * mrSges(5,2);
t114 = t96 * mrSges(5,1) + t135;
t116 = -mrSges(6,1) * t112 + t79 * mrSges(6,2);
t104 = -t114 + t116;
t163 = mrSges(4,3) - t104;
t162 = m(6) * pkin(4);
t149 = pkin(4) * t99;
t161 = m(6) * t149;
t158 = t98 * t112 + t95 * t79;
t157 = t116 * qJD(5);
t100 = cos(qJ(2));
t67 = t79 * t100;
t68 = t112 * t100;
t117 = t68 * mrSges(6,1) - t67 * mrSges(6,2);
t156 = t117 * qJD(5);
t52 = -mrSges(6,1) * t79 - mrSges(6,2) * t112;
t86 = t96 * pkin(4) + qJ(3);
t4 = (-Ifges(6,4) * t79 - (-Ifges(6,1) + Ifges(6,2)) * t112) * t79 + Ifges(6,4) * t112 ^ 2 + t86 * t52;
t97 = sin(qJ(2));
t66 = t79 * t97;
t69 = t112 * t97;
t127 = -t66 * mrSges(6,1) / 0.2e1 - t69 * mrSges(6,2) / 0.2e1;
t94 = t99 ^ 2;
t153 = m(5) / 0.2e1;
t152 = -m(6) / 0.2e1;
t151 = m(6) / 0.2e1;
t150 = t97 / 0.2e1;
t139 = t96 * mrSges(5,2);
t138 = t97 * t52;
t136 = t99 * mrSges(5,1);
t93 = t96 ^ 2;
t125 = t93 + t94;
t87 = t97 * t100;
t24 = m(5) * (-t125 + 0.1e1) * t87 + m(6) * (-t67 * t66 - t68 * t69 + t87);
t122 = t24 * qJD(1);
t120 = t138 / 0.2e1;
t113 = -t112 * t69 - t66 * t79;
t5 = t120 - t127;
t111 = t5 * qJD(1) + t4 * qJD(2);
t81 = (mrSges(6,1) * t95 + mrSges(6,2) * t98) * pkin(4);
t110 = t81 * qJD(4);
t106 = t113 * t152;
t26 = t106 + (t152 + (t93 / 0.2e1 + t94 / 0.2e1 - 0.1e1 / 0.2e1) * m(5)) * t97;
t36 = m(6) * t86 + (m(5) + m(4)) * qJ(3) + t163;
t109 = -qJD(1) * t26 + qJD(2) * t36;
t1 = t86 * t161 - t116 * t149 + qJ(3) * t136 - t94 * Ifges(5,4) + (-qJ(3) * mrSges(5,2) + Ifges(5,4) * t96 + (-Ifges(5,1) + Ifges(5,2)) * t99) * t96 + t4;
t102 = t97 * t161;
t103 = (-t66 * t98 + t69 * t95) * pkin(4) * t151 + t127;
t2 = -t138 / 0.2e1 - t102 / 0.2e1 + t103;
t105 = -t2 * qJD(1) + t1 * qJD(2);
t91 = t100 * qJ(3);
t78 = t81 * qJD(5);
t25 = t106 + (t125 * t153 + m(4)) * t97 + (m(5) + m(6)) * t150;
t6 = t120 + t127;
t3 = (t136 - t139) * t150 + t102 / 0.2e1 + (-t139 / 0.2e1 + t136 / 0.2e1) * t97 + t103 + t120;
t8 = [t24 * qJD(2), t25 * qJD(3) + t3 * qJD(4) + t6 * qJD(5) + t122 + (t113 * mrSges(6,3) + (-t125 * mrSges(5,3) - mrSges(3,1) + mrSges(4,2)) * t97 + (-mrSges(3,2) + t163) * t100 + m(4) * (-t97 * pkin(2) + t91) + 0.2e1 * (t125 * t97 * t101 + t91) * t153 + 0.2e1 * (t86 * t100 - t115 * t66 + t56 * t69) * t151) * qJD(2), qJD(2) * t25, t3 * qJD(2) + (t100 * t114 + t117 + (t67 * t95 + t68 * t98) * t162) * qJD(4) + t156, t6 * qJD(2) + qJD(4) * t117 + t156; -qJD(3) * t26 - qJD(4) * t2 + qJD(5) * t5 - t122, qJD(3) * t36 + qJD(4) * t1 + qJD(5) * t4, t109, t164 + t105 + (-mrSges(5,1) * t124 - Ifges(5,5) * t96 - Ifges(5,6) * t99 - t101 * t135 + (m(6) * (t115 * t95 - t56 * t98) + t158 * mrSges(6,3)) * pkin(4) + t7) * qJD(4), t7 * qJD(4) + t111 + t164; qJD(2) * t26, -t109, 0, (-t158 * t162 + t104) * qJD(4) + t157, qJD(4) * t116 + t157; t2 * qJD(2), -t105, 0, -t78, -t110 - t78; -t5 * qJD(2), -t111, 0, t110, 0;];
Cq = t8;
