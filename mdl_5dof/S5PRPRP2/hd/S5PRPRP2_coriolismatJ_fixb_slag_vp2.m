% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:34
% EndTime: 2019-12-05 15:30:37
% DurationCPUTime: 0.82s
% Computational Cost: add. (1093->125), mult. (2695->183), div. (0->0), fcn. (2096->4), ass. (0->86)
t151 = -Ifges(5,4) - Ifges(6,4);
t70 = cos(qJ(4));
t137 = t70 ^ 2;
t69 = sin(qJ(4));
t138 = t69 ^ 2;
t150 = t137 + t138;
t67 = sin(pkin(8));
t68 = cos(pkin(8));
t53 = -pkin(3) * t68 - pkin(6) * t67 - pkin(2);
t94 = qJ(5) * t67 - t53;
t87 = t94 * t70;
t121 = qJ(3) * t69;
t92 = (pkin(4) + t121) * t68;
t149 = t87 + t92;
t148 = 0.1e1 - t150;
t147 = mrSges(5,1) + mrSges(6,1);
t122 = mrSges(5,2) + mrSges(6,2);
t146 = Ifges(5,1) + Ifges(6,1);
t145 = -Ifges(5,5) - Ifges(6,5);
t144 = -Ifges(5,2) - Ifges(6,2);
t143 = -Ifges(5,6) - Ifges(6,6);
t141 = t151 * t70;
t140 = t151 * t69;
t139 = t67 ^ 2;
t136 = m(6) * pkin(4);
t135 = m(6) * t67;
t124 = t68 * t70;
t109 = qJ(3) * t124;
t88 = t69 * t94;
t23 = -t88 + t109;
t134 = pkin(4) * t23;
t133 = pkin(4) * t70;
t128 = t67 * t69;
t127 = t67 * t70;
t126 = t68 * t67;
t125 = t68 * t69;
t102 = t126 / 0.2e1;
t106 = t128 / 0.2e1;
t113 = t139 / 0.2e1;
t110 = mrSges(6,3) * t127;
t83 = -t68 * mrSges(6,1) - t110;
t72 = t83 * t127;
t111 = mrSges(6,3) * t128;
t81 = t68 * mrSges(6,2) - t111;
t73 = t81 * t128;
t112 = mrSges(5,3) * t128;
t82 = t68 * mrSges(5,2) - t112;
t84 = t68 * mrSges(5,1) + mrSges(5,3) * t127;
t89 = mrSges(6,1) * t70 - mrSges(6,2) * t69;
t90 = mrSges(5,1) * t70 - mrSges(5,2) * t69;
t3 = t82 * t106 - t84 * t127 / 0.2e1 + t73 / 0.2e1 + t72 / 0.2e1 + (t90 + t89) * t102 + (mrSges(5,3) + mrSges(6,3)) * t150 * t113;
t120 = t3 * qJD(2);
t75 = t70 * t81;
t76 = t70 * t82;
t77 = t69 * t83;
t78 = t69 * t84;
t71 = t78 / 0.2e1 - t77 / 0.2e1 + t76 / 0.2e1 + t75 / 0.2e1;
t96 = t125 * t136;
t4 = -t71 - t96 - t147 * t125 / 0.2e1 - t122 * t124 / 0.2e1;
t119 = t4 * qJD(2);
t65 = t139 * qJ(3);
t10 = t72 + t73 - (t149 * t70 - t23 * t69) * t135;
t118 = qJD(2) * t10;
t117 = qJD(4) * t67;
t14 = -0.2e1 * (m(5) / 0.4e1 + m(6) / 0.4e1) * t148 * t126;
t116 = t14 * qJD(2);
t63 = t127 * t136;
t30 = -t67 * t89 - t63;
t115 = t30 * qJD(2);
t32 = (-t138 / 0.2e1 - t137 / 0.2e1 - 0.1e1 / 0.2e1) * t135;
t114 = t32 * qJD(2);
t108 = t68 * t121;
t91 = (pkin(4) * t69 + qJ(3)) * t139;
t38 = t53 * t70 - t108;
t80 = t69 * t53 + t109;
t74 = t70 * t80;
t79 = t139 * (t69 * mrSges(6,1) + t70 * mrSges(6,2));
t1 = -t90 * t65 - t79 * t133 - t89 * t91 + t149 * t111 - (-t87 - t108) * t81 + t67 * mrSges(5,3) * t74 - t80 * t84 - m(6) * (t133 * t91 + t134 * t68) + (t144 * t70 + t140) * t69 * t113 - (-t146 * t69 + t141) * t70 * t139 / 0.2e1 + (-t82 - t112) * t38 + (t83 + t110) * t23 + (t145 * t68 + (t146 * t70 + t140) * t67) * t106 + (t143 * t68 + (t144 * t69 - t141) * t67) * t127 / 0.2e1 + (t143 * t70 + t145 * t69) * t102;
t86 = -t3 * qJD(1) - t1 * qJD(2);
t66 = t68 ^ 2;
t6 = t139 * (t69 * mrSges(5,1) + t70 * mrSges(5,2)) + t79 + m(6) * t91 + m(5) * (-t125 * t38 + t65) + m(4) * (qJ(3) * t66 + t65) + (t66 + t139) * mrSges(4,3) + (-t77 + t75 + t76 + t78 + m(6) * (t69 * t92 + (t23 + t88) * t70) + m(5) * t74) * t68;
t85 = qJD(1) * t14 + qJD(2) * t6;
t31 = t148 * t135 / 0.2e1;
t5 = t96 / 0.2e1 + ((-mrSges(6,2) / 0.2e1 - mrSges(5,2) / 0.2e1) * t70 + (-t136 / 0.2e1 - mrSges(6,1) / 0.2e1 - mrSges(5,1) / 0.2e1) * t69) * t68 + t71;
t2 = qJD(3) * t14 - qJD(4) * t3;
t7 = [0, t2, t116, -t120 - t63 * qJD(4) + (t122 * t69 - t147 * t70) * t117, 0; t2, qJD(3) * t6 - qJD(4) * t1 - qJD(5) * t10, qJD(4) * t5 + qJD(5) * t31 + t85, t5 * qJD(3) + (-m(6) * t134 - mrSges(5,1) * t80 - t23 * mrSges(6,1) - t122 * t38) * qJD(4) + ((mrSges(6,2) * qJ(5) + t143) * t70 + (mrSges(6,3) * pkin(4) + t145) * t69) * t117 + t86, qJD(3) * t31 - t118; -t116, -qJD(4) * t4 + qJD(5) * t32 - t85, 0, -t119 + (-t122 * t70 + (-t147 - t136) * t69) * qJD(4), t114; t120, qJD(3) * t4 + qJD(5) * t30 - t86, t119, 0, t115; 0, -qJD(3) * t32 - qJD(4) * t30 + t118, -t114, -t115, 0;];
Cq = t7;
