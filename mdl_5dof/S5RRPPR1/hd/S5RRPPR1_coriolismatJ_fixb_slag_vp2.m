% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:13
% EndTime: 2019-12-05 18:18:15
% DurationCPUTime: 0.75s
% Computational Cost: add. (2765->112), mult. (5485->155), div. (0->0), fcn. (5087->8), ass. (0->80)
t84 = sin(pkin(9));
t86 = cos(pkin(9));
t88 = sin(qJ(5));
t90 = cos(qJ(5));
t71 = -t88 * t84 + t90 * t86;
t139 = t71 ^ 2;
t91 = cos(qJ(2));
t126 = t91 * pkin(1);
t89 = sin(qJ(2));
t128 = pkin(1) * t89;
t85 = sin(pkin(8));
t76 = t85 * t128;
t87 = cos(pkin(8));
t64 = t87 * t126 - t76;
t72 = t90 * t84 + t88 * t86;
t46 = t72 * t64;
t47 = t71 * t64;
t136 = m(6) * (-t46 * t71 + t47 * t72);
t124 = t139 * Ifges(6,4);
t138 = -Ifges(6,4) * t72 + (Ifges(6,1) - Ifges(6,2)) * t71;
t115 = t84 ^ 2 + t86 ^ 2;
t122 = t87 * t89;
t78 = pkin(2) + t126;
t105 = pkin(1) * t122 + t85 * t78;
t60 = qJ(4) + t105;
t55 = (-pkin(7) - t60) * t84;
t81 = t86 * pkin(7);
t56 = t86 * t60 + t81;
t26 = t90 * t55 - t88 * t56;
t27 = t88 * t55 + t90 * t56;
t119 = -t26 * t72 + t27 * t71;
t77 = t85 * pkin(2) + qJ(4);
t61 = (-pkin(7) - t77) * t84;
t62 = t86 * t77 + t81;
t49 = t90 * t61 - t88 * t62;
t50 = t88 * t61 + t90 * t62;
t117 = -t49 * t72 + t50 * t71;
t133 = t115 * mrSges(5,3);
t53 = t72 * mrSges(6,1) + mrSges(6,2) * t71;
t132 = (qJD(1) + qJD(2)) * t53;
t127 = t86 * pkin(4);
t123 = t72 * mrSges(6,3);
t116 = Ifges(6,5) * t71 - Ifges(6,6) * t72;
t102 = t133 + (t72 ^ 2 + t139) * mrSges(6,3);
t108 = t115 * t60;
t12 = m(5) * t108 + m(6) * t119 + t102;
t113 = qJD(1) * t12;
t51 = t53 * qJD(5);
t111 = qJD(1) * t136;
t110 = -t87 * pkin(2) - pkin(3);
t107 = t115 * t77;
t106 = t87 * t78 - t76;
t104 = -pkin(3) - t106;
t63 = (t85 * t91 + t122) * pkin(1);
t103 = 0.2e1 * (m(6) / 0.4e1 + m(5) / 0.4e1) * t63;
t16 = m(5) * t107 + m(6) * t117 + t102;
t93 = -m(5) * (t107 + t108) / 0.2e1 - m(6) * (t117 + t119) / 0.2e1 - t102;
t7 = t103 + t93;
t101 = -qJD(1) * t7 + qJD(2) * t16;
t100 = -t46 * mrSges(6,1) / 0.2e1 - t47 * mrSges(6,2) / 0.2e1;
t57 = t104 - t127;
t92 = t47 * t71 * mrSges(6,3) - mrSges(3,1) * t128 - mrSges(3,2) * t126 + t46 * t123 + (-mrSges(4,2) + t133) * t64 + (-t86 * mrSges(5,1) - t71 * mrSges(6,1) + t84 * mrSges(5,2) + t72 * mrSges(6,2) - mrSges(4,1)) * t63;
t3 = m(6) * (-t26 * t46 + t27 * t47 + t57 * t63) + m(4) * (t105 * t64 - t106 * t63) + m(5) * (t104 * t63 + t64 * t108) + t92;
t97 = -t3 * qJD(1) - qJD(3) * t136 / 0.2e1;
t22 = t27 * t123;
t35 = t57 * t53;
t4 = -t124 + t22 - t35 + (-t27 * mrSges(6,3) - t138) * t72;
t96 = t4 * qJD(1);
t36 = t50 * t123;
t73 = t110 - t127;
t48 = t73 * t53;
t94 = t124 + (t27 / 0.2e1 + t50 / 0.2e1) * t123 - t22 / 0.2e1 + t35 / 0.2e1 - t36 / 0.2e1 + t48 / 0.2e1 + t138 * t72;
t1 = t94 - t100;
t6 = -t124 + t36 - t48 + (-t50 * mrSges(6,3) - t138) * t72;
t95 = -t1 * qJD(1) + t6 * qJD(2);
t52 = t53 * qJD(4);
t8 = t103 - t93;
t5 = qJD(2) * t136 / 0.2e1;
t2 = t94 + t100;
t9 = [qJD(2) * t3 + qJD(4) * t12 - qJD(5) * t4, (m(6) * (-t49 * t46 + t50 * t47 + t73 * t63) + m(4) * (-t87 * t63 + t64 * t85) * pkin(2) + m(5) * (t64 * t107 + t110 * t63) + t92) * qJD(2) + t8 * qJD(4) + t2 * qJD(5) - t97, t5, qJD(2) * t8 + t113, t2 * qJD(2) + (-mrSges(6,1) * t27 - mrSges(6,2) * t26 + t116) * qJD(5) - t96; -t7 * qJD(4) + t1 * qJD(5) + t97, qJD(4) * t16 - qJD(5) * t6, -t111 / 0.2e1, t101, (-mrSges(6,1) * t50 - t49 * mrSges(6,2) + t116) * qJD(5) - t95; t5, t111 / 0.2e1, 0, 0, -t51; qJD(2) * t7 - t113 + t51, -t101 + t51, 0, 0, t132; -t1 * qJD(2) - t52 + t96, -t52 + t95, 0, -t132, 0;];
Cq = t9;
