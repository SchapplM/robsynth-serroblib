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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:55:46
% EndTime: 2020-01-03 11:55:48
% DurationCPUTime: 0.63s
% Computational Cost: add. (2765->102), mult. (5485->144), div. (0->0), fcn. (5087->8), ass. (0->76)
t85 = cos(qJ(2));
t116 = t85 * pkin(1);
t83 = sin(qJ(2));
t118 = pkin(1) * t83;
t79 = sin(pkin(8));
t70 = t79 * t118;
t81 = cos(pkin(8));
t58 = t81 * t116 - t70;
t78 = sin(pkin(9));
t80 = cos(pkin(9));
t82 = sin(qJ(5));
t84 = cos(qJ(5));
t66 = t84 * t78 + t82 * t80;
t41 = t66 * t58;
t65 = -t82 * t78 + t84 * t80;
t42 = t65 * t58;
t17 = -t41 * t65 + t42 * t66;
t125 = m(6) * t17;
t108 = t78 ^ 2 + t80 ^ 2;
t115 = t81 * t83;
t72 = pkin(2) + t116;
t99 = pkin(1) * t115 + t79 * t72;
t54 = qJ(4) + t99;
t49 = (-pkin(7) - t54) * t78;
t75 = t80 * pkin(7);
t50 = t80 * t54 + t75;
t23 = t84 * t49 - t82 * t50;
t24 = t82 * t49 + t84 * t50;
t112 = -t23 * t66 + t24 * t65;
t71 = t79 * pkin(2) + qJ(4);
t55 = (-pkin(7) - t71) * t78;
t56 = t80 * t71 + t75;
t43 = t84 * t55 - t82 * t56;
t44 = t82 * t55 + t84 * t56;
t110 = -t43 * t66 + t44 * t65;
t124 = t108 * mrSges(5,3);
t47 = t66 * mrSges(6,1) + mrSges(6,2) * t65;
t123 = (qJD(1) + qJD(2)) * t47;
t121 = t66 ^ 2;
t120 = -m(6) / 0.2e1;
t117 = t80 * pkin(4);
t109 = Ifges(6,5) * t65 - Ifges(6,6) * t66;
t102 = t108 * t54;
t96 = t124 + (t65 ^ 2 + t121) * mrSges(6,3);
t12 = m(5) * t102 + m(6) * t112 + t96;
t107 = qJD(1) * t12;
t45 = t47 * qJD(5);
t105 = qJD(1) * t125;
t104 = -t81 * pkin(2) - pkin(3);
t101 = t108 * t71;
t100 = t81 * t72 - t70;
t98 = -pkin(3) - t100;
t57 = (t79 * t85 + t115) * pkin(1);
t97 = 0.2e1 * (m(6) / 0.4e1 + m(5) / 0.4e1) * t57;
t15 = m(5) * t101 + m(6) * t110 + t96;
t87 = -m(5) * (t101 + t102) / 0.2e1 + (t110 + t112) * t120 - t96;
t7 = t97 + t87;
t95 = -qJD(1) * t7 + qJD(2) * t15;
t94 = -t41 * mrSges(6,1) / 0.2e1 - t42 * mrSges(6,2) / 0.2e1;
t51 = t98 - t117;
t86 = -mrSges(3,1) * t118 - mrSges(3,2) * t116 + (t41 * t66 + t42 * t65) * mrSges(6,3) + (-mrSges(4,2) + t124) * t58 + (-t80 * mrSges(5,1) - mrSges(6,1) * t65 + t78 * mrSges(5,2) + mrSges(6,2) * t66 - mrSges(4,1)) * t57;
t3 = m(6) * (-t23 * t41 + t24 * t42 + t51 * t57) + m(4) * (-t100 * t57 + t99 * t58) + m(5) * (t58 * t102 + t98 * t57) + t86;
t92 = qJD(3) * t120 * t17 - t3 * qJD(1);
t89 = (Ifges(6,4) * t65 + (Ifges(6,1) - Ifges(6,2)) * t66) * t65 - t121 * Ifges(6,4);
t4 = t51 * t47 + t89;
t91 = t4 * qJD(1);
t67 = t104 - t117;
t88 = (t51 / 0.2e1 + t67 / 0.2e1) * t47 + t89;
t1 = t88 - t94;
t6 = t67 * t47 + t89;
t90 = -t1 * qJD(1) - t6 * qJD(2);
t46 = t47 * qJD(4);
t8 = t97 - t87;
t5 = qJD(2) * t125 / 0.2e1;
t2 = t88 + t94;
t9 = [qJD(2) * t3 + qJD(4) * t12 + qJD(5) * t4, (m(6) * (-t41 * t43 + t42 * t44 + t57 * t67) + m(4) * (-t81 * t57 + t58 * t79) * pkin(2) + m(5) * (t58 * t101 + t104 * t57) + t86) * qJD(2) + t8 * qJD(4) + t2 * qJD(5) - t92, t5, qJD(2) * t8 + t107, t2 * qJD(2) + (-mrSges(6,1) * t24 - mrSges(6,2) * t23 + t109) * qJD(5) + t91; -t7 * qJD(4) + t1 * qJD(5) + t92, qJD(4) * t15 + qJD(5) * t6, -t105 / 0.2e1, t95, (-mrSges(6,1) * t44 - t43 * mrSges(6,2) + t109) * qJD(5) - t90; t5, t105 / 0.2e1, 0, 0, -t45; qJD(2) * t7 - t107 + t45, -t95 + t45, 0, 0, t123; -t1 * qJD(2) - t46 - t91, -t46 + t90, 0, -t123, 0;];
Cq = t9;
