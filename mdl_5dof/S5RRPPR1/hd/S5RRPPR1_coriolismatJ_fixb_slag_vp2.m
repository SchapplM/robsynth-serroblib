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
% m [6x1]
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:51:22
% EndTime: 2022-01-20 09:51:24
% DurationCPUTime: 0.68s
% Computational Cost: add. (2765->103), mult. (5485->147), div. (0->0), fcn. (5087->8), ass. (0->74)
t87 = cos(qJ(2));
t120 = pkin(1) * t87;
t85 = sin(qJ(2));
t121 = pkin(1) * t85;
t81 = sin(pkin(8));
t74 = t81 * t121;
t83 = cos(pkin(8));
t61 = t83 * t120 - t74;
t80 = sin(pkin(9));
t82 = cos(pkin(9));
t84 = sin(qJ(5));
t86 = cos(qJ(5));
t68 = t80 * t86 + t82 * t84;
t41 = t68 * t61;
t67 = -t80 * t84 + t82 * t86;
t42 = t67 * t61;
t129 = m(6) * (-t41 * t67 + t42 * t68);
t108 = t80 ^ 2 + t82 ^ 2;
t128 = t108 * mrSges(5,3);
t117 = t83 * t85;
t76 = pkin(2) + t120;
t100 = pkin(1) * t117 + t81 * t76;
t57 = qJ(4) + t100;
t50 = (-pkin(7) - t57) * t80;
t77 = t82 * pkin(7);
t51 = t57 * t82 + t77;
t23 = t50 * t86 - t51 * t84;
t24 = t50 * t84 + t51 * t86;
t114 = -t23 * t68 + t24 * t67;
t75 = pkin(2) * t81 + qJ(4);
t58 = (-pkin(7) - t75) * t80;
t59 = t75 * t82 + t77;
t43 = t58 * t86 - t59 * t84;
t44 = t58 * t84 + t59 * t86;
t112 = -t43 * t68 + t44 * t67;
t47 = t68 * mrSges(6,1) + mrSges(6,2) * t67;
t127 = (qJD(1) + qJD(2)) * t47;
t125 = t68 ^ 2;
t119 = pkin(4) * t82;
t111 = t108 * t57;
t110 = Ifges(6,5) * t67 - Ifges(6,6) * t68;
t109 = t108 * t75;
t93 = t128 + (t67 ^ 2 + t125) * mrSges(6,3);
t12 = m(5) * t111 + m(6) * t114 + t93;
t107 = qJD(1) * t12;
t45 = t47 * qJD(5);
t105 = qJD(1) * t129;
t104 = -pkin(2) * t83 - pkin(3);
t102 = t108 * t61;
t101 = t76 * t83 - t74;
t99 = -pkin(3) - t101;
t15 = m(5) * t109 + m(6) * t112 + t93;
t60 = (t81 * t87 + t117) * pkin(1);
t89 = m(5) * (t109 + t111) / 0.2e1 + m(6) * (t112 + t114) / 0.2e1 + t93;
t8 = 0.2e1 * (-m(6) / 0.4e1 - m(5) / 0.4e1) * t60 + t89;
t98 = qJD(1) * t8 + qJD(2) * t15;
t97 = -t41 * mrSges(6,1) / 0.2e1 - t42 * mrSges(6,2) / 0.2e1;
t52 = t99 - t119;
t88 = -mrSges(3,1) * t121 - mrSges(3,2) * t120 + (t41 * t68 + t42 * t67) * mrSges(6,3) + (-mrSges(5,1) * t82 - mrSges(6,1) * t67 + mrSges(5,2) * t80 + mrSges(6,2) * t68 - mrSges(4,1)) * t60 + (-mrSges(4,2) + t128) * t61;
t3 = m(6) * (-t23 * t41 + t24 * t42 + t52 * t60) + m(5) * (t57 * t102 + t99 * t60) + m(4) * (t100 * t61 - t101 * t60) + t88;
t95 = -t3 * qJD(1) - qJD(3) * t129 / 0.2e1;
t91 = (Ifges(6,4) * t67 + (Ifges(6,1) - Ifges(6,2)) * t68) * t67 - t125 * Ifges(6,4);
t4 = t52 * t47 + t91;
t94 = t4 * qJD(1);
t69 = t104 - t119;
t90 = (t52 / 0.2e1 + t69 / 0.2e1) * t47 + t91;
t1 = t90 - t97;
t6 = t69 * t47 + t91;
t92 = -t1 * qJD(1) - t6 * qJD(2);
t46 = t47 * qJD(4);
t7 = t89 + (m(5) + m(6)) * t60 / 0.2e1;
t5 = qJD(2) * t129 / 0.2e1;
t2 = t90 + t97;
t9 = [qJD(2) * t3 + qJD(4) * t12 + qJD(5) * t4, (m(6) * (-t43 * t41 + t44 * t42 + t69 * t60) + m(5) * (t75 * t102 + t104 * t60) + m(4) * (-t83 * t60 + t61 * t81) * pkin(2) + t88) * qJD(2) + t7 * qJD(4) + t2 * qJD(5) - t95, t5, qJD(2) * t7 + t107, t2 * qJD(2) + (-mrSges(6,1) * t24 - mrSges(6,2) * t23 + t110) * qJD(5) + t94; t8 * qJD(4) + t1 * qJD(5) + t95, qJD(4) * t15 + qJD(5) * t6, -t105 / 0.2e1, t98, (-mrSges(6,1) * t44 - mrSges(6,2) * t43 + t110) * qJD(5) - t92; t5, t105 / 0.2e1, 0, 0, -t45; -qJD(2) * t8 - t107 + t45, -t98 + t45, 0, 0, t127; -t1 * qJD(2) - t46 - t94, -t46 + t92, 0, -t127, 0;];
Cq = t9;
