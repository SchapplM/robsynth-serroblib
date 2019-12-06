% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:38
% EndTime: 2019-12-05 16:15:39
% DurationCPUTime: 0.54s
% Computational Cost: add. (1872->94), mult. (3995->133), div. (0->0), fcn. (3653->6), ass. (0->68)
t77 = cos(qJ(3));
t108 = pkin(2) * t77;
t72 = sin(pkin(9));
t73 = cos(pkin(9));
t74 = sin(qJ(5));
t76 = cos(qJ(5));
t58 = t72 * t76 + t73 * t74;
t47 = t58 * t108;
t57 = -t72 * t74 + t73 * t76;
t48 = t57 * t108;
t17 = -t47 * t57 + t48 * t58;
t118 = t17 * qJD(1);
t97 = t72 ^ 2 + t73 ^ 2;
t115 = t97 * mrSges(5,3);
t88 = t97 * qJ(4);
t117 = -mrSges(4,2) + t115;
t41 = mrSges(6,1) * t58 + mrSges(6,2) * t57;
t116 = (qJD(2) + qJD(3)) * t41;
t75 = sin(qJ(3));
t106 = t75 * pkin(2);
t65 = qJ(4) + t106;
t49 = (-pkin(7) - t65) * t72;
t69 = t73 * pkin(7);
t50 = t65 * t73 + t69;
t36 = t49 * t76 - t50 * t74;
t37 = t49 * t74 + t50 * t76;
t100 = -t36 * t58 + t37 * t57;
t60 = (-pkin(7) - qJ(4)) * t72;
t62 = qJ(4) * t73 + t69;
t43 = t60 * t76 - t62 * t74;
t44 = t60 * t74 + t62 * t76;
t99 = -t43 * t58 + t44 * t57;
t114 = (t47 * t58 + t48 * t57) * mrSges(6,3);
t113 = 2 * m(6);
t112 = t58 ^ 2;
t111 = -m(6) / 0.2e1;
t107 = pkin(3) * t75;
t98 = Ifges(6,5) * t57 - Ifges(6,6) * t58;
t42 = -mrSges(6,1) * t57 + mrSges(6,2) * t58;
t66 = -pkin(4) * t73 - pkin(3);
t59 = t66 - t108;
t61 = -mrSges(5,1) * t73 + mrSges(5,2) * t72;
t89 = t97 * t65;
t6 = -mrSges(4,1) * t106 + t42 * t106 + t61 * t106 + m(6) * (t106 * t59 - t36 * t47 + t37 * t48) + m(5) * (-t107 + (t89 - t106) * t77) * pkin(2) + t114 + t117 * t108;
t95 = t6 * qJD(2);
t86 = (t57 ^ 2 + t112) * mrSges(6,3) + t115;
t12 = m(5) * t89 + m(6) * t100 + t86;
t94 = qJD(2) * t12;
t39 = t41 * qJD(5);
t91 = m(6) * t17 * qJD(2);
t87 = (m(6) / 0.2e1 + m(5) / 0.2e1) * t106;
t15 = m(5) * t88 + m(6) * t99 + t86;
t78 = -m(5) * (t89 + t88) / 0.2e1 + (t99 + t100) * t111 - t86;
t7 = t87 + t78;
t85 = -qJD(2) * t7 + qJD(3) * t15;
t84 = -t47 * mrSges(6,1) / 0.2e1 - t48 * mrSges(6,2) / 0.2e1;
t80 = (Ifges(6,4) * t57 + (Ifges(6,1) - Ifges(6,2)) * t58) * t57 - t112 * Ifges(6,4);
t3 = t59 * t41 + t80;
t82 = t3 * qJD(2);
t79 = (t59 / 0.2e1 + t66 / 0.2e1) * t41 + t80;
t1 = t79 - t84;
t4 = t66 * t41 + t80;
t81 = -qJD(2) * t1 - qJD(3) * t4;
t40 = t41 * qJD(4);
t8 = t87 - t78;
t5 = t17 * qJD(3) * t113 / 0.4e1;
t2 = t79 + t84;
t9 = [0, t5, t91 / 0.2e1, 0, -t39; t5, qJD(3) * t6 + qJD(4) * t12 + qJD(5) * t3, t95 + t8 * qJD(4) + t2 * qJD(5) + qJD(3) * t114 + (t118 / 0.4e1 + (-t43 * t47 + t44 * t48) * qJD(3) / 0.2e1) * t113 + (t117 * t77 + m(5) * (t77 * t88 - t107) + (m(6) * t66 - mrSges(4,1) + t42 + t61) * t75) * qJD(3) * pkin(2), qJD(3) * t8 + t94, t2 * qJD(3) + (-mrSges(6,1) * t37 - mrSges(6,2) * t36 + t98) * qJD(5) + t82; -t91 / 0.2e1, -qJD(4) * t7 + qJD(5) * t1 + t111 * t118 - t95, qJD(4) * t15 + qJD(5) * t4, t85, (-mrSges(6,1) * t44 - mrSges(6,2) * t43 + t98) * qJD(5) - t81; 0, qJD(3) * t7 + t39 - t94, -t85 + t39, 0, t116; 0, -qJD(3) * t1 - t40 - t82, -t40 + t81, -t116, 0;];
Cq = t9;
