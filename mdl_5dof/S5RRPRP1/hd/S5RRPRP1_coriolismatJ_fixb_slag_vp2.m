% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:03
% EndTime: 2019-12-05 18:22:05
% DurationCPUTime: 0.84s
% Computational Cost: add. (1465->120), mult. (2959->143), div. (0->0), fcn. (2151->6), ass. (0->81)
t114 = Ifges(6,4) + Ifges(5,4);
t73 = sin(qJ(4));
t157 = t114 * t73;
t132 = pkin(4) * t73;
t75 = cos(qJ(4));
t55 = -t75 * mrSges(6,1) + t73 * mrSges(6,2);
t70 = t75 ^ 2;
t156 = -t114 * t70 - t55 * t132;
t155 = Ifges(5,1) + Ifges(6,1);
t154 = Ifges(5,2) + Ifges(6,2);
t149 = t155 * t75;
t151 = t154 * t75;
t152 = (t157 + t151 / 0.2e1 - t149 / 0.2e1 - (-t154 + t155) * t75 / 0.2e1) * t73 + t156;
t110 = t73 ^ 2 + t70;
t111 = t73 * mrSges(6,1) + t75 * mrSges(6,2);
t131 = t75 * pkin(4);
t71 = sin(pkin(8));
t74 = sin(qJ(2));
t58 = t71 * t74 * pkin(1);
t76 = cos(qJ(2));
t130 = t76 * pkin(1);
t61 = pkin(2) + t130;
t72 = cos(pkin(8));
t97 = t72 * t61 - t58;
t95 = -pkin(3) - t97;
t82 = t95 - t131;
t20 = t82 * t111;
t94 = t73 * mrSges(5,1) + t75 * mrSges(5,2);
t23 = t95 * t94;
t60 = -t72 * pkin(2) - pkin(3);
t54 = t60 - t131;
t62 = m(6) * t132;
t91 = t60 * t94;
t145 = -t20 / 0.2e1 - t23 / 0.2e1 - t54 * t111 / 0.2e1 - t91 / 0.2e1 - (t54 + t82) * t62 / 0.2e1 + t152;
t112 = t110 * mrSges(6,3);
t118 = t72 * t74;
t96 = pkin(1) * t118 + t71 * t61;
t35 = pkin(7) + t96;
t108 = qJ(5) + t35;
t25 = t108 * t75;
t21 = t25 * t75;
t24 = t108 * t73;
t11 = m(6) * (t24 * t73 + t21) + t112;
t141 = t62 + t111;
t144 = (qJD(1) + qJD(2)) * t141;
t134 = m(6) * pkin(4);
t143 = -mrSges(6,1) - t134;
t59 = t71 * pkin(2) + pkin(7);
t107 = qJ(5) + t59;
t43 = t107 * t75;
t36 = t43 * t75;
t42 = t107 * t73;
t142 = t42 * t73 + t36;
t116 = t75 * mrSges(5,1);
t140 = t73 * mrSges(5,2) - mrSges(4,1) - t116 + t55;
t138 = (-t74 * mrSges(3,1) - t76 * mrSges(3,2)) * pkin(1);
t135 = m(6) / 0.2e1;
t113 = -Ifges(5,6) - Ifges(6,6);
t40 = (t71 * t76 + t118) * pkin(1);
t41 = t72 * t130 - t58;
t98 = t110 * t41;
t3 = m(5) * t35 * t98 + t138 + (-m(4) * t97 + m(5) * t95 + m(6) * t82 + t140) * t40 + (m(4) * t96 + t110 * mrSges(5,3) - mrSges(4,2) + t11) * t41;
t109 = t3 * qJD(1);
t106 = qJD(4) * t73;
t105 = t11 * qJD(1);
t103 = t40 * t135;
t93 = -mrSges(6,3) * t131 + (Ifges(5,5) + Ifges(6,5)) * t75;
t14 = m(6) * t142 + t112;
t78 = -(t21 + t36 + (t24 + t42) * t73) * m(6) / 0.2e1 - t112;
t8 = t103 + t78;
t92 = -t8 * qJD(1) + t14 * qJD(2);
t4 = -t20 - t23 + (-t82 * t134 - t149 + t151 + t157) * t73 + t156;
t80 = t4 * qJD(1);
t1 = -((mrSges(5,2) + mrSges(6,2)) * t75 + (mrSges(5,1) - t143) * t73) * t41 / 0.2e1 + t145;
t5 = t141 * t54 - t152 + t91;
t79 = t1 * qJD(1) - t5 * qJD(2);
t39 = t141 * qJD(4);
t38 = t141 * qJD(5);
t9 = t103 - t78;
t2 = ((-mrSges(5,2) / 0.2e1 - mrSges(6,2) / 0.2e1) * t75 + (-mrSges(5,1) / 0.2e1 - t134 / 0.2e1 - mrSges(6,1) / 0.2e1) * t73) * t41 - t145;
t6 = [t3 * qJD(2) - t4 * qJD(4) + t11 * qJD(5), t2 * qJD(4) + t9 * qJD(5) + t109 + (t138 + t140 * t40 + (-mrSges(4,2) + (mrSges(5,3) + mrSges(6,3)) * t110) * t41 + 0.2e1 * (t142 * t41 + t54 * t40) * t135 + m(5) * (t60 * t40 + t59 * t98) + m(4) * (-t40 * t72 + t41 * t71) * pkin(2)) * qJD(2), 0, t2 * qJD(2) + (t24 * mrSges(6,2) - t35 * t116 + t143 * t25 + t93) * qJD(4) + (mrSges(5,2) * t35 + t113) * t106 - t80, t9 * qJD(2) + t105; -t1 * qJD(4) - t8 * qJD(5) - t109, t5 * qJD(4) + t14 * qJD(5), 0, (t42 * mrSges(6,2) - t59 * t116 + t143 * t43 + t93) * qJD(4) + (mrSges(5,2) * t59 + t113) * t106 - t79, t92; 0, 0, 0, (-t94 - t141) * qJD(4), 0; t1 * qJD(2) - t38 + t80, -t38 + t79, 0, 0, -t144; t8 * qJD(2) - t105 + t39, t39 - t92, 0, t144, 0;];
Cq = t6;
