% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:41
% EndTime: 2019-12-31 18:08:42
% DurationCPUTime: 0.66s
% Computational Cost: add. (1509->89), mult. (2935->109), div. (0->0), fcn. (2794->6), ass. (0->52)
t95 = sin(qJ(3));
t83 = t95 * pkin(3);
t107 = m(5) * t83;
t87 = sin(pkin(8));
t88 = cos(pkin(8));
t96 = cos(qJ(3));
t43 = t87 * t95 - t88 * t96;
t45 = -t87 * t96 - t88 * t95;
t82 = -cos(pkin(7)) * pkin(1) - pkin(2);
t48 = -t96 * pkin(3) + t82;
t22 = t43 * pkin(4) + t45 * qJ(5) + t48;
t106 = m(6) * t22 + mrSges(6,1) * t43 + mrSges(6,3) * t45;
t105 = -t45 * mrSges(5,1) - t43 * mrSges(5,2);
t79 = t87 * pkin(3);
t50 = t79 + qJ(5);
t47 = m(6) * t50 + mrSges(6,3);
t102 = -Ifges(6,5) + Ifges(5,4);
t29 = m(6) * t45;
t101 = qJD(5) * t29;
t55 = sin(pkin(7)) * pkin(1) + pkin(6);
t81 = t95 * t55;
t64 = -t95 * qJ(4) - t81;
t49 = t96 * t55;
t90 = t96 * qJ(4) + t49;
t25 = -t88 * t64 + t87 * t90;
t100 = t87 * t64 + t88 * t90;
t98 = m(6) / 0.2e1;
t97 = m(5) * pkin(3);
t94 = t43 * mrSges(6,2);
t62 = (-t87 * t43 + t88 * t45) * t97;
t80 = t88 * pkin(3);
t56 = -t80 - pkin(4);
t68 = m(6) * (-t50 * t43 - t56 * t45);
t60 = t68 / 0.2e1 + t62 / 0.2e1;
t30 = -t45 * pkin(4) + t43 * qJ(5) + t83;
t63 = t30 * t98 + t107 / 0.2e1;
t39 = t43 * mrSges(6,3);
t67 = -t45 * mrSges(6,1) + t105 + t39;
t7 = -t60 + t63 + t67;
t89 = t7 * qJD(1);
t13 = t106 * t45;
t86 = qJD(1) * t13;
t84 = t29 * qJD(1);
t65 = t95 * mrSges(4,1) + t96 * mrSges(4,2);
t1 = t82 * t65 + t22 * t39 + (-t22 * mrSges(6,1) - mrSges(5,2) * t83 - t102 * t45) * t45 + (mrSges(5,1) * t83 + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t45 + t102 * t43) * t43 + (-Ifges(4,2) + Ifges(4,1)) * t96 * t95 + (-t95 ^ 2 + t96 ^ 2) * Ifges(4,4) + (t105 + t107) * t48 + t106 * t30;
t72 = t1 * qJD(1);
t4 = (t43 ^ 2 + t45 ^ 2) * (mrSges(5,3) + mrSges(6,2)) + (m(5) + m(6)) * (-t100 * t43 - t25 * t45);
t70 = qJD(1) * t4;
t69 = qJD(3) * t47;
t11 = 0.2e1 * t100 * t98 - t94;
t8 = t60 + t63;
t2 = [qJD(3) * t1 + qJD(4) * t4 + qJD(5) * t13, 0, t8 * qJD(4) + t11 * qJD(5) + t72 + (-mrSges(4,1) * t49 + mrSges(4,2) * t81 + Ifges(4,5) * t96 - Ifges(4,6) * t95 - t56 * t94 - (t87 * t97 - mrSges(5,2) + t47) * t25 + (m(6) * t56 - t88 * t97 - mrSges(5,1) - mrSges(6,1)) * t100 + (mrSges(5,3) * t80 - Ifges(6,4) - Ifges(5,5)) * t43 + (t50 * mrSges(6,2) + mrSges(5,3) * t79 + Ifges(5,6) - Ifges(6,6)) * t45) * qJD(3), qJD(3) * t8 + t70, qJD(3) * t11 + t86; 0, 0, (t62 + t68 - t65 - t67) * qJD(3) - t101, 0, -qJD(3) * t29; -qJD(4) * t7 - t72, 0, t47 * qJD(5), -t89, t69; qJD(3) * t7 + t101 - t70, 0, t89, 0, t84; -qJD(4) * t29 - t86, 0, -t69, -t84, 0;];
Cq = t2;
