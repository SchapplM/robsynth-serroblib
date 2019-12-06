% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRP6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:23
% EndTime: 2019-12-05 15:40:24
% DurationCPUTime: 0.35s
% Computational Cost: add. (314->81), mult. (767->110), div. (0->0), fcn. (462->4), ass. (0->50)
t40 = sin(qJ(4));
t61 = m(6) * qJD(5);
t42 = cos(qJ(4));
t20 = t40 * mrSges(6,1) - t42 * mrSges(6,3);
t47 = -t40 * mrSges(5,1) - t42 * mrSges(5,2) - t20;
t50 = t40 * pkin(4) - t42 * qJ(5);
t80 = m(6) * t50 - t47;
t82 = t80 * qJD(4) - t40 * t61;
t81 = mrSges(4,3) - t47;
t76 = -t40 ^ 2 - t42 ^ 2;
t41 = sin(qJ(2));
t79 = t41 * t42;
t78 = t41 * (m(5) / 0.4e1 + m(6) / 0.4e1);
t77 = Ifges(6,5) - Ifges(5,4);
t75 = m(6) / 0.2e1;
t74 = t41 / 0.2e1;
t19 = qJ(3) + t50;
t73 = m(6) * t19;
t72 = m(6) * t41;
t71 = pkin(4) * t42;
t70 = t40 * mrSges(5,2);
t69 = t40 * mrSges(6,3);
t67 = t42 * mrSges(5,1);
t66 = t42 * mrSges(6,1);
t44 = -pkin(2) - pkin(6);
t64 = t76 * t41 * t44;
t60 = qJ(5) * t40;
t43 = cos(qJ(2));
t4 = 0.4e1 * (0.1e1 + t76) * t43 * t78;
t59 = t4 * qJD(1);
t52 = t20 + t73;
t8 = t52 * t42;
t58 = t8 * qJD(2);
t57 = qJD(4) * t40;
t30 = m(6) * qJ(5) + mrSges(6,3);
t56 = t30 * qJD(4);
t53 = t72 / 0.2e1;
t51 = -0.2e1 * t76 * t78;
t21 = t60 + t71;
t1 = t52 * t21 + (-qJ(3) * mrSges(5,2) + t19 * mrSges(6,3) - t40 * t77) * t40 + (qJ(3) * mrSges(5,1) + t19 * mrSges(6,1) + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,3)) * t40 + t77 * t42) * t42;
t2 = (t71 / 0.2e1 + t60 / 0.2e1 - t21 / 0.2e1) * t72;
t49 = t2 * qJD(1) - t1 * qJD(2);
t6 = (-m(5) / 0.2e1 - m(6) / 0.2e1) * t41 + t51;
t7 = -t73 + (-m(4) - m(5)) * qJ(3) - t81;
t48 = t6 * qJD(1) + t7 * qJD(2);
t31 = t43 * qJ(3);
t18 = (m(6) * t44 - mrSges(6,2)) * t40;
t5 = m(4) * t41 + m(5) * t74 + t51 + t53;
t3 = t21 * t53 + (-t70 / 0.2e1 + t67 / 0.2e1 + t66 / 0.2e1 + t69 / 0.2e1 + t21 * t75) * t41 + (t67 - t70 + t66 + t69) * t74;
t9 = [t4 * qJD(2), -t61 * t79 + t5 * qJD(3) + t3 * qJD(4) + t59 + ((-mrSges(3,2) + t81) * t43 + m(4) * t31 + m(5) * (t31 - t64) + 0.2e1 * (t19 * t43 - t64) * t75 + (-m(4) * pkin(2) - mrSges(3,1) + mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t76) * t41) * qJD(2), t5 * qJD(2), t3 * qJD(2) + t82 * t43, (-qJD(2) * t79 - t43 * t57) * m(6); -t6 * qJD(3) - t2 * qJD(4) - t59, -t7 * qJD(3) + t1 * qJD(4) - t8 * qJD(5), -t48, t18 * qJD(5) + (pkin(4) * mrSges(6,2) - Ifges(6,4) - Ifges(5,5)) * t57 - t49 + ((-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t42 - t80 * t44) * qJD(4), t18 * qJD(4) - t58; t6 * qJD(2), t48, 0, -t82, m(6) * t57; t2 * qJD(2), t49, 0, t30 * qJD(5), t56; 0, t58, 0, -t56, 0;];
Cq = t9;
