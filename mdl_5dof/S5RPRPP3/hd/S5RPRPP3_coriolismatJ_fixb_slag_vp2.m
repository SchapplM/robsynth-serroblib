% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:22
% EndTime: 2019-12-31 18:12:23
% DurationCPUTime: 0.57s
% Computational Cost: add. (1332->104), mult. (2739->124), div. (0->0), fcn. (2636->4), ass. (0->57)
t88 = m(6) / 0.2e1;
t97 = 0.2e1 * t88;
t56 = sin(pkin(7));
t57 = cos(pkin(7));
t59 = sin(qJ(3));
t84 = cos(qJ(3));
t49 = t84 * t56 + t59 * t57;
t47 = t59 * t56 - t84 * t57;
t77 = qJ(4) * t47;
t28 = pkin(3) * t49 + t77;
t96 = m(5) * t28 - t49 * mrSges(5,2);
t58 = pkin(3) + qJ(5);
t67 = -t57 * pkin(2) - pkin(1);
t62 = -t49 * qJ(4) + t67;
t12 = t58 * t47 + t62;
t66 = m(6) * t12 - t49 * mrSges(6,2) + t47 * mrSges(6,3);
t94 = m(5) + m(6);
t82 = -mrSges(5,1) - mrSges(6,1);
t92 = -mrSges(4,2) + mrSges(5,3);
t16 = t58 * t49 + t77;
t91 = Ifges(6,6) - Ifges(4,4) - Ifges(5,6);
t90 = 0.2e1 * t49;
t89 = -m(6) / 0.2e1;
t81 = pkin(6) + qJ(2);
t78 = m(6) * qJD(3);
t21 = t47 * pkin(3) + t62;
t29 = -t47 * mrSges(5,2) - t49 * mrSges(5,3);
t40 = t47 * mrSges(6,2);
t42 = t49 * mrSges(4,1);
t1 = t12 * t40 + t67 * t42 + t28 * t29 + (-t67 * mrSges(4,2) - t91 * t47) * t47 + (t12 * mrSges(6,3) + (-Ifges(4,1) + Ifges(4,2) - Ifges(5,2) + Ifges(6,2) + Ifges(5,3) - Ifges(6,3)) * t47 + t91 * t49) * t49 + (t47 * mrSges(5,3) + t96) * t21 + t66 * t16;
t76 = t1 * qJD(1);
t50 = t81 * t57;
t65 = t81 * t56;
t30 = t59 * t50 + t84 * t65;
t17 = t49 * pkin(4) + t30;
t31 = t84 * t50 - t59 * t65;
t19 = -t47 * pkin(4) + t31;
t2 = m(6) * (t17 * t49 - t19 * t47) + (m(3) * qJ(2) + mrSges(3,3)) * (t56 ^ 2 + t57 ^ 2) + (t47 ^ 2 + t49 ^ 2) * (mrSges(4,3) - t82) + (m(4) + m(5)) * (t30 * t49 - t31 * t47);
t75 = t2 * qJD(1);
t3 = t49 * mrSges(6,3) + t16 * t97 + t92 * t47 + t40 + t42 + t96;
t74 = t3 * qJD(1);
t4 = (-m(5) * t21 - t29 - t66) * t49;
t73 = t4 * qJD(1);
t7 = t66 * t47;
t72 = t7 * qJD(1);
t11 = (t89 - m(5) / 0.2e1) * t90;
t71 = t11 * qJD(1);
t23 = m(6) * t47;
t70 = t23 * qJD(1);
t52 = t94 * qJ(4) + mrSges(6,2) + mrSges(5,3);
t69 = t52 * qJD(3);
t51 = m(6) * t58 + mrSges(6,3);
t63 = t51 * qJD(3);
t10 = (-m(6) / 0.4e1 - m(5) / 0.4e1) * t90 + t94 * t49 / 0.2e1;
t8 = -t49 * mrSges(6,1) + (-t88 + t89) * t17;
t6 = m(5) * t31 + m(6) * t19 + t82 * t47;
t5 = [t2 * qJD(2) + t1 * qJD(3) + t4 * qJD(4) + t7 * qJD(5), t10 * qJD(4) + t75, t6 * qJD(4) + t8 * qJD(5) + t76 + (-t17 * mrSges(6,2) - t19 * mrSges(6,3) + (-qJ(4) * t17 - t58 * t19) * t97 + (t82 * qJ(4) + Ifges(6,4) + Ifges(5,5) - Ifges(4,6)) * t49 + (pkin(3) * mrSges(5,1) + t58 * mrSges(6,1) + Ifges(5,4) - Ifges(4,5) - Ifges(6,5)) * t47 + (-m(5) * pkin(3) - mrSges(4,1) + mrSges(5,2)) * t31 + (-m(5) * qJ(4) - t92) * t30) * qJD(3), t10 * qJD(2) + t6 * qJD(3) + t73, t8 * qJD(3) + t72; t3 * qJD(3) + t11 * qJD(4) + t23 * qJD(5) - t75, 0, t74, t71, t70; -t3 * qJD(2) - t76, -t74, t52 * qJD(4) + t51 * qJD(5), t69, t63; -t11 * qJD(2) - t73, -t71, -m(6) * qJD(5) - t69, 0, -t78; -t23 * qJD(2) - t72, -t70, m(6) * qJD(4) - t63, t78, 0;];
Cq = t5;
