% Calculate vector of cutting forces with Newton-Euler
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:22
% EndTime: 2019-12-05 15:40:23
% DurationCPUTime: 0.31s
% Computational Cost: add. (1610->100), mult. (2807->126), div. (0->0), fcn. (1279->6), ass. (0->50)
t48 = sin(pkin(7));
t71 = cos(pkin(7));
t36 = -t71 * g(1) - t48 * g(2);
t47 = -g(3) + qJDD(1);
t50 = sin(qJ(2));
t52 = cos(qJ(2));
t72 = t52 * t36 + t50 * t47;
t80 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t72;
t79 = -pkin(2) - pkin(6);
t54 = qJD(2) ^ 2;
t63 = -t50 * t36 + t52 * t47;
t58 = -t54 * qJ(3) + qJDD(3) - t63;
t17 = t79 * qJDD(2) + t58;
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t28 = (pkin(4) * t49 - qJ(5) * t51) * qJD(2);
t53 = qJD(4) ^ 2;
t35 = t48 * g(1) - t71 * g(2);
t77 = t49 * t35;
t78 = m(6) * (-qJDD(4) * pkin(4) - t53 * qJ(5) - t77 + qJDD(5) + (qJD(2) * t28 - t17) * t51);
t76 = mrSges(3,1) - mrSges(4,2);
t75 = (-mrSges(3,2) + mrSges(4,3));
t74 = -mrSges(5,3) - mrSges(6,2);
t73 = t49 * t17 - t51 * t35;
t70 = qJD(2) * t49;
t69 = qJD(2) * t51;
t68 = qJD(2) * qJD(4);
t31 = t49 * qJDD(2) + t51 * t68;
t38 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t69;
t29 = (mrSges(6,1) * t49 - mrSges(6,3) * t51) * qJD(2);
t62 = qJD(2) * (-t29 - (mrSges(5,1) * t49 + mrSges(5,2) * t51) * qJD(2));
t39 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t69;
t65 = m(6) * (-t53 * pkin(4) + qJDD(4) * qJ(5) + (2 * qJD(5) * qJD(4)) - t28 * t70 + t73) + qJD(4) * t39 + qJDD(4) * mrSges(6,3);
t8 = m(5) * t73 - qJDD(4) * mrSges(5,2) - qJD(4) * t38 + t74 * t31 + t49 * t62 + t65;
t32 = t51 * qJDD(2) - t49 * t68;
t37 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t70;
t40 = -mrSges(6,2) * t70 + qJD(4) * mrSges(6,3);
t9 = m(5) * (t51 * t17 + t77) - t78 + t74 * t32 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + (t37 + t40) * qJD(4) + t51 * t62;
t60 = -m(4) * (-qJDD(2) * pkin(2) + t58) - t49 * t8 - t51 * t9;
t3 = m(3) * t63 + t76 * qJDD(2) + (t75 * t54) + t60;
t59 = t79 * t54 - t80;
t57 = -t32 * mrSges(6,3) - t39 * t69 + m(6) * (t31 * pkin(4) - t32 * qJ(5) + (-0.2e1 * qJD(5) * t51 + (pkin(4) * t51 + qJ(5) * t49) * qJD(4)) * qJD(2) + t59) + t40 * t70 + t31 * mrSges(6,1);
t56 = m(5) * t59 + t31 * mrSges(5,1) + t32 * mrSges(5,2) + t37 * t70 + t38 * t69 + t57;
t55 = -m(4) * (t54 * pkin(2) + t80) + t56;
t6 = m(3) * t72 + t75 * qJDD(2) - t76 * t54 + t55;
t66 = m(2) * t47 + t52 * t3 + t50 * t6;
t61 = m(4) * t35 + t49 * t9 - t51 * t8;
t4 = (m(2) + m(3)) * t35 + t61;
t1 = m(2) * t36 - t50 * t3 + t52 * t6;
t2 = [-m(1) * g(1) + t71 * t1 - t48 * t4, t1, t6, -t61, t8, -t31 * mrSges(6,2) - t29 * t70 + t65; -m(1) * g(2) + t48 * t1 + t71 * t4, t4, t3, -(t54 * mrSges(4,2)) - qJDD(2) * mrSges(4,3) - t55, t9, t57; -m(1) * g(3) + t66, t66, -m(3) * t35 - t61, qJDD(2) * mrSges(4,2) - t54 * mrSges(4,3) - t60, t56, -qJDD(4) * mrSges(6,1) + t32 * mrSges(6,2) - qJD(4) * t40 + t29 * t69 + t78;];
f_new = t2;
