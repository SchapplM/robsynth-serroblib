% Calculate vector of cutting forces with Newton-Euler
% S5PRPRP4
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
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRPRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:06
% EndTime: 2019-12-05 15:35:08
% DurationCPUTime: 0.45s
% Computational Cost: add. (3009->101), mult. (5140->133), div. (0->0), fcn. (2719->8), ass. (0->52)
t59 = qJD(2) ^ 2;
t51 = sin(pkin(7));
t53 = cos(pkin(7));
t38 = -t53 * g(1) - t51 * g(2);
t49 = -g(3) + qJDD(1);
t55 = sin(qJ(2));
t57 = cos(qJ(2));
t62 = -t55 * t38 + t57 * t49;
t23 = qJDD(2) * pkin(2) + t62;
t71 = t57 * t38 + t55 * t49;
t24 = -t59 * pkin(2) + t71;
t50 = sin(pkin(8));
t52 = cos(pkin(8));
t63 = t52 * t23 - t50 * t24;
t18 = -qJDD(2) * pkin(3) - t59 * pkin(6) - t63;
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t66 = qJD(2) * qJD(4);
t33 = t54 * qJDD(2) + t56 * t66;
t34 = t56 * qJDD(2) - t54 * t66;
t68 = qJD(2) * t54;
t39 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t68;
t40 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t68;
t67 = qJD(2) * t56;
t42 = mrSges(6,2) * t67 + qJD(4) * mrSges(6,3);
t69 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t67 + t42;
t74 = m(6) * (-t34 * pkin(4) - t33 * qJ(5) + (-0.2e1 * qJD(5) * t54 + (pkin(4) * t54 - qJ(5) * t56) * qJD(4)) * qJD(2) + t18) - t34 * mrSges(6,1);
t81 = (-t69 * t56 + (t39 - t40) * t54) * qJD(2) + m(5) * t18 - t34 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t33 + t74;
t72 = t50 * t23 + t52 * t24;
t19 = -t59 * pkin(3) + qJDD(2) * pkin(6) + t72;
t30 = (-pkin(4) * t56 - qJ(5) * t54) * qJD(2);
t58 = qJD(4) ^ 2;
t37 = t51 * g(1) - t53 * g(2);
t36 = qJDD(3) - t37;
t77 = t56 * t36;
t78 = m(6) * (-qJDD(4) * pkin(4) - t58 * qJ(5) - t77 + qJDD(5) + (qJD(2) * t30 + t19) * t54);
t75 = mrSges(5,3) + mrSges(6,2);
t73 = t56 * t19 + t54 * t36;
t32 = (-mrSges(5,1) * t56 + mrSges(5,2) * t54) * qJD(2);
t31 = (-mrSges(6,1) * t56 - mrSges(6,3) * t54) * qJD(2);
t61 = m(6) * (-t58 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t30 * t67 + t73) + t31 * t67 + qJD(4) * t40 + qJDD(4) * mrSges(6,3);
t11 = m(5) * t73 - qJDD(4) * mrSges(5,2) - qJD(4) * t39 + t32 * t67 + t75 * t34 + t61;
t12 = m(5) * (-t54 * t19 + t77) - t78 - t75 * t33 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + t69 * qJD(4) + (-t31 - t32) * t68;
t6 = m(4) * t72 - t59 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t56 * t11 - t54 * t12;
t8 = m(4) * t63 + qJDD(2) * mrSges(4,1) - t59 * mrSges(4,2) - t81;
t4 = m(3) * t62 + qJDD(2) * mrSges(3,1) - t59 * mrSges(3,2) + t50 * t6 + t52 * t8;
t5 = m(3) * t71 - t59 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t50 * t8 + t52 * t6;
t65 = m(2) * t49 + t57 * t4 + t55 * t5;
t64 = m(4) * t36 + t54 * t11 + t56 * t12;
t7 = (m(2) + m(3)) * t37 - t64;
t1 = m(2) * t38 - t55 * t4 + t57 * t5;
t2 = [-m(1) * g(1) + t53 * t1 - t51 * t7, t1, t5, t6, t11, t34 * mrSges(6,2) + t61; -m(1) * g(2) + t51 * t1 + t53 * t7, t7, t4, t8, t12, -t33 * mrSges(6,3) + (-t54 * t40 - t56 * t42) * qJD(2) + t74; -m(1) * g(3) + t65, t65, -m(3) * t37 + t64, t64, t81, -qJDD(4) * mrSges(6,1) + t33 * mrSges(6,2) - qJD(4) * t42 + t31 * t68 + t78;];
f_new = t2;
