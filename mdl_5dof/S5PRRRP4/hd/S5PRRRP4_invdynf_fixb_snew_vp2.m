% Calculate vector of cutting forces with Newton-Euler
% S5PRRRP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:39
% EndTime: 2019-12-05 16:45:41
% DurationCPUTime: 0.48s
% Computational Cost: add. (4117->103), mult. (5140->134), div. (0->0), fcn. (2719->8), ass. (0->54)
t49 = qJD(2) + qJD(3);
t47 = t49 ^ 2;
t48 = qJDD(2) + qJDD(3);
t52 = sin(pkin(8));
t68 = cos(pkin(8));
t41 = -g(1) * t68 - g(2) * t52;
t51 = -g(3) + qJDD(1);
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t63 = -t41 * t55 + t58 * t51;
t23 = qJDD(2) * pkin(2) + t63;
t60 = qJD(2) ^ 2;
t71 = t58 * t41 + t55 * t51;
t24 = -pkin(2) * t60 + t71;
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t64 = t57 * t23 - t54 * t24;
t18 = -t48 * pkin(3) - t47 * pkin(7) - t64;
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t67 = qJD(4) * t49;
t30 = t48 * t53 + t56 * t67;
t31 = t48 * t56 - t53 * t67;
t79 = t49 * t53;
t35 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t79;
t36 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t79;
t78 = t49 * t56;
t38 = mrSges(6,2) * t78 + qJD(4) * mrSges(6,3);
t69 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t78 + t38;
t74 = m(6) * (-t31 * pkin(4) - t30 * qJ(5) + (-0.2e1 * qJD(5) * t53 + (pkin(4) * t53 - qJ(5) * t56) * qJD(4)) * t49 + t18) - t31 * mrSges(6,1);
t83 = (-t69 * t56 + (t35 - t36) * t53) * t49 + m(5) * t18 - t31 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t30 + t74;
t72 = t54 * t23 + t57 * t24;
t19 = -pkin(3) * t47 + pkin(7) * t48 + t72;
t27 = (-pkin(4) * t56 - qJ(5) * t53) * t49;
t59 = qJD(4) ^ 2;
t40 = g(1) * t52 - g(2) * t68;
t77 = t56 * t40;
t80 = m(6) * (-qJDD(4) * pkin(4) - t59 * qJ(5) + t77 + qJDD(5) + (t27 * t49 + t19) * t53);
t75 = mrSges(5,3) + mrSges(6,2);
t73 = t56 * t19 - t53 * t40;
t29 = (-mrSges(5,1) * t56 + mrSges(5,2) * t53) * t49;
t28 = (-mrSges(6,1) * t56 - mrSges(6,3) * t53) * t49;
t62 = m(6) * (-pkin(4) * t59 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t27 * t78 + t73) + t28 * t78 + qJD(4) * t36 + qJDD(4) * mrSges(6,3);
t11 = m(5) * t73 - qJDD(4) * mrSges(5,2) - qJD(4) * t35 + t29 * t78 + t31 * t75 + t62;
t12 = m(5) * (-t19 * t53 - t77) - t80 + (-t28 - t29) * t79 - t75 * t30 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + t69 * qJD(4);
t6 = m(4) * t72 - mrSges(4,1) * t47 - mrSges(4,2) * t48 + t11 * t56 - t12 * t53;
t8 = m(4) * t64 + t48 * mrSges(4,1) - t47 * mrSges(4,2) - t83;
t4 = m(3) * t63 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t60 + t54 * t6 + t57 * t8;
t5 = m(3) * t71 - mrSges(3,1) * t60 - qJDD(2) * mrSges(3,2) - t54 * t8 + t57 * t6;
t66 = m(2) * t51 + t4 * t58 + t5 * t55;
t65 = m(4) * t40 - t11 * t53 - t12 * t56;
t7 = (m(2) + m(3)) * t40 + t65;
t1 = m(2) * t41 - t4 * t55 + t5 * t58;
t2 = [-m(1) * g(1) + t1 * t68 - t52 * t7, t1, t5, t6, t11, t31 * mrSges(6,2) + t62; -m(1) * g(2) + t1 * t52 + t68 * t7, t7, t4, t8, t12, -t30 * mrSges(6,3) + (-t53 * t36 - t56 * t38) * t49 + t74; -m(1) * g(3) + t66, t66, -m(3) * t40 - t65, -t65, t83, -qJDD(4) * mrSges(6,1) + t30 * mrSges(6,2) - qJD(4) * t38 + t28 * t79 + t80;];
f_new = t2;
