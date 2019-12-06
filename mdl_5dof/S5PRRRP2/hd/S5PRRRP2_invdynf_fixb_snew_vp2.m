% Calculate vector of cutting forces with Newton-Euler
% S5PRRRP2
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP2_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP2_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP2_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:41
% EndTime: 2019-12-05 16:41:42
% DurationCPUTime: 0.51s
% Computational Cost: add. (4461->101), mult. (5994->134), div. (0->0), fcn. (3373->8), ass. (0->55)
t47 = qJD(2) + qJD(3);
t45 = t47 ^ 2;
t46 = qJDD(2) + qJDD(3);
t50 = sin(pkin(8));
t51 = cos(pkin(8));
t38 = t50 * g(1) - t51 * g(2);
t39 = -t51 * g(1) - t50 * g(2);
t54 = sin(qJ(2));
t57 = cos(qJ(2));
t63 = t57 * t38 - t54 * t39;
t21 = qJDD(2) * pkin(2) + t63;
t59 = qJD(2) ^ 2;
t70 = t54 * t38 + t57 * t39;
t22 = -t59 * pkin(2) + t70;
t53 = sin(qJ(3));
t56 = cos(qJ(3));
t64 = t56 * t21 - t53 * t22;
t16 = -t46 * pkin(3) - t45 * pkin(7) - t64;
t52 = sin(qJ(4));
t55 = cos(qJ(4));
t67 = qJD(4) * t47;
t28 = t52 * t46 + t55 * t67;
t29 = t55 * t46 - t52 * t67;
t78 = t47 * t52;
t34 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t78;
t35 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t78;
t77 = t47 * t55;
t37 = mrSges(6,2) * t77 + qJD(4) * mrSges(6,3);
t68 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t77 + t37;
t73 = m(6) * (-t29 * pkin(4) - t28 * qJ(5) + (-0.2e1 * qJD(5) * t52 + (pkin(4) * t52 - qJ(5) * t55) * qJD(4)) * t47 + t16) - t29 * mrSges(6,1);
t82 = (-t68 * t55 + (t34 - t35) * t52) * t47 + m(5) * t16 - t29 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t28 + t73;
t71 = t53 * t21 + t56 * t22;
t17 = -t45 * pkin(3) + t46 * pkin(7) + t71;
t25 = (-pkin(4) * t55 - qJ(5) * t52) * t47;
t58 = qJD(4) ^ 2;
t49 = -g(3) + qJDD(1);
t76 = t55 * t49;
t79 = m(6) * (-qJDD(4) * pkin(4) - t58 * qJ(5) - t76 + qJDD(5) + (t25 * t47 + t17) * t52);
t74 = mrSges(5,3) + mrSges(6,2);
t72 = t55 * t17 + t52 * t49;
t26 = (-mrSges(6,1) * t55 - mrSges(6,3) * t52) * t47;
t27 = (-mrSges(5,1) * t55 + mrSges(5,2) * t52) * t47;
t10 = m(5) * (-t52 * t17 + t76) - t79 + (-t26 - t27) * t78 - t74 * t28 + (mrSges(5,1) + mrSges(6,1)) * qJDD(4) + t68 * qJD(4);
t62 = m(6) * (-t58 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t25 * t77 + t72) + t26 * t77 + qJD(4) * t35 + qJDD(4) * mrSges(6,3);
t9 = m(5) * t72 - qJDD(4) * mrSges(5,2) - qJD(4) * t34 + t27 * t77 + t74 * t29 + t62;
t66 = m(4) * t49 + t55 * t10 + t52 * t9;
t65 = m(3) * t49 + t66;
t61 = m(2) * t49 + t65;
t6 = m(4) * t64 + t46 * mrSges(4,1) - t45 * mrSges(4,2) - t82;
t5 = m(4) * t71 - t45 * mrSges(4,1) - t46 * mrSges(4,2) - t52 * t10 + t55 * t9;
t4 = m(3) * t70 - t59 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t56 * t5 - t53 * t6;
t3 = m(3) * t63 + qJDD(2) * mrSges(3,1) - t59 * mrSges(3,2) + t53 * t5 + t56 * t6;
t2 = m(2) * t39 - t54 * t3 + t57 * t4;
t1 = m(2) * t38 + t57 * t3 + t54 * t4;
t7 = [-m(1) * g(1) - t50 * t1 + t51 * t2, t2, t4, t5, t9, t29 * mrSges(6,2) + t62; -m(1) * g(2) + t51 * t1 + t50 * t2, t1, t3, t6, t10, -t28 * mrSges(6,3) + (-t52 * t35 - t55 * t37) * t47 + t73; -m(1) * g(3) + t61, t61, t65, t66, t82, -qJDD(4) * mrSges(6,1) + t28 * mrSges(6,2) - qJD(4) * t37 + t26 * t78 + t79;];
f_new = t7;
