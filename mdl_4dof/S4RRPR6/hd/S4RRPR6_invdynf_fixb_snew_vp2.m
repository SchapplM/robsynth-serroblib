% Calculate vector of cutting forces with Newton-Euler
% S4RRPR6
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x5]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR6_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:25
% EndTime: 2019-12-31 17:04:26
% DurationCPUTime: 0.76s
% Computational Cost: add. (6748->125), mult. (15638->172), div. (0->0), fcn. (10061->8), ass. (0->64)
t61 = sin(qJ(2));
t64 = cos(qJ(2));
t78 = qJD(1) * qJD(2);
t49 = t61 * qJDD(1) + t64 * t78;
t50 = t64 * qJDD(1) - t61 * t78;
t80 = qJD(1) * t61;
t52 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t80;
t79 = qJD(1) * t64;
t53 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t79;
t66 = qJD(1) ^ 2;
t58 = sin(pkin(7));
t59 = cos(pkin(7));
t34 = -t58 * t49 + t59 * t50;
t35 = t59 * t49 + t58 * t50;
t43 = (-t58 * t61 + t59 * t64) * qJD(1);
t36 = -qJD(2) * mrSges(4,2) + t43 * mrSges(4,3);
t44 = (t58 * t64 + t59 * t61) * qJD(1);
t37 = qJD(2) * mrSges(4,1) - t44 * mrSges(4,3);
t51 = qJD(2) * pkin(2) - qJ(3) * t80;
t57 = t64 ^ 2;
t62 = sin(qJ(1));
t65 = cos(qJ(1));
t76 = t62 * g(1) - t65 * g(2);
t71 = -qJDD(1) * pkin(1) - t76;
t68 = -t50 * pkin(2) + qJDD(3) + t51 * t80 + (-qJ(3) * t57 - pkin(5)) * t66 + t71;
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t30 = t60 * t43 + t63 * t44;
t16 = -t30 * qJD(4) + t63 * t34 - t60 * t35;
t29 = t63 * t43 - t60 * t44;
t17 = t29 * qJD(4) + t60 * t34 + t63 * t35;
t56 = qJD(2) + qJD(4);
t27 = -t56 * mrSges(5,2) + t29 * mrSges(5,3);
t28 = t56 * mrSges(5,1) - t30 * mrSges(5,3);
t38 = qJD(2) * pkin(3) - t44 * pkin(6);
t42 = t43 ^ 2;
t70 = t16 * mrSges(5,1) + t29 * t27 - m(5) * (-t34 * pkin(3) - t42 * pkin(6) + t44 * t38 + t68) - t17 * mrSges(5,2) - t30 * t28;
t69 = -m(4) * t68 + t34 * mrSges(4,1) - t35 * mrSges(4,2) + t43 * t36 - t44 * t37 + t70;
t84 = (t61 * t52 - t64 * t53) * qJD(1) + m(3) * (-t66 * pkin(5) + t71) - t50 * mrSges(3,1) + t49 * mrSges(3,2) - t69;
t48 = (-mrSges(3,1) * t64 + mrSges(3,2) * t61) * qJD(1);
t73 = -t65 * g(1) - t62 * g(2);
t46 = -t66 * pkin(1) + qJDD(1) * pkin(5) + t73;
t81 = t61 * t46;
t82 = pkin(2) * t66;
t23 = qJDD(2) * pkin(2) - t49 * qJ(3) - t81 + (qJ(3) * t78 + t61 * t82 - g(3)) * t64;
t75 = -t61 * g(3) + t64 * t46;
t24 = t50 * qJ(3) - qJD(2) * t51 - t57 * t82 + t75;
t74 = -0.2e1 * qJD(3) * t44 + t59 * t23 - t58 * t24;
t11 = (qJD(2) * t43 - t35) * pkin(6) + (t43 * t44 + qJDD(2)) * pkin(3) + t74;
t77 = 0.2e1 * qJD(3) * t43 + t58 * t23 + t59 * t24;
t12 = -t42 * pkin(3) + t34 * pkin(6) - qJD(2) * t38 + t77;
t19 = -t29 * mrSges(5,1) + t30 * mrSges(5,2);
t55 = qJDD(2) + qJDD(4);
t10 = m(5) * (t60 * t11 + t63 * t12) + t16 * mrSges(5,3) - t55 * mrSges(5,2) + t29 * t19 - t56 * t28;
t32 = -t43 * mrSges(4,1) + t44 * mrSges(4,2);
t9 = m(5) * (t63 * t11 - t60 * t12) - t17 * mrSges(5,3) + t55 * mrSges(5,1) - t30 * t19 + t56 * t27;
t6 = m(4) * t74 + qJDD(2) * mrSges(4,1) - t35 * mrSges(4,3) + qJD(2) * t36 + t60 * t10 - t44 * t32 + t63 * t9;
t7 = m(4) * t77 - qJDD(2) * mrSges(4,2) + t34 * mrSges(4,3) - qJD(2) * t37 + t63 * t10 + t43 * t32 - t60 * t9;
t4 = m(3) * (-t64 * g(3) - t81) - t49 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t48 * t80 + qJD(2) * t53 + t58 * t7 + t59 * t6;
t5 = m(3) * t75 - qJDD(2) * mrSges(3,2) + t50 * mrSges(3,3) - qJD(2) * t52 + t48 * t79 - t58 * t6 + t59 * t7;
t83 = t64 * t4 + t61 * t5;
t8 = m(2) * t76 + qJDD(1) * mrSges(2,1) - t66 * mrSges(2,2) - t84;
t1 = m(2) * t73 - t66 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t61 * t4 + t64 * t5;
t2 = [-m(1) * g(1) + t65 * t1 - t62 * t8, t1, t5, t7, t10; -m(1) * g(2) + t62 * t1 + t65 * t8, t8, t4, t6, t9; (-m(1) - m(2)) * g(3) + t83, -m(2) * g(3) + t83, t84, -t69, -t70;];
f_new = t2;
