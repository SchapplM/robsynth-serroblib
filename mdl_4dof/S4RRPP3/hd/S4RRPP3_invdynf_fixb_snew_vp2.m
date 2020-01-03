% Calculate vector of cutting forces with Newton-Euler
% S4RRPP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRPP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:23
% EndTime: 2019-12-31 16:57:24
% DurationCPUTime: 0.43s
% Computational Cost: add. (2690->120), mult. (6114->157), div. (0->0), fcn. (3505->6), ass. (0->56)
t52 = sin(qJ(2));
t54 = cos(qJ(2));
t70 = qJD(1) * qJD(2);
t42 = t52 * qJDD(1) + t54 * t70;
t43 = t54 * qJDD(1) - t52 * t70;
t72 = qJD(1) * t52;
t45 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t72;
t71 = qJD(1) * t54;
t46 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t71;
t57 = qJD(1) ^ 2;
t51 = sin(pkin(6));
t73 = cos(pkin(6));
t26 = t51 * t42 - t73 * t43;
t27 = t73 * t42 + t51 * t43;
t35 = t51 * t72 - t73 * t71;
t29 = -qJD(2) * mrSges(4,2) - t35 * mrSges(4,3);
t36 = (t51 * t54 + t73 * t52) * qJD(1);
t30 = qJD(2) * mrSges(4,1) - t36 * mrSges(4,3);
t44 = qJD(2) * pkin(2) - qJ(3) * t72;
t50 = t54 ^ 2;
t53 = sin(qJ(1));
t55 = cos(qJ(1));
t67 = t53 * g(1) - t55 * g(2);
t63 = -qJDD(1) * pkin(1) - t67;
t59 = -t43 * pkin(2) + qJDD(3) + t44 * t72 + (-qJ(3) * t50 - pkin(5)) * t57 + t63;
t31 = -qJD(2) * mrSges(5,1) + t36 * mrSges(5,2);
t32 = -t35 * mrSges(5,2) + qJD(2) * mrSges(5,3);
t61 = t27 * mrSges(5,3) + t36 * t31 - m(5) * (-0.2e1 * qJD(4) * t36 + (qJD(2) * t35 - t27) * qJ(4) + (qJD(2) * t36 + t26) * pkin(3) + t59) - t35 * t32 - t26 * mrSges(5,1);
t60 = m(4) * t59 + t26 * mrSges(4,1) + t27 * mrSges(4,2) + t35 * t29 + t36 * t30 - t61;
t81 = (t52 * t45 - t54 * t46) * qJD(1) + m(3) * (-t57 * pkin(5) + t63) - t43 * mrSges(3,1) + t42 * mrSges(3,2) + t60;
t80 = -2 * qJD(3);
t41 = (-mrSges(3,1) * t54 + mrSges(3,2) * t52) * qJD(1);
t65 = -t55 * g(1) - t53 * g(2);
t39 = -t57 * pkin(1) + qJDD(1) * pkin(5) + t65;
t76 = t52 * t39;
t77 = pkin(2) * t57;
t15 = qJDD(2) * pkin(2) - t42 * qJ(3) - t76 + (qJ(3) * t70 + t52 * t77 - g(3)) * t54;
t66 = -t52 * g(3) + t54 * t39;
t16 = t43 * qJ(3) - qJD(2) * t44 - t50 * t77 + t66;
t68 = t51 * t15 + t73 * t16 + t35 * t80;
t21 = t35 * pkin(3) - t36 * qJ(4);
t56 = qJD(2) ^ 2;
t69 = qJD(2) * t31 + qJDD(2) * mrSges(5,3) + m(5) * (-t56 * pkin(3) + qJDD(2) * qJ(4) + 0.2e1 * qJD(4) * qJD(2) - t35 * t21 + t68);
t22 = t35 * mrSges(5,1) - t36 * mrSges(5,3);
t74 = -t35 * mrSges(4,1) - t36 * mrSges(4,2) - t22;
t75 = -mrSges(4,3) - mrSges(5,2);
t7 = m(4) * t68 - qJDD(2) * mrSges(4,2) - qJD(2) * t30 + t75 * t26 + t74 * t35 + t69;
t62 = t73 * t15 - t51 * t16;
t78 = m(5) * (-qJDD(2) * pkin(3) - t56 * qJ(4) + qJDD(4) + ((2 * qJD(3)) + t21) * t36 - t62);
t8 = m(4) * t62 - t78 + (m(4) * t80 + t74) * t36 + t75 * t27 + (mrSges(4,1) + mrSges(5,1)) * qJDD(2) + (t29 + t32) * qJD(2);
t4 = m(3) * (-t54 * g(3) - t76) - t42 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t41 * t72 + qJD(2) * t46 + t51 * t7 + t73 * t8;
t5 = m(3) * t66 - qJDD(2) * mrSges(3,2) + t43 * mrSges(3,3) - qJD(2) * t45 + t41 * t71 - t51 * t8 + t73 * t7;
t79 = t54 * t4 + t52 * t5;
t6 = m(2) * t67 + qJDD(1) * mrSges(2,1) - t57 * mrSges(2,2) - t81;
t1 = m(2) * t65 - t57 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t52 * t4 + t54 * t5;
t2 = [-m(1) * g(1) + t55 * t1 - t53 * t6, t1, t5, t7, -t26 * mrSges(5,2) - t35 * t22 + t69; -m(1) * g(2) + t53 * t1 + t55 * t6, t6, t4, t8, -t61; (-m(1) - m(2)) * g(3) + t79, -m(2) * g(3) + t79, t81, t60, -qJDD(2) * mrSges(5,1) + t27 * mrSges(5,2) - qJD(2) * t32 + t36 * t22 + t78;];
f_new = t2;
