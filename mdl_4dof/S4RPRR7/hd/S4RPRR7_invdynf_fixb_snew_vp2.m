% Calculate vector of cutting forces with Newton-Euler
% S4RPRR7
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR7_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR7_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:43
% EndTime: 2019-12-31 16:53:44
% DurationCPUTime: 0.56s
% Computational Cost: add. (4538->112), mult. (10701->148), div. (0->0), fcn. (7179->8), ass. (0->65)
t59 = qJD(1) ^ 2;
t51 = cos(pkin(7));
t49 = t51 ^ 2;
t50 = sin(pkin(7));
t78 = t50 ^ 2 + t49;
t82 = t78 * mrSges(3,3);
t53 = sin(qJ(3));
t56 = cos(qJ(3));
t65 = t50 * t53 - t51 * t56;
t41 = t65 * qJD(1);
t66 = t50 * t56 + t51 * t53;
t42 = t66 * qJD(1);
t75 = t42 * qJD(3);
t32 = -t65 * qJDD(1) - t75;
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t70 = -t57 * g(1) - t54 * g(2);
t43 = -t59 * pkin(1) + qJDD(1) * qJ(2) + t70;
t68 = -t51 * mrSges(3,1) + t50 * mrSges(3,2);
t64 = qJDD(1) * mrSges(3,3) + t59 * t68;
t30 = t41 * pkin(3) - t42 * pkin(6);
t58 = qJD(3) ^ 2;
t74 = qJD(1) * qJD(2);
t72 = -t51 * g(3) - 0.2e1 * t50 * t74;
t77 = pkin(5) * qJDD(1);
t80 = pkin(2) * t59;
t23 = (t51 * t80 - t43 - t77) * t50 + t72;
t71 = -t50 * g(3) + (t43 + 0.2e1 * t74) * t51;
t24 = -t49 * t80 + t51 * t77 + t71;
t79 = t53 * t23 + t56 * t24;
t14 = -t58 * pkin(3) + qJDD(3) * pkin(6) - t41 * t30 + t79;
t76 = t41 * qJD(3);
t33 = t66 * qJDD(1) - t76;
t73 = t54 * g(1) - t57 * g(2);
t69 = qJDD(2) - t73;
t60 = (-pkin(2) * t51 - pkin(1)) * qJDD(1) + (-t78 * pkin(5) - qJ(2)) * t59 + t69;
t15 = (-t33 + t76) * pkin(6) + (-t32 + t75) * pkin(3) + t60;
t52 = sin(qJ(4));
t55 = cos(qJ(4));
t34 = t55 * qJD(3) - t52 * t42;
t17 = t34 * qJD(4) + t52 * qJDD(3) + t55 * t33;
t35 = t52 * qJD(3) + t55 * t42;
t18 = -t34 * mrSges(5,1) + t35 * mrSges(5,2);
t39 = qJD(4) + t41;
t20 = -t39 * mrSges(5,2) + t34 * mrSges(5,3);
t29 = qJDD(4) - t32;
t11 = m(5) * (-t52 * t14 + t55 * t15) - t17 * mrSges(5,3) + t29 * mrSges(5,1) - t35 * t18 + t39 * t20;
t16 = -t35 * qJD(4) + t55 * qJDD(3) - t52 * t33;
t21 = t39 * mrSges(5,1) - t35 * mrSges(5,3);
t12 = m(5) * (t55 * t14 + t52 * t15) + t16 * mrSges(5,3) - t29 * mrSges(5,2) + t34 * t18 - t39 * t21;
t27 = t41 * mrSges(4,1) + t42 * mrSges(4,2);
t37 = qJD(3) * mrSges(4,1) - t42 * mrSges(4,3);
t7 = m(4) * t79 - qJDD(3) * mrSges(4,2) + t32 * mrSges(4,3) - qJD(3) * t37 - t52 * t11 + t55 * t12 - t41 * t27;
t36 = -qJD(3) * mrSges(4,2) - t41 * mrSges(4,3);
t67 = t56 * t23 - t53 * t24;
t61 = m(5) * (-qJDD(3) * pkin(3) - t58 * pkin(6) + t42 * t30 - t67) - t16 * mrSges(5,1) + t17 * mrSges(5,2) - t34 * t20 + t35 * t21;
t8 = m(4) * t67 + qJDD(3) * mrSges(4,1) - t33 * mrSges(4,3) + qJD(3) * t36 - t42 * t27 - t61;
t4 = m(3) * t72 + t53 * t7 + t56 * t8 + (-m(3) * t43 - t64) * t50;
t5 = m(3) * t71 + t64 * t51 - t53 * t8 + t56 * t7;
t81 = t51 * t4 + t50 * t5;
t63 = -m(4) * t60 + t32 * mrSges(4,1) - t33 * mrSges(4,2) - t55 * t11 - t52 * t12 - t41 * t36 - t42 * t37;
t62 = m(3) * (-qJDD(1) * pkin(1) - t59 * qJ(2) + t69) - t63;
t6 = m(2) * t73 + (-mrSges(2,2) + t82) * t59 + (mrSges(2,1) - t68) * qJDD(1) - t62;
t1 = m(2) * t70 - t59 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t50 * t4 + t51 * t5;
t2 = [-m(1) * g(1) + t57 * t1 - t54 * t6, t1, t5, t7, t12; -m(1) * g(2) + t54 * t1 + t57 * t6, t6, t4, t8, t11; (-m(1) - m(2)) * g(3) + t81, -m(2) * g(3) + t81, t68 * qJDD(1) - t59 * t82 + t62, -t63, t61;];
f_new = t2;
