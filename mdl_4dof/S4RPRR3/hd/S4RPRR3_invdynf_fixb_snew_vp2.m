% Calculate vector of cutting forces with Newton-Euler
% S4RPRR3
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR3_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR3_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:06
% EndTime: 2019-12-31 16:49:07
% DurationCPUTime: 0.39s
% Computational Cost: add. (3157->97), mult. (6130->134), div. (0->0), fcn. (3457->8), ass. (0->54)
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t65 = qJD(1) * qJD(3);
t62 = t52 * t65;
t33 = t49 * qJDD(1) + t62;
t34 = t52 * qJDD(1) - t49 * t65;
t67 = qJD(1) * t49;
t35 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t67;
t66 = qJD(1) * t52;
t36 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t66;
t54 = qJD(1) ^ 2;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t29 = (t48 * t52 + t49 * t51) * qJD(1);
t16 = -t29 * qJD(4) - t48 * t33 + t51 * t34;
t28 = (-t48 * t49 + t51 * t52) * qJD(1);
t17 = t28 * qJD(4) + t51 * t33 + t48 * t34;
t43 = qJD(3) + qJD(4);
t23 = -t43 * mrSges(5,2) + t28 * mrSges(5,3);
t24 = t43 * mrSges(5,1) - t29 * mrSges(5,3);
t37 = qJD(3) * pkin(3) - pkin(6) * t67;
t44 = t52 ^ 2;
t50 = sin(qJ(1));
t53 = cos(qJ(1));
t63 = t50 * g(1) - t53 * g(2);
t30 = qJDD(1) * pkin(1) + t63;
t59 = -t53 * g(1) - t50 * g(2);
t32 = -t54 * pkin(1) + t59;
t46 = sin(pkin(7));
t47 = cos(pkin(7));
t60 = t47 * t30 - t46 * t32;
t57 = -qJDD(1) * pkin(2) - t60;
t56 = t16 * mrSges(5,1) + t28 * t23 - m(5) * (t37 * t67 - t34 * pkin(3) + (-pkin(6) * t44 - pkin(5)) * t54 + t57) - t17 * mrSges(5,2) - t29 * t24;
t70 = (t49 * t35 - t52 * t36) * qJD(1) + m(4) * (-t54 * pkin(5) + t57) - t34 * mrSges(4,1) + t33 * mrSges(4,2) - t56;
t68 = t46 * t30 + t47 * t32;
t21 = -t54 * pkin(2) + qJDD(1) * pkin(5) + t68;
t45 = -g(3) + qJDD(2);
t69 = t52 * t21 + t49 * t45;
t61 = -t49 * t21 + t52 * t45;
t11 = (-t33 + t62) * pkin(6) + (t49 * t52 * t54 + qJDD(3)) * pkin(3) + t61;
t12 = -t44 * t54 * pkin(3) + t34 * pkin(6) - qJD(3) * t37 + t69;
t22 = -t28 * mrSges(5,1) + t29 * mrSges(5,2);
t42 = qJDD(3) + qJDD(4);
t10 = m(5) * (t48 * t11 + t51 * t12) + t16 * mrSges(5,3) - t42 * mrSges(5,2) + t28 * t22 - t43 * t24;
t31 = (-mrSges(4,1) * t52 + mrSges(4,2) * t49) * qJD(1);
t9 = m(5) * (t51 * t11 - t48 * t12) - t17 * mrSges(5,3) + t42 * mrSges(5,1) - t29 * t22 + t43 * t23;
t6 = m(4) * t61 + qJDD(3) * mrSges(4,1) - t33 * mrSges(4,3) + qJD(3) * t36 + t48 * t10 - t31 * t67 + t51 * t9;
t7 = m(4) * t69 - qJDD(3) * mrSges(4,2) + t34 * mrSges(4,3) - qJD(3) * t35 + t51 * t10 + t31 * t66 - t48 * t9;
t64 = m(3) * t45 + t49 * t7 + t52 * t6;
t8 = m(3) * t60 + qJDD(1) * mrSges(3,1) - t54 * mrSges(3,2) - t70;
t3 = m(3) * t68 - t54 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t49 * t6 + t52 * t7;
t2 = m(2) * t59 - t54 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t47 * t3 - t46 * t8;
t1 = m(2) * t63 + qJDD(1) * mrSges(2,1) - t54 * mrSges(2,2) + t46 * t3 + t47 * t8;
t4 = [-m(1) * g(1) - t50 * t1 + t53 * t2, t2, t3, t7, t10; -m(1) * g(2) + t53 * t1 + t50 * t2, t1, t8, t6, t9; (-m(1) - m(2)) * g(3) + t64, -m(2) * g(3) + t64, t64, t70, -t56;];
f_new = t4;
