% Calculate vector of cutting forces with Newton-Euler
% S4PRRR7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR7_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR7_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:05
% EndTime: 2019-12-31 16:36:07
% DurationCPUTime: 0.48s
% Computational Cost: add. (4144->95), mult. (7818->135), div. (0->0), fcn. (5071->10), ass. (0->60)
t44 = sin(pkin(8));
t46 = cos(pkin(8));
t37 = -g(1) * t46 - g(2) * t44;
t43 = -g(3) + qJDD(1);
t45 = sin(pkin(4));
t50 = sin(qJ(2));
t53 = cos(qJ(2));
t36 = g(1) * t44 - g(2) * t46;
t47 = cos(pkin(4));
t70 = t36 * t47;
t73 = -t50 * t37 + (t43 * t45 + t70) * t53;
t49 = sin(qJ(3));
t52 = cos(qJ(3));
t33 = (-pkin(3) * t52 - pkin(7) * t49) * qJD(2);
t54 = qJD(3) ^ 2;
t65 = t52 * qJD(2);
t55 = qJD(2) ^ 2;
t69 = t45 * t50;
t63 = t53 * t37 + t43 * t69 + t50 * t70;
t18 = -pkin(2) * t55 + qJDD(2) * pkin(6) + t63;
t25 = -t36 * t45 + t43 * t47;
t67 = t52 * t18 + t49 * t25;
t14 = -pkin(3) * t54 + qJDD(3) * pkin(7) + t33 * t65 + t67;
t17 = -qJDD(2) * pkin(2) - t55 * pkin(6) - t73;
t64 = qJD(2) * qJD(3);
t60 = t52 * t64;
t34 = qJDD(2) * t49 + t60;
t61 = t49 * t64;
t35 = qJDD(2) * t52 - t61;
t15 = (-t34 - t60) * pkin(7) + (-t35 + t61) * pkin(3) + t17;
t48 = sin(qJ(4));
t51 = cos(qJ(4));
t66 = qJD(2) * t49;
t30 = qJD(3) * t51 - t48 * t66;
t20 = qJD(4) * t30 + qJDD(3) * t48 + t34 * t51;
t31 = qJD(3) * t48 + t51 * t66;
t21 = -mrSges(5,1) * t30 + mrSges(5,2) * t31;
t41 = qJD(4) - t65;
t23 = -mrSges(5,2) * t41 + mrSges(5,3) * t30;
t27 = qJDD(4) - t35;
t11 = m(5) * (-t14 * t48 + t15 * t51) - t20 * mrSges(5,3) + t27 * mrSges(5,1) - t31 * t21 + t41 * t23;
t19 = -qJD(4) * t31 + qJDD(3) * t51 - t34 * t48;
t24 = mrSges(5,1) * t41 - mrSges(5,3) * t31;
t12 = m(5) * (t14 * t51 + t15 * t48) + t19 * mrSges(5,3) - t27 * mrSges(5,2) + t30 * t21 - t41 * t24;
t38 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t66;
t39 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t65;
t72 = m(4) * t17 - t35 * mrSges(4,1) + t34 * mrSges(4,2) + t51 * t11 + t48 * t12 + (t38 * t49 - t39 * t52) * qJD(2);
t8 = m(3) * t73 + qJDD(2) * mrSges(3,1) - t55 * mrSges(3,2) - t72;
t71 = t53 * t8;
t68 = t52 * t25;
t32 = (-mrSges(4,1) * t52 + mrSges(4,2) * t49) * qJD(2);
t56 = m(5) * (-qJDD(3) * pkin(3) - t54 * pkin(7) - t68 + (qJD(2) * t33 + t18) * t49) - t19 * mrSges(5,1) + t20 * mrSges(5,2) - t30 * t23 + t31 * t24;
t10 = m(4) * (-t18 * t49 + t68) - t34 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t32 * t66 + qJD(3) * t39 - t56;
t9 = m(4) * t67 - qJDD(3) * mrSges(4,2) + t35 * mrSges(4,3) - qJD(3) * t38 - t48 * t11 + t51 * t12 + t32 * t65;
t4 = m(3) * t63 - t55 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t49 * t10 + t52 * t9;
t6 = m(3) * t25 + t10 * t52 + t49 * t9;
t62 = m(2) * t43 + t4 * t69 + t45 * t71 + t47 * t6;
t2 = m(2) * t37 + t4 * t53 - t50 * t8;
t1 = m(2) * t36 - t45 * t6 + (t4 * t50 + t71) * t47;
t3 = [-m(1) * g(1) - t1 * t44 + t2 * t46, t2, t4, t9, t12; -m(1) * g(2) + t1 * t46 + t2 * t44, t1, t8, t10, t11; -m(1) * g(3) + t62, t62, t6, t72, t56;];
f_new = t3;
