% Calculate vector of cutting forces with Newton-Euler
% S4PRRP6
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP6_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:14
% EndTime: 2019-12-31 16:30:14
% DurationCPUTime: 0.28s
% Computational Cost: add. (1160->87), mult. (2173->116), div. (0->0), fcn. (1054->6), ass. (0->43)
t45 = qJD(2) ^ 2;
t39 = sin(pkin(6));
t54 = cos(pkin(6));
t28 = -t54 * g(1) - t39 * g(2);
t38 = -g(3) + qJDD(1);
t41 = sin(qJ(2));
t43 = cos(qJ(2));
t49 = -t41 * t28 + t43 * t38;
t14 = -qJDD(2) * pkin(2) - t45 * pkin(5) - t49;
t40 = sin(qJ(3));
t42 = cos(qJ(3));
t51 = qJD(2) * qJD(3);
t24 = t40 * qJDD(2) + t42 * t51;
t25 = t42 * qJDD(2) - t40 * t51;
t53 = qJD(2) * t40;
t29 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t53;
t30 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t53;
t52 = qJD(2) * t42;
t32 = mrSges(5,2) * t52 + qJD(3) * mrSges(5,3);
t55 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t52 + t32;
t62 = m(5) * (-t25 * pkin(3) - t24 * qJ(4) + (-0.2e1 * qJD(4) * t40 + (pkin(3) * t40 - qJ(4) * t42) * qJD(3)) * qJD(2) + t14) - t25 * mrSges(5,1);
t66 = ((t29 - t30) * t40 - t55 * t42) * qJD(2) + m(4) * t14 - t25 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t24 + t62;
t57 = t43 * t28 + t41 * t38;
t15 = -t45 * pkin(2) + qJDD(2) * pkin(5) + t57;
t21 = (-pkin(3) * t42 - qJ(4) * t40) * qJD(2);
t44 = qJD(3) ^ 2;
t27 = t39 * g(1) - t54 * g(2);
t61 = t42 * t27;
t63 = m(5) * (-qJDD(3) * pkin(3) - t44 * qJ(4) + t61 + qJDD(4) + (qJD(2) * t21 + t15) * t40);
t59 = mrSges(4,3) + mrSges(5,2);
t58 = t42 * t15 - t40 * t27;
t23 = (-mrSges(4,1) * t42 + mrSges(4,2) * t40) * qJD(2);
t22 = (-mrSges(5,1) * t42 - mrSges(5,3) * t40) * qJD(2);
t48 = m(5) * (-t44 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t21 * t52 + t58) + t22 * t52 + qJD(3) * t30 + qJDD(3) * mrSges(5,3);
t7 = m(4) * t58 - qJDD(3) * mrSges(4,2) - qJD(3) * t29 + t23 * t52 + t59 * t25 + t48;
t8 = m(4) * (-t40 * t15 - t61) - t63 - t59 * t24 + (mrSges(4,1) + mrSges(5,1)) * qJDD(3) + t55 * qJD(3) + (-t22 - t23) * t53;
t3 = m(3) * t57 - t45 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t40 * t8 + t42 * t7;
t6 = m(3) * t49 + qJDD(2) * mrSges(3,1) - t45 * mrSges(3,2) - t66;
t50 = m(2) * t38 + t41 * t3 + t43 * t6;
t47 = -t40 * t7 - t42 * t8;
t4 = (m(2) + m(3)) * t27 + t47;
t1 = m(2) * t28 + t43 * t3 - t41 * t6;
t2 = [-m(1) * g(1) + t54 * t1 - t39 * t4, t1, t3, t7, t25 * mrSges(5,2) + t48; -m(1) * g(2) + t39 * t1 + t54 * t4, t4, t6, t8, -t24 * mrSges(5,3) + (-t40 * t30 - t42 * t32) * qJD(2) + t62; -m(1) * g(3) + t50, t50, -m(3) * t27 - t47, t66, -qJDD(3) * mrSges(5,1) + t24 * mrSges(5,2) - qJD(3) * t32 + t22 * t53 + t63;];
f_new = t2;
