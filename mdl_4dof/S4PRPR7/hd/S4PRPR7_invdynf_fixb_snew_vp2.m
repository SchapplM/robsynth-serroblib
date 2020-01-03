% Calculate vector of cutting forces with Newton-Euler
% S4PRPR7
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR7_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:37
% EndTime: 2019-12-31 16:25:38
% DurationCPUTime: 0.18s
% Computational Cost: add. (790->63), mult. (1320->82), div. (0->0), fcn. (632->6), ass. (0->38)
t53 = -pkin(2) - pkin(5);
t52 = mrSges(3,1) - mrSges(4,2);
t51 = -mrSges(3,2) + mrSges(4,3);
t32 = sin(pkin(6));
t49 = cos(pkin(6));
t25 = -t49 * g(1) - t32 * g(2);
t31 = -g(3) + qJDD(1);
t34 = sin(qJ(2));
t36 = cos(qJ(2));
t50 = t36 * t25 + t34 * t31;
t33 = sin(qJ(4));
t48 = qJD(2) * t33;
t35 = cos(qJ(4));
t47 = qJD(2) * t35;
t46 = qJD(2) * qJD(4);
t37 = qJD(2) ^ 2;
t44 = -t34 * t25 + t36 * t31;
t40 = -t37 * qJ(3) + qJDD(3) - t44;
t12 = t53 * qJDD(2) + t40;
t20 = (mrSges(5,1) * t33 + mrSges(5,2) * t35) * qJD(2);
t22 = t35 * qJDD(2) - t33 * t46;
t24 = t32 * g(1) - t49 * g(2);
t26 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t48;
t8 = m(5) * (t35 * t12 + t33 * t24) - t22 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t20 * t47 + qJD(4) * t26;
t21 = -t33 * qJDD(2) - t35 * t46;
t27 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t47;
t9 = m(5) * (t33 * t12 - t35 * t24) + t21 * mrSges(5,3) - qJDD(4) * mrSges(5,2) - t20 * t48 - qJD(4) * t27;
t41 = -m(4) * (-qJDD(2) * pkin(2) + t40) - t33 * t9 - t35 * t8;
t3 = m(3) * t44 + t52 * qJDD(2) + t51 * t37 + t41;
t39 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t50;
t42 = -t21 * mrSges(5,1) + m(5) * (t53 * t37 + t39) + t26 * t48 + t27 * t47 + t22 * mrSges(5,2);
t38 = -m(4) * (t37 * pkin(2) - t39) + t42;
t6 = m(3) * t50 + t51 * qJDD(2) - t52 * t37 + t38;
t45 = m(2) * t31 + t36 * t3 + t34 * t6;
t43 = m(4) * t24 + t33 * t8 - t35 * t9;
t4 = (m(2) + m(3)) * t24 + t43;
t1 = m(2) * t25 - t34 * t3 + t36 * t6;
t2 = [-m(1) * g(1) + t49 * t1 - t32 * t4, t1, t6, -t43, t9; -m(1) * g(2) + t32 * t1 + t49 * t4, t4, t3, -t37 * mrSges(4,2) - qJDD(2) * mrSges(4,3) - t38, t8; -m(1) * g(3) + t45, t45, -m(3) * t24 - t43, qJDD(2) * mrSges(4,2) - t37 * mrSges(4,3) - t41, t42;];
f_new = t2;
