% Calculate vector of cutting forces with Newton-Euler
% S4RPRR5
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR5_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR5_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR5_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:33
% EndTime: 2019-12-31 16:51:34
% DurationCPUTime: 0.23s
% Computational Cost: add. (1525->73), mult. (2006->92), div. (0->0), fcn. (714->6), ass. (0->42)
t24 = -qJDD(1) + qJDD(3);
t28 = sin(qJ(4));
t31 = cos(qJ(4));
t25 = -qJD(1) + qJD(3);
t46 = qJD(4) * t25;
t17 = t24 * t28 + t31 * t46;
t18 = t24 * t31 - t28 * t46;
t50 = t25 * t28;
t19 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t50;
t49 = t25 * t31;
t20 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t49;
t23 = t25 ^ 2;
t34 = qJD(1) ^ 2;
t30 = sin(qJ(1));
t33 = cos(qJ(1));
t43 = -t33 * g(1) - t30 * g(2);
t38 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t43;
t51 = -pkin(1) - pkin(2);
t12 = t34 * t51 + t38;
t44 = t30 * g(1) - t33 * g(2);
t35 = -t34 * qJ(2) + qJDD(2) - t44;
t14 = qJDD(1) * t51 + t35;
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t41 = -t12 * t29 + t14 * t32;
t53 = (t19 * t28 - t20 * t31) * t25 + m(5) * (-pkin(3) * t24 - pkin(6) * t23 - t41) - t18 * mrSges(5,1) + t17 * mrSges(5,2);
t52 = -m(3) - m(4);
t48 = mrSges(2,1) + mrSges(3,1);
t47 = t32 * t12 + t29 * t14;
t45 = -m(2) + t52;
t16 = (-mrSges(5,1) * t31 + mrSges(5,2) * t28) * t25;
t9 = -pkin(3) * t23 + pkin(6) * t24 + t47;
t6 = m(5) * (g(3) * t31 - t28 * t9) - t17 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t16 * t50 + qJD(4) * t20;
t7 = m(5) * (g(3) * t28 + t31 * t9) + t18 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t16 * t49 - qJD(4) * t19;
t42 = -t28 * t7 - t31 * t6;
t4 = m(4) * t47 - t23 * mrSges(4,1) - t24 * mrSges(4,2) - t28 * t6 + t31 * t7;
t5 = m(4) * t41 + t24 * mrSges(4,1) - t23 * mrSges(4,2) - t53;
t39 = -t29 * t5 + t32 * t4 + m(3) * (-pkin(1) * t34 + t38) + qJDD(1) * mrSges(3,3);
t37 = -m(3) * (-qJDD(1) * pkin(1) + t35) - t29 * t4 - t32 * t5;
t2 = m(2) * t44 + (-mrSges(2,2) + mrSges(3,3)) * t34 + t48 * qJDD(1) + t37;
t1 = m(2) * t43 - qJDD(1) * mrSges(2,2) - t34 * t48 + t39;
t3 = [-m(1) * g(1) + t1 * t33 - t2 * t30, t1, -t34 * mrSges(3,1) + t39, t4, t7; -m(1) * g(2) + t1 * t30 + t2 * t33, t2, g(3) * t52 + t42, t5, t6; (-m(1) + t45) * g(3) + t42, g(3) * t45 + t42, -qJDD(1) * mrSges(3,1) - t34 * mrSges(3,3) - t37, m(4) * g(3) - t42, t53;];
f_new = t3;
