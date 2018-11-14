% Calculate vector of inverse dynamics joint torques for
% S4PRPP3
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:01
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PRPP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_invdynJ_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:01:19
% EndTime: 2018-11-14 14:01:20
% DurationCPUTime: 0.38s
% Computational Cost: add. (160->74), mult. (305->70), div. (0->0), fcn. (106->2), ass. (0->35)
t21 = cos(qJ(2));
t26 = -t21 * qJD(1) + qJD(3);
t48 = mrSges(4,3) + mrSges(5,2);
t20 = sin(qJ(2));
t33 = t20 * qJD(1);
t13 = qJD(2) * qJ(3) + t33;
t31 = qJD(1) * qJD(2);
t15 = t21 * t31;
t10 = t20 * qJDD(1) + t15;
t46 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t4 = t10 + t46;
t47 = t4 * qJ(3) + t26 * t13;
t45 = -g(1) * t20 + g(2) * t21;
t14 = t20 * t31;
t44 = t14 - t45;
t43 = g(1) * t21 + g(2) * t20 + t15;
t22 = -pkin(2) - pkin(3);
t34 = t13 * qJD(2);
t41 = t4 * t20 + t21 * t34;
t39 = -mrSges(4,1) - mrSges(5,1);
t37 = t21 * pkin(2) + t20 * qJ(3);
t36 = qJD(2) * t20;
t35 = qJDD(2) * pkin(2);
t28 = mrSges(3,1) - t39;
t27 = -mrSges(3,2) + t48;
t9 = t21 * qJDD(1) - t14;
t25 = qJDD(3) - t9;
t24 = -t34 + t45;
t23 = qJD(2) ^ 2;
t18 = t21 * qJ(3);
t12 = -qJD(2) * pkin(2) + t26;
t8 = t22 * qJD(2) + t26;
t5 = t25 - t35;
t3 = t22 * qJDD(2) + t25;
t1 = [m(2) * qJDD(1) + m(3) * (t10 * t20 + t9 * t21) + m(4) * (t12 * t36 - t5 * t21 + t41) + m(5) * (-t3 * t21 + t36 * t8 + t41) + (-m(2) - m(3) - m(4) - m(5)) * g(2) + (-t20 * t28 + t21 * t27) * t23 + (t20 * t27 + t21 * t28) * qJDD(2); (-t10 + t43) * mrSges(3,2) + (t44 + t9) * mrSges(3,1) + (-t3 + t44) * mrSges(5,1) + (t3 * t22 - g(1) * (t22 * t20 + t18) - t33 * t8 - g(2) * (t21 * pkin(3) + t37) + t47) * m(5) + (-t5 * pkin(2) - t12 * t33 - g(2) * t37 - g(1) * (-t20 * pkin(2) + t18) + t47) * m(4) + (t35 + t44 - t5) * mrSges(4,1) + (-t22 * mrSges(5,1) + Ifges(4,2) + Ifges(3,3) + Ifges(5,3)) * qJDD(2) + t48 * (t4 - t43 + t46); -t48 * t23 + t39 * qJDD(2) + (t24 + t3) * m(5) + (t24 + t5) * m(4); (g(3) + qJDD(4)) * m(5);];
tau  = t1;
