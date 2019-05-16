% Calculate vector of cutting forces with Newton-Euler
% S4RPPR2
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-05-04 19:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR2_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR2_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR2_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:10:44
% EndTime: 2019-05-04 19:10:44
% DurationCPUTime: 0.18s
% Computational Cost: add. (1483->58), mult. (2201->67), div. (0->0), fcn. (724->6), ass. (0->35)
t44 = -m(2) - m(3);
t43 = -pkin(1) - pkin(2);
t42 = mrSges(2,1) + mrSges(3,1);
t32 = qJD(1) ^ 2;
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t37 = -t31 * g(1) - t29 * g(2);
t35 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t37;
t13 = t43 * t32 + t35;
t39 = t29 * g(1) - t31 * g(2);
t33 = -t32 * qJ(2) + qJDD(2) - t39;
t15 = t43 * qJDD(1) + t33;
t26 = sin(pkin(6));
t27 = cos(pkin(6));
t41 = t27 * t13 + t26 * t15;
t25 = g(3) + qJDD(3);
t17 = m(5) * t25;
t40 = m(4) * t25 + t17;
t38 = -t26 * t13 + t27 * t15;
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t22 = -qJD(1) + qJD(4);
t20 = t22 ^ 2;
t21 = -qJDD(1) + qJDD(4);
t8 = -qJDD(1) * pkin(3) + t38;
t9 = -t32 * pkin(3) + t41;
t6 = m(5) * (-t28 * t9 + t30 * t8) + t21 * mrSges(5,1) - (t20 * mrSges(5,2));
t7 = m(5) * (t28 * t8 + t30 * t9) - t21 * mrSges(5,2) - t20 * mrSges(5,1);
t4 = m(4) * t38 - qJDD(1) * mrSges(4,1) - t32 * mrSges(4,2) + t28 * t7 + t30 * t6;
t5 = m(4) * t41 - t32 * mrSges(4,1) + qJDD(1) * mrSges(4,2) - t28 * t6 + t30 * t7;
t36 = -t26 * t4 + t27 * t5 + m(3) * (-t32 * pkin(1) + t35) + qJDD(1) * mrSges(3,3);
t34 = -m(3) * (-qJDD(1) * pkin(1) + t33) - t26 * t5 - t27 * t4;
t2 = m(2) * t39 + (-mrSges(2,2) + mrSges(3,3)) * t32 + t42 * qJDD(1) + t34;
t1 = m(2) * t37 - qJDD(1) * mrSges(2,2) - t42 * t32 + t36;
t3 = [-m(1) * g(1) + t31 * t1 - t29 * t2, t1, -t32 * mrSges(3,1) + t36, t5, t7; -m(1) * g(2) + t29 * t1 + t31 * t2, t2, -m(3) * g(3) - t40, t4, t6; (-m(1) + t44) * g(3) - t40, t44 * g(3) - t40, -qJDD(1) * mrSges(3,1) - t32 * mrSges(3,3) - t34, t40, t17;];
f_new  = t3;
