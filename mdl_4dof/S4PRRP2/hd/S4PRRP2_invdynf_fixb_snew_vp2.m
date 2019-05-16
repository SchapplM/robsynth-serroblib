% Calculate vector of cutting forces with Newton-Euler
% S4PRRP2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2019-05-04 19:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP2_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP2_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP2_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:03:23
% EndTime: 2019-05-04 19:03:23
% DurationCPUTime: 0.10s
% Computational Cost: add. (663->46), mult. (743->49), div. (0->0), fcn. (324->4), ass. (0->28)
t40 = -m(3) - m(4);
t22 = qJDD(2) + qJDD(3);
t25 = -g(2) + qJDD(1);
t27 = sin(qJ(2));
t29 = cos(qJ(2));
t36 = t27 * g(1) + t29 * t25;
t13 = qJDD(2) * pkin(2) + t36;
t30 = qJD(2) ^ 2;
t32 = -t29 * g(1) + t27 * t25;
t14 = -t30 * pkin(2) + t32;
t26 = sin(qJ(3));
t28 = cos(qJ(3));
t31 = t28 * t13 - t26 * t14;
t39 = t22 * mrSges(5,1) + m(5) * (t22 * pkin(3) + t31);
t38 = -mrSges(4,2) - mrSges(5,2);
t37 = t26 * t13 + t28 * t14;
t35 = m(2) - t40;
t23 = qJD(2) + qJD(3);
t21 = t23 ^ 2;
t6 = m(4) * t31 + t22 * mrSges(4,1) + t38 * t21 + t39;
t9 = m(5) * (-t21 * pkin(3) + t37);
t7 = m(4) * t37 + t9 + t38 * t22 + (-mrSges(4,1) - mrSges(5,1)) * t21;
t4 = m(3) * t36 + qJDD(2) * mrSges(3,1) - t30 * mrSges(3,2) + t26 * t7 + t28 * t6;
t5 = m(3) * t32 - t30 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t26 * t6 + t28 * t7;
t34 = m(2) * t25 + t27 * t5 + t29 * t4;
t33 = -t27 * t4 + t29 * t5;
t18 = m(5) * (-g(3) + qJDD(4));
t1 = [(-m(1) - m(2)) * g(1) + t33, -m(2) * g(1) + t33, t5, t7, -t21 * mrSges(5,1) - t22 * mrSges(5,2) + t9; -m(1) * g(2) + t34, t35 * g(3) - t18, t4, t6, -t21 * mrSges(5,2) + t39; t18 + (-m(1) - t35) * g(3), t34, t40 * g(3) + t18, -m(4) * g(3) + t18, t18;];
f_new  = t1;
