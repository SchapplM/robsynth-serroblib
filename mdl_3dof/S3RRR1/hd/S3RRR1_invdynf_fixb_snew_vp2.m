% Calculate vector of cutting forces with Newton-Euler
% S3RRR1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x4]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 18:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S3RRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_invdynf_fixb_snew_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_invdynf_fixb_snew_vp2: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRR1_invdynf_fixb_snew_vp2: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_invdynf_fixb_snew_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_invdynf_fixb_snew_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRR1_invdynf_fixb_snew_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:32:31
% EndTime: 2019-05-04 18:32:32
% DurationCPUTime: 0.11s
% Computational Cost: add. (825->42), mult. (1084->55), div. (0->0), fcn. (560->6), ass. (0->30)
t33 = -m(3) - m(4);
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t30 = t23 * g(1) - t26 * g(2);
t12 = qJDD(1) * pkin(1) + t30;
t27 = qJD(1) ^ 2;
t28 = -t26 * g(1) - t23 * g(2);
t13 = -t27 * pkin(1) + t28;
t22 = sin(qJ(2));
t25 = cos(qJ(2));
t32 = t22 * t12 + t25 * t13;
t31 = -m(2) + t33;
t20 = qJD(1) + qJD(2);
t19 = qJDD(1) + qJDD(2);
t29 = t25 * t12 - t22 * t13;
t24 = cos(qJ(3));
t21 = sin(qJ(3));
t18 = t20 ^ 2;
t16 = qJD(3) + t20;
t15 = qJDD(3) + t19;
t14 = t16 ^ 2;
t8 = -t18 * pkin(2) + t32;
t7 = t19 * pkin(2) + t29;
t6 = m(4) * (t21 * t7 + t24 * t8) - t15 * mrSges(4,2) - t14 * mrSges(4,1);
t5 = m(4) * (-t21 * t8 + t24 * t7) + t15 * mrSges(4,1) - t14 * mrSges(4,2);
t4 = m(3) * t32 - t18 * mrSges(3,1) - t19 * mrSges(3,2) - t21 * t5 + t24 * t6;
t3 = m(3) * t29 + t19 * mrSges(3,1) - t18 * mrSges(3,2) + t21 * t6 + t24 * t5;
t2 = m(2) * t28 - t27 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t22 * t3 + t25 * t4;
t1 = m(2) * t30 + qJDD(1) * mrSges(2,1) - t27 * mrSges(2,2) + t22 * t4 + t25 * t3;
t9 = [-m(1) * g(1) - t23 * t1 + t26 * t2, t2, t4, t6; -m(1) * g(2) + t26 * t1 + t23 * t2, t1, t3, t5; (-m(1) + t31) * g(3), t31 * g(3), t33 * g(3), -m(4) * g(3);];
f_new  = t9;
