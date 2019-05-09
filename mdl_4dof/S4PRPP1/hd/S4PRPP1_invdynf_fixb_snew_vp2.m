% Calculate vector of cutting forces with Newton-Euler
% S4PRPP1
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
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2019-05-04 18:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRPP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP1_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_invdynf_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP1_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP1_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:53:21
% EndTime: 2019-05-04 18:53:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (454->49), mult. (712->54), div. (0->0), fcn. (316->4), ass. (0->29)
t26 = qJD(2) ^ 2;
t22 = sin(pkin(5));
t23 = cos(pkin(5));
t14 = t22 * g(1) - t23 * g(2);
t15 = -t23 * g(1) - t22 * g(2);
t24 = sin(qJ(2));
t25 = cos(qJ(2));
t36 = t24 * t14 + t25 * t15;
t27 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t36;
t37 = -pkin(2) - qJ(4);
t40 = qJDD(2) * mrSges(5,2) + m(5) * (t37 * t26 + qJDD(4) + t27);
t39 = mrSges(4,2) - mrSges(5,3);
t38 = -mrSges(5,2) - mrSges(4,3);
t21 = -g(3) + qJDD(1);
t16 = m(5) * t21;
t35 = m(4) * t21 + t16;
t34 = mrSges(3,1) - t39;
t33 = m(3) * t21 + t35;
t31 = t25 * t14 - t24 * t15;
t28 = -t26 * qJ(3) + qJDD(3) - t31;
t5 = m(5) * (-(2 * qJD(4) * qJD(2)) + t37 * qJDD(2) + t28);
t32 = m(4) * (-qJDD(2) * pkin(2) + t28) + t5;
t30 = m(2) * t21 + t33;
t29 = m(4) * (t26 * pkin(2) - t27) - t40;
t4 = m(3) * t31 + (-mrSges(3,2) - t38) * t26 + t34 * qJDD(2) - t32;
t3 = m(3) * t36 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t34 * t26 - t29;
t2 = m(2) * t15 - t24 * t4 + t25 * t3;
t1 = m(2) * t14 + t24 * t3 + t25 * t4;
t6 = [-m(1) * g(1) - t22 * t1 + t23 * t2, t2, t3, t35, t16; -m(1) * g(2) + t23 * t1 + t22 * t2, t1, t4, -qJDD(2) * mrSges(4,3) - t39 * t26 + t29, -t26 * mrSges(5,2) - qJDD(2) * mrSges(5,3) + t5; -m(1) * g(3) + t30, t30, t33, t39 * qJDD(2) + t38 * t26 + t32, -t26 * mrSges(5,3) + t40;];
f_new  = t6;
