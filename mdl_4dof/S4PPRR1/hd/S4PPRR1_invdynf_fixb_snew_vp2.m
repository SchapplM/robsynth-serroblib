% Calculate vector of cutting forces with Newton-Euler
% S4PPRR1
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-05-04 18:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PPRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR1_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR1_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR1_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:47:36
% EndTime: 2019-05-04 18:47:36
% DurationCPUTime: 0.13s
% Computational Cost: add. (782->38), mult. (1093->51), div. (0->0), fcn. (724->6), ass. (0->31)
t23 = sin(pkin(6));
t24 = cos(pkin(6));
t32 = t23 * g(1) - t24 * g(2);
t13 = qJDD(2) - t32;
t15 = -t24 * g(1) - t23 * g(2);
t26 = sin(qJ(3));
t28 = cos(qJ(3));
t37 = t26 * t13 + t28 * t15;
t22 = g(3) - qJDD(1);
t16 = m(5) * t22;
t36 = m(4) * t22 + t16;
t35 = t28 * t13 - t26 * t15;
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t29 = qJD(3) ^ 2;
t21 = qJD(3) + qJD(4);
t19 = t21 ^ 2;
t20 = qJDD(3) + qJDD(4);
t8 = qJDD(3) * pkin(3) + t35;
t9 = -t29 * pkin(3) + t37;
t6 = m(5) * (-t25 * t9 + t27 * t8) + t20 * mrSges(5,1) - t19 * mrSges(5,2);
t7 = m(5) * (t25 * t8 + t27 * t9) - t20 * mrSges(5,2) - t19 * mrSges(5,1);
t4 = m(4) * t35 + qJDD(3) * mrSges(4,1) - t29 * mrSges(4,2) + t25 * t7 + t27 * t6;
t5 = m(4) * t37 - t29 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t25 * t6 + t27 * t7;
t34 = m(3) * t15 - t26 * t4 + t28 * t5;
t33 = -m(3) * t22 - t36;
t31 = -m(2) * t22 + t33;
t30 = m(3) * t13 + t26 * t5 + t28 * t4;
t2 = m(2) * t15 + t34;
t1 = m(2) * t32 - t30;
t3 = [-m(1) * g(1) - t23 * t1 + t24 * t2, t2, t34, t5, t7; -m(1) * g(2) + t24 * t1 + t23 * t2, t1, t33, t4, t6; -m(1) * g(3) + t31, t31, t30, t36, t16;];
f_new  = t3;
