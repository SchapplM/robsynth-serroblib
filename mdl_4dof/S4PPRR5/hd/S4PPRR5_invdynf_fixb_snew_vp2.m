% Calculate vector of cutting forces with Newton-Euler
% S4PPRR5
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PPRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR5_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR5_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:42
% EndTime: 2019-12-31 16:19:42
% DurationCPUTime: 0.17s
% Computational Cost: add. (735->52), mult. (1227->75), div. (0->0), fcn. (646->6), ass. (0->34)
t27 = sin(qJ(4));
t29 = cos(qJ(4));
t40 = qJD(3) * qJD(4);
t14 = t27 * qJDD(3) + t29 * t40;
t15 = t29 * qJDD(3) - t27 * t40;
t42 = qJD(3) * t27;
t19 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t42;
t41 = qJD(3) * t29;
t20 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t41;
t31 = qJD(3) ^ 2;
t25 = sin(pkin(6));
t26 = cos(pkin(6));
t37 = t25 * g(1) - t26 * g(2);
t16 = qJDD(2) - t37;
t24 = -g(3) + qJDD(1);
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t36 = t30 * t16 - t28 * t24;
t44 = (t27 * t19 - t29 * t20) * qJD(3) + m(5) * (-qJDD(3) * pkin(3) - t31 * pkin(5) - t36) - t15 * mrSges(5,1) + t14 * mrSges(5,2);
t43 = t28 * t16 + t30 * t24;
t18 = t26 * g(1) + t25 * g(2);
t11 = -t31 * pkin(3) + qJDD(3) * pkin(5) + t43;
t13 = (-mrSges(5,1) * t29 + mrSges(5,2) * t27) * qJD(3);
t8 = m(5) * (-t27 * t11 - t29 * t18) - t14 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t13 * t42 + qJD(4) * t20;
t9 = m(5) * (t29 * t11 - t27 * t18) + t15 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t13 * t41 - qJD(4) * t19;
t39 = -m(4) * t18 + t27 * t9 + t29 * t8;
t4 = m(4) * t43 - t31 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t27 * t8 + t29 * t9;
t5 = m(4) * t36 + qJDD(3) * mrSges(4,1) - t31 * mrSges(4,2) - t44;
t38 = m(3) * t24 - t28 * t5 + t30 * t4;
t34 = m(2) * t24 + t38;
t33 = m(3) * t16 + t28 * t4 + t30 * t5;
t2 = (-m(2) - m(3)) * t18 + t39;
t1 = m(2) * t37 - t33;
t3 = [-m(1) * g(1) - t25 * t1 + t26 * t2, t2, t38, t4, t9; -m(1) * g(2) + t26 * t1 + t25 * t2, t1, m(3) * t18 - t39, t5, t8; -m(1) * g(3) + t34, t34, t33, t39, t44;];
f_new = t3;
