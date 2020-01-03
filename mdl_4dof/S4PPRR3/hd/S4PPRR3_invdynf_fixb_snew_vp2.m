% Calculate vector of cutting forces with Newton-Euler
% S4PPRR3
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PPRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR3_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:22
% EndTime: 2019-12-31 16:17:23
% DurationCPUTime: 0.18s
% Computational Cost: add. (762->51), mult. (1316->75), div. (0->0), fcn. (714->6), ass. (0->35)
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t39 = qJD(3) * qJD(4);
t13 = t25 * qJDD(3) + t27 * t39;
t14 = t27 * qJDD(3) - t25 * t39;
t41 = qJD(3) * t25;
t18 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t41;
t40 = qJD(3) * t27;
t19 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t40;
t29 = qJD(3) ^ 2;
t23 = sin(pkin(6));
t24 = cos(pkin(6));
t37 = t23 * g(1) - t24 * g(2);
t15 = qJDD(2) - t37;
t17 = -t24 * g(1) - t23 * g(2);
t26 = sin(qJ(3));
t28 = cos(qJ(3));
t36 = t28 * t15 - t26 * t17;
t43 = (t25 * t18 - t27 * t19) * qJD(3) + m(5) * (-qJDD(3) * pkin(3) - t29 * pkin(5) - t36) - t14 * mrSges(5,1) + t13 * mrSges(5,2);
t42 = t26 * t15 + t28 * t17;
t12 = (-mrSges(5,1) * t27 + mrSges(5,2) * t25) * qJD(3);
t22 = g(3) - qJDD(1);
t9 = -t29 * pkin(3) + qJDD(3) * pkin(5) + t42;
t6 = m(5) * (t27 * t22 - t25 * t9) - t13 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t12 * t41 + qJD(4) * t19;
t7 = m(5) * (t25 * t22 + t27 * t9) + t14 * mrSges(5,3) - qJDD(4) * mrSges(5,2) + t12 * t40 - qJD(4) * t18;
t4 = m(4) * t42 - t29 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t25 * t6 + t27 * t7;
t5 = m(4) * t36 + qJDD(3) * mrSges(4,1) - t29 * mrSges(4,2) - t43;
t38 = m(3) * t17 - t26 * t5 + t28 * t4;
t34 = m(3) * t15 + t26 * t4 + t28 * t5;
t33 = m(4) * t22 + t25 * t7 + t27 * t6;
t31 = -m(3) * t22 - t33;
t30 = -m(2) * t22 + t31;
t2 = m(2) * t17 + t38;
t1 = m(2) * t37 - t34;
t3 = [-m(1) * g(1) - t23 * t1 + t24 * t2, t2, t38, t4, t7; -m(1) * g(2) + t24 * t1 + t23 * t2, t1, t31, t5, t6; -m(1) * g(3) + t30, t30, t34, t33, t43;];
f_new = t3;
