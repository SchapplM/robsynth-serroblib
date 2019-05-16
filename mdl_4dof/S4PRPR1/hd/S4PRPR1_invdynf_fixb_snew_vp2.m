% Calculate vector of cutting forces with Newton-Euler
% S4PRPR1
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
% Datum: 2019-05-04 18:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR1_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR1_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR1_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:56:58
% EndTime: 2019-05-04 18:56:58
% DurationCPUTime: 0.13s
% Computational Cost: add. (936->50), mult. (1433->61), div. (0->0), fcn. (716->6), ass. (0->34)
t45 = -pkin(2) - pkin(3);
t26 = g(3) - qJDD(1);
t44 = m(5) * t26;
t43 = mrSges(3,1) + mrSges(4,1);
t27 = sin(pkin(6));
t28 = cos(pkin(6));
t15 = t27 * g(1) - t28 * g(2);
t16 = -t28 * g(1) - t27 * g(2);
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t42 = t30 * t15 + t32 * t16;
t41 = -m(4) * t26 - t44;
t40 = t32 * t15 - t30 * t16;
t39 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t42;
t38 = -m(3) * t26 + t41;
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t33 = qJD(2) ^ 2;
t34 = -t33 * qJ(3) + qJDD(3) - t40;
t10 = t45 * qJDD(2) + t34;
t23 = -qJD(2) + qJD(4);
t21 = t23 ^ 2;
t22 = -qJDD(2) + qJDD(4);
t8 = t45 * t33 + t39;
t6 = m(5) * (t31 * t10 - t29 * t8) + t22 * mrSges(5,1) - (t21 * mrSges(5,2));
t7 = m(5) * (t29 * t10 + t31 * t8) - t22 * mrSges(5,2) - t21 * mrSges(5,1);
t37 = -t29 * t6 + t31 * t7 + qJDD(2) * mrSges(4,3) + m(4) * (-(t33 * pkin(2)) + t39);
t36 = -m(2) * t26 + t38;
t35 = -m(4) * (-qJDD(2) * pkin(2) + t34) - t29 * t7 - t31 * t6;
t4 = m(3) * t40 + (-mrSges(3,2) + mrSges(4,3)) * t33 + t43 * qJDD(2) + t35;
t3 = m(3) * t42 - qJDD(2) * mrSges(3,2) - t43 * t33 + t37;
t2 = m(2) * t16 + t32 * t3 - t30 * t4;
t1 = m(2) * t15 + t30 * t3 + t32 * t4;
t5 = [-m(1) * g(1) - t27 * t1 + t28 * t2, t2, t3, -(t33 * mrSges(4,1)) + t37, t7; -m(1) * g(2) + t28 * t1 + t27 * t2, t1, t4, t41, t6; -m(1) * g(3) + t36, t36, t38, -qJDD(2) * mrSges(4,1) - t33 * mrSges(4,3) - t35, t44;];
f_new  = t5;
