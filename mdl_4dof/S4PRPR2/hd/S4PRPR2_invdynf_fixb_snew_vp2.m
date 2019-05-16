% Calculate vector of cutting forces with Newton-Euler
% S4PRPR2
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
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2019-05-04 19:00
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4PRPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR2_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR2_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR2_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR2_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:59:47
% EndTime: 2019-05-04 18:59:47
% DurationCPUTime: 0.13s
% Computational Cost: add. (1299->47), mult. (1681->57), div. (0->0), fcn. (868->6), ass. (0->32)
t42 = m(2) + m(3);
t27 = -g(2) + qJDD(1);
t31 = sin(qJ(2));
t33 = cos(qJ(2));
t39 = t31 * g(1) + t33 * t27;
t15 = qJDD(2) * pkin(2) + t39;
t34 = qJD(2) ^ 2;
t36 = -t33 * g(1) + t31 * t27;
t16 = -t34 * pkin(2) + t36;
t28 = sin(pkin(6));
t29 = cos(pkin(6));
t41 = t28 * t15 + t29 * t16;
t26 = -g(3) + qJDD(3);
t19 = m(5) * t26;
t40 = m(4) * t26 + t19;
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t35 = t29 * t15 - t28 * t16;
t10 = qJDD(2) * pkin(3) + t35;
t11 = -t34 * pkin(3) + t41;
t25 = qJD(2) + qJD(4);
t23 = t25 ^ 2;
t24 = qJDD(2) + qJDD(4);
t8 = m(5) * (t32 * t10 - t30 * t11) + t24 * mrSges(5,1) - t23 * mrSges(5,2);
t9 = m(5) * (t30 * t10 + t32 * t11) - t24 * mrSges(5,2) - t23 * mrSges(5,1);
t6 = m(4) * t35 + qJDD(2) * mrSges(4,1) - t34 * mrSges(4,2) + t30 * t9 + t32 * t8;
t7 = m(4) * t41 - t34 * mrSges(4,1) - qJDD(2) * mrSges(4,2) - t30 * t8 + t32 * t9;
t4 = m(3) * t39 + qJDD(2) * mrSges(3,1) - t34 * mrSges(3,2) + t28 * t7 + t29 * t6;
t5 = m(3) * t36 - t34 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t28 * t6 + t29 * t7;
t38 = m(2) * t27 + t31 * t5 + t33 * t4;
t37 = -t31 * t4 + t33 * t5;
t1 = [(-m(1) - m(2)) * g(1) + t37, -m(2) * g(1) + t37, t5, t7, t9; -m(1) * g(2) + t38, t42 * g(3) - t40, t4, t6, t8; (-m(1) - t42) * g(3) + t40, t38, -m(3) * g(3) + t40, t40, t19;];
f_new  = t1;
