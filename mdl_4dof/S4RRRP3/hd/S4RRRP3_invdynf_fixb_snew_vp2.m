% Calculate vector of cutting forces with Newton-Euler
% S4RRRP3
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynf_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP3_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP3_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:02
% EndTime: 2019-12-31 17:14:02
% DurationCPUTime: 0.32s
% Computational Cost: add. (2021->95), mult. (2562->123), div. (0->0), fcn. (1149->6), ass. (0->47)
t37 = qJD(1) + qJD(2);
t35 = t37 ^ 2;
t36 = qJDD(1) + qJDD(2);
t41 = sin(qJ(1));
t44 = cos(qJ(1));
t52 = t41 * g(1) - t44 * g(2);
t27 = qJDD(1) * pkin(1) + t52;
t46 = qJD(1) ^ 2;
t48 = -t44 * g(1) - t41 * g(2);
t28 = -t46 * pkin(1) + t48;
t40 = sin(qJ(2));
t43 = cos(qJ(2));
t50 = t43 * t27 - t40 * t28;
t14 = -t36 * pkin(2) - t35 * pkin(6) - t50;
t39 = sin(qJ(3));
t42 = cos(qJ(3));
t53 = qJD(3) * t37;
t21 = t39 * t36 + t42 * t53;
t22 = t42 * t36 - t39 * t53;
t60 = t37 * t39;
t29 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t60;
t30 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t60;
t59 = t37 * t42;
t32 = mrSges(5,2) * t59 + qJD(3) * mrSges(5,3);
t54 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t59 + t32;
t61 = m(5) * (-t22 * pkin(3) - t21 * qJ(4) + (-0.2e1 * qJD(4) * t39 + (pkin(3) * t39 - qJ(4) * t42) * qJD(3)) * t37 + t14) - t22 * mrSges(5,1);
t68 = ((t29 - t30) * t39 - t42 * t54) * t37 + m(4) * t14 - t22 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t21 + t61;
t65 = -m(2) - m(3);
t20 = (-mrSges(4,1) * t42 + mrSges(4,2) * t39) * t37;
t18 = (-pkin(3) * t42 - qJ(4) * t39) * t37;
t19 = (-mrSges(5,1) * t42 - mrSges(5,3) * t39) * t37;
t45 = qJD(3) ^ 2;
t56 = t40 * t27 + t43 * t28;
t15 = -t35 * pkin(2) + t36 * pkin(6) + t56;
t51 = -t39 * g(3) + t42 * t15;
t49 = m(5) * (-t45 * pkin(3) + qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t18 * t59 + t51) + t19 * t59 + qJD(3) * t30 + qJDD(3) * mrSges(5,3);
t57 = mrSges(4,3) + mrSges(5,2);
t7 = m(4) * t51 - qJDD(3) * mrSges(4,2) - qJD(3) * t29 + t20 * t59 + t22 * t57 + t49;
t62 = t42 * g(3);
t63 = m(5) * (-qJDD(3) * pkin(3) + t62 - t45 * qJ(4) + qJDD(4) + (t18 * t37 + t15) * t39);
t8 = m(4) * (-t39 * t15 - t62) - t63 + (-t19 - t20) * t60 - t57 * t21 + (mrSges(4,1) + mrSges(5,1)) * qJDD(3) + t54 * qJD(3);
t64 = t39 * t7 + t42 * t8;
t4 = m(3) * t50 + t36 * mrSges(3,1) - t35 * mrSges(3,2) - t68;
t3 = m(3) * t56 - t35 * mrSges(3,1) - t36 * mrSges(3,2) - t39 * t8 + t42 * t7;
t2 = m(2) * t48 - t46 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t43 * t3 - t40 * t4;
t1 = m(2) * t52 + qJDD(1) * mrSges(2,1) - t46 * mrSges(2,2) + t40 * t3 + t43 * t4;
t5 = [-m(1) * g(1) - t41 * t1 + t44 * t2, t2, t3, t7, t22 * mrSges(5,2) + t49; -m(1) * g(2) + t44 * t1 + t41 * t2, t1, t4, t8, -t21 * mrSges(5,3) + (-t39 * t30 - t42 * t32) * t37 + t61; (-m(1) + t65) * g(3) + t64, g(3) * t65 + t64, -m(3) * g(3) + t64, t68, -qJDD(3) * mrSges(5,1) + t21 * mrSges(5,2) - qJD(3) * t32 + t19 * t60 + t63;];
f_new = t5;
