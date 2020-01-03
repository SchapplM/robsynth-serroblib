% Calculate vector of cutting forces with Newton-Euler
% S4RPRR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S4RPRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynf_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynf_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynf_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_invdynf_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_invdynf_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR4_invdynf_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:15
% EndTime: 2019-12-31 16:50:17
% DurationCPUTime: 0.38s
% Computational Cost: add. (2882->96), mult. (5416->130), div. (0->0), fcn. (2978->8), ass. (0->55)
t50 = qJD(1) ^ 2;
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t58 = t45 * g(1) - t48 * g(2);
t28 = qJDD(1) * pkin(1) + t58;
t54 = -t48 * g(1) - t45 * g(2);
t30 = -t50 * pkin(1) + t54;
t41 = sin(pkin(7));
t42 = cos(pkin(7));
t55 = t42 * t28 - t41 * t30;
t15 = -qJDD(1) * pkin(2) - t50 * pkin(5) - t55;
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t60 = qJD(1) * qJD(3);
t56 = t47 * t60;
t32 = t44 * qJDD(1) + t56;
t57 = t44 * t60;
t33 = t47 * qJDD(1) - t57;
t11 = (-t32 - t56) * pkin(6) + (-t33 + t57) * pkin(3) + t15;
t31 = (-pkin(3) * t47 - pkin(6) * t44) * qJD(1);
t49 = qJD(3) ^ 2;
t61 = t47 * qJD(1);
t63 = t41 * t28 + t42 * t30;
t16 = -t50 * pkin(2) + qJDD(1) * pkin(5) + t63;
t40 = -g(3) + qJDD(2);
t64 = t47 * t16 + t44 * t40;
t13 = -t49 * pkin(3) + qJDD(3) * pkin(6) + t31 * t61 + t64;
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t62 = qJD(1) * t44;
t27 = t43 * qJD(3) + t46 * t62;
t17 = -t27 * qJD(4) + t46 * qJDD(3) - t43 * t32;
t26 = t46 * qJD(3) - t43 * t62;
t19 = -t26 * mrSges(5,1) + t27 * mrSges(5,2);
t36 = qJD(4) - t61;
t21 = t36 * mrSges(5,1) - t27 * mrSges(5,3);
t25 = qJDD(4) - t33;
t10 = m(5) * (t43 * t11 + t46 * t13) + t17 * mrSges(5,3) - t25 * mrSges(5,2) + t26 * t19 - t36 * t21;
t34 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t62;
t35 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t61;
t18 = t26 * qJD(4) + t43 * qJDD(3) + t46 * t32;
t20 = -t36 * mrSges(5,2) + t26 * mrSges(5,3);
t9 = m(5) * (t46 * t11 - t43 * t13) - t18 * mrSges(5,3) + t25 * mrSges(5,1) - t27 * t19 + t36 * t20;
t66 = m(4) * t15 - t33 * mrSges(4,1) + t32 * mrSges(4,2) + t43 * t10 + t46 * t9 + (t44 * t34 - t47 * t35) * qJD(1);
t65 = t47 * t40;
t29 = (-mrSges(4,1) * t47 + mrSges(4,2) * t44) * qJD(1);
t6 = m(4) * t64 - qJDD(3) * mrSges(4,2) + t33 * mrSges(4,3) - qJD(3) * t34 + t46 * t10 + t29 * t61 - t43 * t9;
t51 = m(5) * (-qJDD(3) * pkin(3) - t49 * pkin(6) - t65 + (qJD(1) * t31 + t16) * t44) - t17 * mrSges(5,1) + t18 * mrSges(5,2) - t26 * t20 + t27 * t21;
t8 = m(4) * (-t44 * t16 + t65) - t32 * mrSges(4,3) + qJDD(3) * mrSges(4,1) - t29 * t62 + qJD(3) * t35 - t51;
t59 = m(3) * t40 + t44 * t6 + t47 * t8;
t4 = m(3) * t55 + qJDD(1) * mrSges(3,1) - t50 * mrSges(3,2) - t66;
t3 = m(3) * t63 - t50 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t44 * t8 + t47 * t6;
t2 = m(2) * t54 - t50 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t42 * t3 - t41 * t4;
t1 = m(2) * t58 + qJDD(1) * mrSges(2,1) - t50 * mrSges(2,2) + t41 * t3 + t42 * t4;
t5 = [-m(1) * g(1) - t45 * t1 + t48 * t2, t2, t3, t6, t10; -m(1) * g(2) + t48 * t1 + t45 * t2, t1, t4, t8, t9; (-m(1) - m(2)) * g(3) + t59, -m(2) * g(3) + t59, t59, t66, t51;];
f_new = t5;
