% Calculate vector of inverse dynamics joint torques for
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:23:59
% EndTime: 2019-03-08 18:24:00
% DurationCPUTime: 0.38s
% Computational Cost: add. (384->90), mult. (753->116), div. (0->0), fcn. (454->6), ass. (0->44)
t54 = mrSges(4,1) + mrSges(5,1);
t41 = cos(qJ(2));
t26 = qJD(2) * pkin(2) + t41 * qJD(1);
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t39 = sin(qJ(2));
t51 = qJD(1) * t39;
t15 = t40 * t26 - t38 * t51;
t16 = t38 * t26 + t40 * t51;
t22 = -t38 * t39 + t40 * t41;
t19 = t22 * qJD(1);
t52 = pkin(2) * qJD(3);
t49 = qJD(1) * qJD(2);
t46 = t39 * t49;
t23 = t41 * qJDD(1) - t46;
t20 = qJDD(2) * pkin(2) + t23;
t45 = t41 * t49;
t24 = t39 * qJDD(1) + t45;
t6 = t15 * qJD(3) + t38 * t20 + t40 * t24;
t60 = t6 * t38 * pkin(2) + (t40 * t52 - t19) * t16;
t36 = qJD(2) + qJD(3);
t57 = m(4) + m(5);
t59 = pkin(2) * t57 + mrSges(3,1);
t58 = m(5) * pkin(3);
t11 = t36 * t22;
t21 = t38 * t41 + t40 * t39;
t56 = t16 * t11 + t6 * t21;
t53 = mrSges(4,2) + mrSges(5,2);
t50 = qJD(3) * t38;
t37 = qJ(2) + qJ(3);
t34 = cos(t37);
t47 = t53 * t34;
t33 = sin(t37);
t44 = t33 * t53 - t54 * t34;
t35 = qJDD(2) + qJDD(3);
t7 = -qJD(3) * t16 + t40 * t20 - t38 * t24;
t3 = t35 * pkin(3) + t7;
t43 = t7 * mrSges(4,1) + t3 * mrSges(5,1) - t53 * t6 + (Ifges(4,3) + Ifges(5,3)) * t35;
t42 = qJD(2) ^ 2;
t29 = t40 * pkin(2) + pkin(3);
t18 = t21 * qJD(1);
t14 = t36 * pkin(3) + t15;
t12 = t36 * t21;
t1 = [m(2) * qJDD(1) + (-t39 * qJDD(2) - t42 * t41) * mrSges(3,2) + (t41 * qJDD(2) - t42 * t39) * mrSges(3,1) + m(3) * (t23 * t41 + t24 * t39) + m(4) * (-t15 * t12 + t7 * t22 + t56) + m(5) * (-t14 * t12 + t3 * t22 + t56) + (-t11 * t53 - t12 * t54) * t36 + (-t21 * t53 + t22 * t54) * t35 + (-m(2) - m(3) - t57) * g(2); Ifges(3,3) * qJDD(2) + (-t24 + t45) * mrSges(3,2) + (t23 + t46) * mrSges(3,1) + (t39 * mrSges(3,2) - t34 * t58 - t59 * t41 + t44) * g(2) + (t41 * mrSges(3,2) + t47 + t59 * t39 + (t54 + t58) * t33) * g(1) + (t29 * mrSges(5,1) + (mrSges(4,1) * t40 - t38 * t53) * pkin(2)) * t35 + (t53 * t19 + t54 * t18 + (-t38 * t54 - t40 * t53) * t52) * t36 + t43 + (t3 * t29 + (-pkin(2) * t50 + t18) * t14 + t60) * m(5) + (t15 * t18 + (-t15 * t50 + t40 * t7) * pkin(2) + t60) * m(4); t44 * g(2) + (t54 * t33 + t47) * g(1) + (t35 * mrSges(5,1) + (g(1) * t33 - g(2) * t34 + t3) * m(5)) * pkin(3) + t53 * t15 * t36 + t43 + (-m(5) * (-t14 + t15) + t54 * t36) * t16; (-g(3) + qJDD(4)) * m(5);];
tau  = t1;
