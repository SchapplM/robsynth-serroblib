% Calculate vector of inverse dynamics joint torques for
% S4PRPR3
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:11:25
% EndTime: 2018-11-14 14:11:25
% DurationCPUTime: 0.38s
% Computational Cost: add. (447->106), mult. (921->142), div. (0->0), fcn. (646->10), ass. (0->57)
t64 = m(4) + m(5);
t47 = cos(pkin(6));
t37 = t47 * pkin(2) + pkin(3);
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t46 = sin(pkin(6));
t63 = pkin(2) * t46;
t21 = t50 * t37 - t48 * t63;
t49 = sin(qJ(2));
t51 = cos(qJ(2));
t29 = t46 * t51 + t47 * t49;
t23 = t29 * qJD(1);
t28 = -t46 * t49 + t47 * t51;
t25 = t28 * qJD(1);
t67 = t21 * qJD(4) + t48 * t23 - t50 * t25;
t22 = t48 * t37 + t50 * t63;
t66 = -t22 * qJD(4) + t50 * t23 + t48 * t25;
t65 = pkin(2) * t64 + mrSges(3,1);
t45 = qJ(2) + pkin(6);
t42 = qJ(4) + t45;
t35 = sin(t42);
t62 = g(2) * t35;
t56 = qJD(1) * qJD(2);
t31 = t51 * qJDD(1) - t49 * t56;
t27 = qJDD(2) * pkin(2) + t31;
t32 = t49 * qJDD(1) + t51 * t56;
t11 = t47 * t27 - t46 * t32;
t10 = qJDD(2) * pkin(3) + t11;
t12 = t46 * t27 + t47 * t32;
t33 = qJD(2) * pkin(2) + t51 * qJD(1);
t58 = qJD(1) * t49;
t16 = t47 * t33 - t46 * t58;
t15 = qJD(2) * pkin(3) + t16;
t17 = t46 * t33 + t47 * t58;
t7 = t48 * t15 + t50 * t17;
t3 = -qJD(4) * t7 + t50 * t10 - t48 * t12;
t43 = qJDD(2) + qJDD(4);
t61 = t3 * mrSges(5,1) + Ifges(5,3) * t43;
t60 = t49 * mrSges(3,2);
t59 = t51 * mrSges(3,2);
t55 = m(5) * pkin(3) + mrSges(4,1);
t54 = t49 * mrSges(3,1) + t59;
t6 = t50 * t15 - t48 * t17;
t13 = t50 * t28 - t48 * t29;
t14 = t48 * t28 + t50 * t29;
t2 = qJD(4) * t6 + t48 * t10 + t50 * t12;
t36 = cos(t42);
t53 = g(1) * t35 + g(2) * t36 - t2;
t44 = qJD(2) + qJD(4);
t40 = cos(t45);
t39 = sin(t45);
t34 = t36 * mrSges(5,1);
t26 = t28 * qJD(2);
t24 = t29 * qJD(2);
t5 = -qJD(4) * t14 - t50 * t24 - t48 * t26;
t4 = qJD(4) * t13 - t48 * t24 + t50 * t26;
t1 = [m(2) * qJDD(1) + (-t14 * t43 - t4 * t44) * mrSges(5,2) + (t13 * t43 + t5 * t44) * mrSges(5,1) + m(3) * (t31 * t51 + t32 * t49) + m(4) * (t11 * t28 + t12 * t29 - t16 * t24 + t17 * t26) + m(5) * (t3 * t13 + t2 * t14 + t7 * t4 + t6 * t5) + (t51 * mrSges(3,1) + t28 * mrSges(4,1) - t29 * mrSges(4,2) - t60) * qJDD(2) + (-m(2) - m(3) - t64) * g(1) + (-t24 * mrSges(4,1) - t26 * mrSges(4,2) - qJD(2) * t54) * qJD(2); t31 * mrSges(3,1) + t11 * mrSges(4,1) - t32 * mrSges(3,2) - t12 * mrSges(4,2) + (t40 * mrSges(4,2) + t55 * t39 + t65 * t49 + t59) * g(2) + (t39 * mrSges(4,2) - t55 * t40 - t65 * t51 - t34 + t60) * g(1) + (t21 * t43 + t66 * t44 + t62) * mrSges(5,1) + (Ifges(3,3) + Ifges(4,3) + (mrSges(4,1) * t47 - mrSges(4,2) * t46) * pkin(2)) * qJDD(2) + (t23 * mrSges(4,1) + t25 * mrSges(4,2) + qJD(1) * t54) * qJD(2) + (-t22 * t43 - t67 * t44 + t53) * mrSges(5,2) + t61 + (t2 * t22 + t3 * t21 + t66 * t6 + t67 * t7) * m(5) + (t16 * t23 - t17 * t25 + (t11 * t47 + t12 * t46) * pkin(2)) * m(4); t64 * (g(3) + qJDD(3)); -g(1) * t34 + (t7 * t44 + t62) * mrSges(5,1) + (t6 * t44 + t53) * mrSges(5,2) + t61;];
tau  = t1;
