% Calculate vector of inverse dynamics joint torques for
% S4RPRR1
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4RPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:29
% EndTime: 2018-11-14 13:50:30
% DurationCPUTime: 0.48s
% Computational Cost: add. (655->119), mult. (1347->152), div. (0->0), fcn. (726->14), ass. (0->65)
t76 = m(3) * pkin(1);
t50 = cos(pkin(7));
t37 = t50 * pkin(1) + pkin(2);
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t49 = sin(pkin(7));
t74 = pkin(1) * t49;
t22 = t55 * t37 - t52 * t74;
t75 = m(4) + m(5);
t48 = qJ(1) + pkin(7);
t44 = qJ(3) + t48;
t38 = qJ(4) + t44;
t31 = sin(t38);
t73 = g(1) * t31;
t35 = sin(t44);
t72 = g(1) * t35;
t26 = t37 * qJDD(1);
t62 = qJD(1) * qJD(3);
t27 = t37 * qJD(1);
t66 = qJD(3) * t27;
t12 = t55 * t66 + t52 * t26 + (qJDD(1) * t55 - t52 * t62) * t74;
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t60 = qJD(1) * t74;
t17 = t55 * t27 - t52 * t60;
t47 = qJD(1) + qJD(3);
t16 = t47 * pkin(3) + t17;
t18 = t52 * t27 + t55 * t60;
t69 = t54 * t18;
t7 = t51 * t16 + t69;
t13 = -t52 * t66 + t55 * t26 + (-qJDD(1) * t52 - t55 * t62) * t74;
t46 = qJDD(1) + qJDD(3);
t8 = t46 * pkin(3) + t13;
t3 = -t7 * qJD(4) - t51 * t12 + t54 * t8;
t40 = qJDD(4) + t46;
t71 = t3 * mrSges(5,1) + Ifges(5,3) * t40;
t70 = t51 * t18;
t32 = cos(t38);
t29 = t32 * mrSges(5,1);
t36 = cos(t44);
t68 = -t36 * mrSges(4,1) - t29;
t42 = cos(t48);
t56 = cos(qJ(1));
t67 = t56 * pkin(1) + pkin(2) * t42;
t65 = qJD(4) * t51;
t64 = qJD(4) * t54;
t63 = m(3) + t75;
t59 = t13 * mrSges(4,1) + Ifges(4,3) * t46 + t71;
t6 = t54 * t16 - t70;
t21 = pkin(3) + t22;
t23 = t52 * t37 + t55 * t74;
t14 = t54 * t21 - t51 * t23;
t15 = t51 * t21 + t54 * t23;
t2 = t6 * qJD(4) + t54 * t12 + t51 * t8;
t58 = g(1) * t32 + g(2) * t31 - t2;
t53 = sin(qJ(1));
t43 = qJD(4) + t47;
t41 = sin(t48);
t20 = t23 * qJD(3);
t19 = t22 * qJD(3);
t10 = t54 * t17 - t70;
t9 = -t51 * t17 - t69;
t5 = -t15 * qJD(4) - t51 * t19 - t54 * t20;
t4 = t14 * qJD(4) + t54 * t19 - t51 * t20;
t1 = [(t14 * t40 + t5 * t43) * mrSges(5,1) + (-t20 * t47 + t22 * t46) * mrSges(4,1) + m(5) * (t3 * t14 + t2 * t15 + t7 * t4 + t6 * t5) + m(4) * (t12 * t23 + t13 * t22 - t17 * t20 + t18 * t19) + (-t15 * t40 - t4 * t43 - t2) * mrSges(5,2) + (-t19 * t47 - t23 * t46 - t12) * mrSges(4,2) + (-t42 * mrSges(3,1) + t41 * mrSges(3,2) + t53 * mrSges(2,2) - m(4) * t67 + t35 * mrSges(4,2) - m(5) * (pkin(3) * t36 + t67) + t31 * mrSges(5,2) + (-mrSges(2,1) - t76) * t56 + t68) * g(2) + (t31 * mrSges(5,1) + t56 * mrSges(2,2) + t42 * mrSges(3,2) + t36 * mrSges(4,2) + t32 * mrSges(5,2) + (m(5) * pkin(3) + mrSges(4,1)) * t35 + (t75 * pkin(2) + mrSges(3,1)) * t41 + (t63 * pkin(1) + mrSges(2,1)) * t53) * g(1) + (Ifges(2,3) + Ifges(3,3) + (0.2e1 * t50 * mrSges(3,1) - 0.2e1 * t49 * mrSges(3,2) + (t49 ^ 2 + t50 ^ 2) * t76) * pkin(1)) * qJDD(1) + t59; (-g(3) + qJDD(2)) * t63; -m(5) * (t7 * t10 + t6 * t9) + t68 * g(2) + (-t9 * t43 + t73) * mrSges(5,1) + (t10 * t43 + t58) * mrSges(5,2) + (t18 * t47 + t72) * mrSges(4,1) + (g(1) * t36 + g(2) * t35 + t17 * t47 - t12) * mrSges(4,2) + ((-t51 * t40 - t43 * t64) * mrSges(5,2) + (t54 * t40 - t43 * t65) * mrSges(5,1) + (-g(2) * t36 + t2 * t51 + t3 * t54 - t6 * t65 + t7 * t64 + t72) * m(5)) * pkin(3) + t59; -g(2) * t29 + (t7 * t43 + t73) * mrSges(5,1) + (t6 * t43 + t58) * mrSges(5,2) + t71;];
tau  = t1;
