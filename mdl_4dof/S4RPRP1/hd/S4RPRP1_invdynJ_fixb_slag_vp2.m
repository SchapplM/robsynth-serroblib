% Calculate vector of inverse dynamics joint torques for
% S4RPRP1
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4RPRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:27
% EndTime: 2018-11-14 13:48:28
% DurationCPUTime: 0.38s
% Computational Cost: add. (374->99), mult. (654->115), div. (0->0), fcn. (304->10), ass. (0->55)
t68 = mrSges(4,1) + mrSges(5,1);
t74 = mrSges(4,2) - mrSges(5,3);
t73 = m(3) * pkin(1);
t50 = cos(qJ(3));
t46 = sin(pkin(6));
t71 = pkin(1) * t46;
t32 = t50 * t71;
t47 = cos(pkin(6));
t36 = t47 * pkin(1) + pkin(2);
t48 = sin(qJ(3));
t17 = t48 * t36 + t32;
t22 = t36 * qJD(1);
t60 = qJD(1) * t71;
t57 = t48 * t60;
t9 = t50 * t22 - t57;
t64 = qJD(4) - t9;
t21 = t36 * qJDD(1);
t63 = qJD(3) * t50;
t69 = t48 * t22;
t6 = -(qJD(1) * t63 + qJDD(1) * t48) * t71 - qJD(3) * t69 + t50 * t21;
t72 = m(4) + m(5);
t43 = qJDD(1) + qJDD(3);
t70 = t43 * pkin(3);
t45 = qJ(1) + pkin(6);
t41 = qJ(3) + t45;
t34 = sin(t41);
t35 = cos(t41);
t67 = t35 * pkin(3) + t34 * qJ(4);
t40 = cos(t45);
t51 = cos(qJ(1));
t66 = t51 * pkin(1) + pkin(2) * t40;
t65 = t43 * qJ(4);
t62 = m(3) + t72;
t61 = t48 * t71;
t3 = qJDD(4) - t6 - t70;
t59 = g(1) * t34 - t3;
t58 = t74 * t35;
t56 = t74 * t34 - t68 * t35;
t12 = -qJD(3) * t61 + t36 * t63;
t16 = t50 * t36 - t61;
t44 = qJD(1) + qJD(3);
t5 = -qJD(3) * t57 + qJDD(1) * t32 + t48 * t21 + t22 * t63;
t2 = t44 * qJD(4) + t5 + t65;
t55 = t6 * mrSges(4,1) - t5 * mrSges(4,2) + t2 * mrSges(5,3) + (Ifges(5,2) + Ifges(4,3)) * t43;
t10 = t50 * t60 + t69;
t49 = sin(qJ(1));
t39 = sin(t45);
t26 = t35 * qJ(4);
t15 = -pkin(3) - t16;
t14 = qJ(4) + t17;
t13 = t17 * qJD(3);
t11 = qJD(4) + t12;
t8 = t44 * qJ(4) + t10;
t7 = -t44 * pkin(3) + t64;
t1 = [-t3 * mrSges(5,1) + m(5) * (t8 * t11 + t7 * t13 + t2 * t14 + t3 * t15) + m(4) * (t10 * t12 - t9 * t13 + t6 * t16 + t5 * t17) + (-t12 * mrSges(4,2) + t11 * mrSges(5,3) - t13 * t68) * t44 + (t16 * mrSges(4,1) - t15 * mrSges(5,1) - t17 * mrSges(4,2) + t14 * mrSges(5,3)) * t43 + (-m(5) * (t66 + t67) - m(4) * t66 - t40 * mrSges(3,1) + t39 * mrSges(3,2) + t49 * mrSges(2,2) + (-mrSges(2,1) - t73) * t51 + t56) * g(2) + (-m(5) * t26 + t51 * mrSges(2,2) + t40 * mrSges(3,2) + (pkin(2) * t72 + mrSges(3,1)) * t39 + (pkin(1) * t62 + mrSges(2,1)) * t49 + (m(5) * pkin(3) + t68) * t34 + t58) * g(1) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t47 * mrSges(3,1) - 0.2e1 * t46 * mrSges(3,2) + (t46 ^ 2 + t47 ^ 2) * t73) * pkin(1)) * qJDD(1) + t55; (-g(3) + qJDD(2)) * t62; mrSges(5,3) * t65 + t56 * g(2) + (t34 * mrSges(4,1) + t58) * g(1) + (t59 + t70) * mrSges(5,1) + (t9 * mrSges(4,2) + mrSges(5,3) * t64 + t10 * t68) * t44 + t55 + (-t67 * g(2) + (t34 * pkin(3) - t26) * g(1) - t7 * t10 - t3 * pkin(3) + t2 * qJ(4) + t64 * t8) * m(5); -t44 ^ 2 * mrSges(5,3) - t43 * mrSges(5,1) + (g(2) * t35 - t8 * t44 - t59) * m(5);];
tau  = t1;
