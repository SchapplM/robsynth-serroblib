% Calculate vector of inverse dynamics joint torques for
% S4RPPR2
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4RPPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:26
% EndTime: 2018-11-14 13:47:27
% DurationCPUTime: 0.67s
% Computational Cost: add. (516->119), mult. (817->145), div. (0->0), fcn. (404->8), ass. (0->53)
t69 = qJD(1) - qJD(4);
t68 = m(4) + m(5);
t44 = sin(pkin(6));
t45 = cos(pkin(6));
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t23 = t44 * t48 + t45 * t46;
t67 = t69 * t23;
t21 = -t44 * t46 + t45 * t48;
t66 = t69 * t21;
t59 = m(3) + t68;
t50 = -pkin(1) - pkin(2);
t42 = -qJDD(1) + qJDD(4);
t65 = Ifges(5,3) * t42;
t64 = mrSges(2,1) + mrSges(3,1);
t62 = qJ(2) * t44;
t60 = qJ(2) * qJD(1);
t58 = qJD(1) * qJD(2);
t57 = pkin(6) + qJ(4);
t30 = t50 * qJDD(1) + qJDD(2);
t35 = qJ(2) * qJDD(1) + t58;
t10 = t45 * t30 - t35 * t44;
t27 = t45 * t50 - t62;
t56 = cos(t57);
t55 = sin(t57);
t34 = t50 * qJD(1) + qJD(2);
t29 = t45 * t34;
t12 = t29 + (-pkin(3) - t62) * qJD(1);
t14 = t34 * t44 + t45 * t60;
t5 = t12 * t48 - t14 * t46;
t6 = t12 * t46 + t14 * t48;
t54 = -(-t44 * t60 + t29) * t44 + t14 * t45;
t25 = -pkin(3) + t27;
t28 = qJ(2) * t45 + t44 * t50;
t7 = t25 * t48 - t28 * t46;
t8 = t25 * t46 + t28 * t48;
t53 = -m(5) * pkin(3) * t44 + mrSges(2,2) - mrSges(3,3);
t52 = mrSges(4,1) * t44 + mrSges(4,2) * t45 + mrSges(3,3);
t49 = cos(qJ(1));
t47 = sin(qJ(1));
t38 = -pkin(1) * qJDD(1) + qJDD(2);
t36 = pkin(3) * t45 + pkin(2);
t24 = t44 * t47 + t45 * t49;
t22 = t44 * t49 - t45 * t47;
t16 = -t47 * t56 + t49 * t55;
t15 = -t47 * t55 - t49 * t56;
t11 = t30 * t44 + t35 * t45;
t9 = -qJDD(1) * pkin(3) + t10;
t4 = -t23 * qJD(2) - t8 * qJD(4);
t3 = t21 * qJD(2) + t7 * qJD(4);
t2 = -t6 * qJD(4) - t46 * t11 + t48 * t9;
t1 = t5 * qJD(4) + t48 * t11 + t46 * t9;
t13 = [-t38 * mrSges(3,1) - t10 * mrSges(4,1) + t11 * mrSges(4,2) + t35 * mrSges(3,3) - t65 + t52 * t58 + (t3 * t69 - t42 * t8 + t1) * mrSges(5,2) + (-t4 * t69 + t42 * t7 - t2) * mrSges(5,1) + m(5) * (t1 * t8 + t2 * t7 + t3 * t6 + t4 * t5) + m(3) * (-t38 * pkin(1) + (t35 + t58) * qJ(2)) + m(4) * (t54 * qJD(2) + t10 * t27 + t11 * t28) + (t15 * mrSges(5,1) + t16 * mrSges(5,2) - t24 * mrSges(4,1) + t22 * mrSges(4,2) + (-m(4) * pkin(2) - m(5) * t36 - t64) * t49 + t53 * t47 - t59 * (t49 * pkin(1) + t47 * qJ(2))) * g(2) + (-t22 * mrSges(4,1) - t16 * mrSges(5,1) - t24 * mrSges(4,2) + t15 * mrSges(5,2) + (-m(4) * t50 - m(5) * (-pkin(1) - t36) + m(3) * pkin(1) + t64) * t47 + (-t59 * qJ(2) + t53) * t49) * g(1) + (mrSges(3,1) * pkin(1) - mrSges(4,1) * t27 + mrSges(4,2) * t28 + mrSges(3,3) * qJ(2) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1); (mrSges(5,1) * t21 - mrSges(5,2) * t23) * t42 + (-mrSges(4,1) * t45 + mrSges(4,2) * t44 - mrSges(3,1)) * qJDD(1) + (-m(3) * qJ(2) - t52) * qJD(1) ^ 2 - (t67 * mrSges(5,1) + t66 * mrSges(5,2)) * t69 + m(3) * t38 + (-g(1) * t47 + g(2) * t49) * t59 + (t1 * t23 + t2 * t21 + t67 * t5 - t66 * t6) * m(5) + (-t54 * qJD(1) + t10 * t45 + t11 * t44) * m(4); t68 * (g(3) + qJDD(3)); t65 + (-g(1) * t15 - g(2) * t16 - t5 * t69 - t1) * mrSges(5,2) + (g(1) * t16 - g(2) * t15 - t6 * t69 + t2) * mrSges(5,1);];
tau  = t13;
