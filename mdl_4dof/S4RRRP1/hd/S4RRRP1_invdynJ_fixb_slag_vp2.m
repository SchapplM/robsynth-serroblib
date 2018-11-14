% Calculate vector of inverse dynamics joint torques for
% S4RRRP1
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4RRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:25
% EndTime: 2018-11-14 13:54:26
% DurationCPUTime: 0.50s
% Computational Cost: add. (640->116), mult. (1129->150), div. (0->0), fcn. (526->10), ass. (0->62)
t75 = mrSges(4,1) + mrSges(5,1);
t50 = qJD(1) + qJD(2);
t56 = cos(qJ(2));
t72 = pkin(1) * qJD(1);
t67 = t56 * t72;
t28 = pkin(2) * t50 + t67;
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t53 = sin(qJ(2));
t68 = t53 * t72;
t15 = t55 * t28 - t52 * t68;
t16 = t28 * t52 + t55 * t68;
t77 = t52 * t53;
t63 = t55 * t56 - t77;
t21 = t63 * t72;
t80 = pkin(1) * t56;
t26 = -qJD(2) * t68 + qJDD(1) * t80;
t49 = qJDD(1) + qJDD(2);
t18 = pkin(2) * t49 + t26;
t71 = qJD(2) * t56;
t27 = (qJD(1) * t71 + qJDD(1) * t53) * pkin(1);
t6 = t15 * qJD(3) + t18 * t52 + t27 * t55;
t69 = qJD(3) * t55;
t83 = -t16 * t21 + (t16 * t69 + t52 * t6) * pkin(2);
t82 = m(4) + m(5);
t40 = pkin(2) + t80;
t70 = qJD(3) * t52;
t11 = t40 * t69 + (t63 * qJD(2) - t53 * t70) * pkin(1);
t76 = t53 * t55;
t23 = pkin(1) * t76 + t40 * t52;
t81 = t16 * t11 + t6 * t23;
t51 = qJ(1) + qJ(2);
t46 = cos(t51);
t36 = pkin(2) * t46;
t54 = sin(qJ(1));
t79 = g(1) * t54;
t74 = mrSges(4,2) + mrSges(5,2);
t57 = cos(qJ(1));
t73 = t57 * pkin(1) + t36;
t47 = qJ(3) + t51;
t38 = cos(t47);
t66 = t74 * t38;
t22 = -pkin(1) * t77 + t55 * t40;
t64 = -t52 * t56 - t76;
t37 = sin(t47);
t62 = t74 * t37 - t75 * t38;
t43 = qJDD(3) + t49;
t7 = -t16 * qJD(3) + t55 * t18 - t27 * t52;
t3 = pkin(3) * t43 + t7;
t61 = t7 * mrSges(4,1) + t3 * mrSges(5,1) - t74 * t6 + (Ifges(4,3) + Ifges(5,3)) * t43;
t60 = t26 * mrSges(3,1) + Ifges(3,3) * t49 + t61;
t45 = sin(t51);
t59 = -t46 * mrSges(3,1) + t45 * mrSges(3,2) + t62;
t58 = t46 * mrSges(3,2) + t66 + (t82 * pkin(2) + mrSges(3,1)) * t45 + (m(5) * pkin(3) + t75) * t37;
t44 = qJD(3) + t50;
t39 = pkin(2) * t55 + pkin(3);
t32 = pkin(3) * t38;
t20 = t64 * t72;
t19 = pkin(3) + t22;
t13 = pkin(3) * t44 + t15;
t12 = -t40 * t70 + (t64 * qJD(2) - t53 * t69) * pkin(1);
t1 = [-t27 * mrSges(3,2) + Ifges(2,3) * qJDD(1) + m(5) * (t12 * t13 + t19 * t3 + t81) + m(4) * (t12 * t15 + t22 * t7 + t81) + (-t74 * t11 + t75 * t12) * t44 + (t22 * mrSges(4,1) + t19 * mrSges(5,1) - t74 * t23) * t43 + (-m(4) * t73 - m(5) * (t32 + t73) - t57 * mrSges(2,1) + t54 * mrSges(2,2) + t59) * g(2) + (t54 * mrSges(2,1) + t57 * mrSges(2,2) + t58) * g(1) + (t82 * t79 + (-t49 * t53 - t50 * t71) * mrSges(3,2) + (-qJD(2) * t50 * t53 + t49 * t56) * mrSges(3,1) + (-g(2) * t57 + t26 * t56 + t27 * t53 + t79) * m(3)) * pkin(1) + t60; t59 * g(2) + t58 * g(1) + (t74 * t21 - t75 * t20 + (-t75 * t52 - t74 * t55) * qJD(3) * pkin(2)) * t44 + t60 + t50 * mrSges(3,1) * t68 + (t50 * t67 - t27) * mrSges(3,2) + (mrSges(5,1) * t39 + (mrSges(4,1) * t55 - t74 * t52) * pkin(2)) * t43 + (t3 * t39 + (-t32 - t36) * g(2) + (-pkin(2) * t70 - t20) * t13 + t83) * m(5) + (-t15 * t20 - t36 * g(2) + (-t15 * t70 + t55 * t7) * pkin(2) + t83) * m(4); t62 * g(2) + (t75 * t37 + t66) * g(1) + (t43 * mrSges(5,1) + (g(1) * t37 - g(2) * t38 + t3) * m(5)) * pkin(3) + t74 * t15 * t44 + t61 + (-m(5) * (-t13 + t15) + t75 * t44) * t16; (-g(3) + qJDD(4)) * m(5);];
tau  = t1;
