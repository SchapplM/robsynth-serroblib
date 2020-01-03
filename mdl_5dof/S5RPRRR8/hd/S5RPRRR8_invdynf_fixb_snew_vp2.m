% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR8_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR8_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:45
% EndTime: 2019-12-31 19:05:46
% DurationCPUTime: 0.60s
% Computational Cost: add. (7142->116), mult. (9250->148), div. (0->0), fcn. (4361->8), ass. (0->64)
t45 = -qJDD(1) + qJDD(3);
t52 = sin(qJ(4));
t56 = cos(qJ(4));
t47 = -qJD(1) + qJD(3);
t75 = qJD(4) * t47;
t73 = t56 * t75;
t34 = t52 * t45 + t73;
t35 = t56 * t45 - t52 * t75;
t80 = t47 * t52;
t36 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t80;
t79 = t47 * t56;
t37 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t79;
t43 = t47 ^ 2;
t51 = sin(qJ(5));
t55 = cos(qJ(5));
t31 = (t51 * t56 + t52 * t55) * t47;
t16 = -t31 * qJD(5) - t51 * t34 + t55 * t35;
t30 = (-t51 * t52 + t55 * t56) * t47;
t17 = t30 * qJD(5) + t55 * t34 + t51 * t35;
t46 = qJD(4) + qJD(5);
t22 = -t46 * mrSges(6,2) + t30 * mrSges(6,3);
t23 = t46 * mrSges(6,1) - t31 * mrSges(6,3);
t38 = qJD(4) * pkin(4) - pkin(8) * t80;
t50 = t56 ^ 2;
t59 = qJD(1) ^ 2;
t54 = sin(qJ(1));
t58 = cos(qJ(1));
t69 = -t58 * g(1) - t54 * g(2);
t64 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t69;
t81 = -pkin(1) - pkin(2);
t27 = t59 * t81 + t64;
t72 = t54 * g(1) - t58 * g(2);
t62 = -t59 * qJ(2) + qJDD(2) - t72;
t29 = qJDD(1) * t81 + t62;
t53 = sin(qJ(3));
t57 = cos(qJ(3));
t70 = -t53 * t27 + t57 * t29;
t65 = -t45 * pkin(3) - t70;
t61 = t16 * mrSges(6,1) + t30 * t22 - m(6) * (t38 * t80 - t35 * pkin(4) + (-pkin(8) * t50 - pkin(7)) * t43 + t65) - t17 * mrSges(6,2) - t31 * t23;
t83 = (t52 * t36 - t56 * t37) * t47 + m(5) * (-t43 * pkin(7) + t65) - t35 * mrSges(5,1) + t34 * mrSges(5,2) - t61;
t82 = -m(3) - m(4);
t78 = mrSges(2,1) + mrSges(3,1);
t76 = t57 * t27 + t53 * t29;
t19 = -t43 * pkin(3) + t45 * pkin(7) + t76;
t77 = t52 * g(3) + t56 * t19;
t74 = -m(2) + t82;
t71 = t56 * g(3) - t52 * t19;
t33 = (-mrSges(5,1) * t56 + mrSges(5,2) * t52) * t47;
t10 = (-t34 + t73) * pkin(8) + (t43 * t52 * t56 + qJDD(4)) * pkin(4) + t71;
t11 = -t50 * t43 * pkin(4) + t35 * pkin(8) - qJD(4) * t38 + t77;
t21 = -t30 * mrSges(6,1) + t31 * mrSges(6,2);
t44 = qJDD(4) + qJDD(5);
t8 = m(6) * (t55 * t10 - t51 * t11) - t17 * mrSges(6,3) + t44 * mrSges(6,1) - t31 * t21 + t46 * t22;
t9 = m(6) * (t51 * t10 + t55 * t11) + t16 * mrSges(6,3) - t44 * mrSges(6,2) + t30 * t21 - t46 * t23;
t5 = m(5) * t71 + qJDD(4) * mrSges(5,1) - t34 * mrSges(5,3) + qJD(4) * t37 - t33 * t80 + t51 * t9 + t55 * t8;
t6 = m(5) * t77 - qJDD(4) * mrSges(5,2) + t35 * mrSges(5,3) - qJD(4) * t36 + t33 * t79 - t51 * t8 + t55 * t9;
t68 = -t56 * t5 - t52 * t6;
t4 = m(4) * t76 - t43 * mrSges(4,1) - t45 * mrSges(4,2) - t52 * t5 + t56 * t6;
t7 = m(4) * t70 + t45 * mrSges(4,1) - t43 * mrSges(4,2) - t83;
t66 = t57 * t4 - t53 * t7 + m(3) * (-t59 * pkin(1) + t64) + qJDD(1) * mrSges(3,3);
t63 = -m(3) * (-qJDD(1) * pkin(1) + t62) - t53 * t4 - t57 * t7;
t2 = m(2) * t72 + (-mrSges(2,2) + mrSges(3,3)) * t59 + t78 * qJDD(1) + t63;
t1 = m(2) * t69 - qJDD(1) * mrSges(2,2) - t59 * t78 + t66;
t3 = [-m(1) * g(1) + t58 * t1 - t54 * t2, t1, -t59 * mrSges(3,1) + t66, t4, t6, t9; -m(1) * g(2) + t54 * t1 + t58 * t2, t2, g(3) * t82 + t68, t7, t5, t8; (-m(1) + t74) * g(3) + t68, g(3) * t74 + t68, -qJDD(1) * mrSges(3,1) - t59 * mrSges(3,3) - t63, m(4) * g(3) - t68, t83, -t61;];
f_new = t3;
