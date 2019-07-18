% Calculate vector of inverse dynamics joint torques for
% S4RRPR2
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
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:31
% EndTime: 2019-07-18 18:16:32
% DurationCPUTime: 0.69s
% Computational Cost: add. (804->147), mult. (1023->182), div. (0->0), fcn. (434->8), ass. (0->64)
t84 = mrSges(3,1) + mrSges(4,1);
t91 = mrSges(3,2) - mrSges(4,3);
t61 = sin(qJ(2));
t90 = t84 * t61;
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t66 = -pkin(2) - pkin(3);
t30 = -t60 * qJ(3) + t63 * t66;
t64 = cos(qJ(2));
t82 = pkin(1) * qJD(1);
t89 = t63 * qJD(3) + t30 * qJD(4) - (t60 * t61 + t63 * t64) * t82;
t31 = t63 * qJ(3) + t60 * t66;
t88 = -t60 * qJD(3) - t31 * qJD(4) - (-t60 * t64 + t61 * t63) * t82;
t87 = -m(4) - m(5);
t57 = qJDD(1) + qJDD(2);
t86 = pkin(2) * t57;
t58 = qJD(1) + qJD(2);
t72 = -t64 * t82 + qJD(3);
t85 = (-pkin(2) * t58 + t72) * t61;
t81 = pkin(1) * qJD(2);
t76 = qJD(1) * t81;
t80 = pkin(1) * qJDD(1);
t27 = t61 * t80 + t64 * t76;
t59 = qJ(1) + qJ(2);
t54 = sin(t59);
t55 = cos(t59);
t83 = t55 * pkin(2) + t54 * qJ(3);
t53 = qJD(4) - t58;
t79 = qJD(4) * t53;
t65 = cos(qJ(1));
t78 = t65 * pkin(1) + t83;
t77 = t61 * t81;
t48 = -pkin(1) * t64 - pkin(2);
t71 = qJ(3) * t57 + qJD(3) * t58;
t10 = t71 + t27;
t20 = t66 * t58 + t72;
t29 = qJ(3) * t58 + t61 * t82;
t6 = t20 * t63 - t29 * t60;
t26 = -t61 * t76 + t64 * t80;
t70 = qJDD(3) - t26;
t9 = t66 * t57 + t70;
t2 = qJD(4) * t6 + t10 * t63 + t60 * t9;
t52 = qJDD(4) - t57;
t75 = -t2 * mrSges(5,2) + Ifges(5,3) * t52;
t74 = -g(1) * t54 + g(2) * t55;
t7 = t20 * t60 + t29 * t63;
t73 = -t6 * t60 + t63 * t7;
t38 = -pkin(3) + t48;
t45 = pkin(1) * t61 + qJ(3);
t12 = t38 * t63 - t45 * t60;
t13 = t38 * t60 + t45 * t63;
t23 = -t54 * t60 - t55 * t63;
t24 = -t54 * t63 + t55 * t60;
t69 = t23 * mrSges(5,1) + t24 * mrSges(5,2) + t91 * t54 - t84 * t55;
t68 = t26 * mrSges(3,1) - t27 * mrSges(3,2) + t10 * mrSges(4,3) - t75 + (Ifges(4,2) + Ifges(3,3)) * t57;
t67 = -t24 * mrSges(5,1) + t23 * mrSges(5,2) + (m(4) * pkin(2) - m(5) * t66 + t84) * t54 + (qJ(3) * t87 + t91) * t55;
t62 = sin(qJ(1));
t46 = t55 * pkin(3);
t35 = t64 * t81 + qJD(3);
t11 = t70 - t86;
t5 = -t13 * qJD(4) - t35 * t60 + t63 * t77;
t4 = t12 * qJD(4) + t35 * t63 + t60 * t77;
t3 = -qJD(4) * t7 - t10 * t60 + t63 * t9;
t1 = [m(5) * (t12 * t3 + t13 * t2 + t4 * t7 + t5 * t6) + (t35 * t58 + t45 * t57) * mrSges(4,3) + (-t13 * t52 - t4 * t53) * mrSges(5,2) + (t12 * t52 + t5 * t53 - t3) * mrSges(5,1) + ((t64 * mrSges(3,1) - t61 * mrSges(3,2)) * t57 + (-g(2) * t65 + t26 * t64 + t27 * t61) * m(3) + (m(4) * t85 + (-t64 * mrSges(3,2) - t90) * t58) * qJD(2)) * pkin(1) + (-t65 * mrSges(2,1) + t62 * mrSges(2,2) - m(4) * t78 - m(5) * (t46 + t78) + t69) * g(2) + t68 + m(4) * (t10 * t45 + t11 * t48 + t29 * t35) + Ifges(2,3) * qJDD(1) + (-t48 * t57 - t11) * mrSges(4,1) + (t65 * mrSges(2,2) + (mrSges(2,1) + (m(3) - t87) * pkin(1)) * t62 + t67) * g(1); (t30 * t52 + t53 * t88 - t3) * mrSges(5,1) + (-t31 * t52 - t53 * t89) * mrSges(5,2) + (-t11 + t86) * mrSges(4,1) + (-m(4) * t83 + t69) * g(2) + t71 * mrSges(4,3) + m(4) * (-pkin(2) * t11 + qJ(3) * t10 + qJD(3) * t29) + t68 + (-m(4) * (t29 * t64 + t85) + (t91 * t64 + t90) * t58) * t82 + t67 * g(1) + ((-t46 - t83) * g(2) + t2 * t31 + t3 * t30 + t89 * t7 + t88 * t6) * m(5); -t57 * mrSges(4,1) + (-t52 * t60 - t63 * t79) * mrSges(5,2) + (t52 * t63 - t60 * t79) * mrSges(5,1) + (t73 * qJD(4) + t2 * t60 + t3 * t63 + t74) * m(5) + (t11 + t74) * m(4) + (-m(4) * t29 - m(5) * t73 - mrSges(4,3) * t58 + (mrSges(5,1) * t60 + mrSges(5,2) * t63) * t53) * t58; (-g(1) * t23 - g(2) * t24 + t53 * t6) * mrSges(5,2) + (g(1) * t24 - g(2) * t23 + t53 * t7 + t3) * mrSges(5,1) + t75;];
tau  = t1;
