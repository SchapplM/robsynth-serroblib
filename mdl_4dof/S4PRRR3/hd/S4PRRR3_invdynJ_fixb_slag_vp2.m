% Calculate vector of inverse dynamics joint torques for
% S4PRRR3
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:34
% EndTime: 2019-12-31 16:31:36
% DurationCPUTime: 0.56s
% Computational Cost: add. (609->123), mult. (958->178), div. (0->0), fcn. (433->8), ass. (0->63)
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t37 = -t60 * mrSges(5,1) + t58 * mrSges(5,2);
t103 = mrSges(4,1) - t37;
t97 = t58 / 0.2e1;
t57 = qJD(2) + qJD(3);
t59 = sin(qJ(3));
t82 = qJD(2) * t59;
t33 = pkin(2) * t82 + t57 * pkin(6);
t17 = t60 * qJD(1) - t58 * t33;
t61 = cos(qJ(3));
t83 = qJD(2) * pkin(2);
t78 = t61 * t83;
t93 = t59 * pkin(2);
t29 = qJD(3) * t78 + qJDD(2) * t93;
t55 = qJDD(2) + qJDD(3);
t22 = t55 * pkin(6) + t29;
t2 = t17 * qJD(4) + t58 * qJDD(1) + t60 * t22;
t18 = t58 * qJD(1) + t60 * t33;
t3 = -t18 * qJD(4) + t60 * qJDD(1) - t58 * t22;
t102 = t2 * t60 - t3 * t58;
t101 = -m(5) * pkin(6) + mrSges(4,2) - mrSges(5,3);
t100 = m(5) * pkin(3) + t103;
t34 = -t57 * pkin(3) - t78;
t73 = mrSges(5,1) * t58 + mrSges(5,2) * t60;
t99 = t34 * t73 + qJD(4) * (Ifges(5,5) * t60 - Ifges(5,6) * t58) / 0.2e1 - (t17 * t60 + t18 * t58) * mrSges(5,3);
t28 = (qJD(3) * t82 - qJDD(2) * t61) * pkin(2);
t89 = t57 * t58;
t31 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t89;
t88 = t57 * t60;
t32 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t88;
t84 = t60 * t61;
t86 = t58 * t61;
t98 = t31 * t86 - m(5) * (-t17 * t86 + t18 * t84 + t34 * t59) - t32 * t84 + (t61 * mrSges(4,2) + t103 * t59) * t57;
t92 = Ifges(5,4) * t58;
t91 = Ifges(5,4) * t60;
t90 = Ifges(5,2) * t60;
t81 = qJD(4) * t58;
t80 = qJD(4) * t60;
t79 = m(2) + m(3) + m(4);
t56 = pkin(7) + qJ(2);
t21 = -t55 * pkin(3) + t28;
t24 = t60 * t55 - t57 * t81;
t25 = t58 * t55 + t57 * t80;
t77 = m(5) * t21 - t24 * mrSges(5,1) + t25 * mrSges(5,2);
t52 = sin(t56);
t53 = cos(t56);
t75 = g(1) * t52 - g(2) * t53;
t72 = t90 + t92;
t68 = t58 * (Ifges(5,1) * t60 - t92);
t54 = qJ(3) + t56;
t46 = sin(t54);
t47 = cos(t54);
t66 = -t100 * t47 + t101 * t46;
t65 = t100 * t46 + t101 * t47;
t13 = -qJDD(4) * mrSges(5,2) + t24 * mrSges(5,3);
t14 = qJDD(4) * mrSges(5,1) - t25 * mrSges(5,3);
t63 = -t31 * t80 - t32 * t81 - t58 * t14 + t60 * t13 + m(5) * (-t17 * t80 - t18 * t81 + t102);
t19 = Ifges(5,6) * qJD(4) + t72 * t57;
t38 = Ifges(5,4) * t88;
t20 = Ifges(5,1) * t89 + Ifges(5,5) * qJD(4) + t38;
t62 = -t19 * t81 / 0.2e1 - t28 * mrSges(4,1) + (Ifges(5,1) * t25 + Ifges(5,4) * t24) * t97 + t60 * (Ifges(5,4) * t25 + Ifges(5,2) * t24) / 0.2e1 + Ifges(4,3) * t55 + t21 * t37 + t24 * t72 / 0.2e1 + t25 * (Ifges(5,1) * t58 + t91) / 0.2e1 - t29 * mrSges(4,2) + (t20 + t57 * (-Ifges(5,2) * t58 + t91)) * t80 / 0.2e1 + t102 * mrSges(5,3) + (0.2e1 * Ifges(5,5) * t97 + Ifges(5,6) * t60) * qJDD(4) + (t68 * t57 / 0.2e1 + t99) * qJD(4);
t1 = [t32 * t80 - t31 * t81 + t58 * t13 + t60 * t14 + m(5) * (t2 * t58 + t3 * t60 + (-t17 * t58 + t18 * t60) * qJD(4)) + t79 * qJDD(1) + (-m(5) - t79) * g(3); t62 + t77 * (-t61 * pkin(2) - pkin(3)) + t63 * (pkin(6) + t93) + ((t61 * mrSges(4,1) - t59 * mrSges(4,2)) * t55 + t75 * m(5) + (-t28 * t61 + t29 * t59 + t75) * m(4) - t98 * qJD(3)) * pkin(2) + (t52 * mrSges(3,1) + t53 * mrSges(3,2) + t65) * g(1) + (-t53 * mrSges(3,1) + t52 * mrSges(3,2) + t66) * g(2) + Ifges(3,3) * qJDD(2); -t77 * pkin(3) + t63 * pkin(6) + t65 * g(1) + t66 * g(2) + t98 * t83 + t62; t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t25 + Ifges(5,6) * t24 + Ifges(5,3) * qJDD(4) + g(3) * t37 - t17 * t32 + t18 * t31 + (t19 * t97 + (-t68 / 0.2e1 + t90 * t97) * t57 - (t20 + t38) * t60 / 0.2e1 - t99) * t57 + (g(1) * t47 + g(2) * t46) * t73;];
tau = t1;
