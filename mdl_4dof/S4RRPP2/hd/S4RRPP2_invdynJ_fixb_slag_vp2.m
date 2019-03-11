% Calculate vector of inverse dynamics joint torques for
% S4RRPP2
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
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:59
% EndTime: 2019-03-08 18:34:00
% DurationCPUTime: 0.54s
% Computational Cost: add. (356->98), mult. (473->107), div. (0->0), fcn. (150->6), ass. (0->57)
t90 = -mrSges(4,1) - mrSges(5,1);
t95 = mrSges(3,1) - t90;
t91 = -mrSges(4,3) - mrSges(5,2);
t45 = cos(qJ(2));
t70 = pkin(1) * qJD(1);
t61 = t45 * t70;
t52 = qJD(3) - t61;
t93 = mrSges(3,2) + t91;
t41 = qJD(1) + qJD(2);
t92 = t95 * t41;
t40 = qJDD(1) + qJDD(2);
t89 = t40 * qJ(3) + t41 * qJD(3);
t43 = sin(qJ(2));
t62 = t43 * t70;
t16 = t41 * qJ(3) + t62;
t69 = pkin(1) * qJD(2);
t60 = qJD(1) * t69;
t67 = pkin(1) * qJDD(1);
t14 = t43 * t67 + t45 * t60;
t6 = t14 + t89;
t87 = t6 * qJ(3) + t52 * t16;
t42 = qJ(1) + qJ(2);
t37 = sin(t42);
t38 = cos(t42);
t47 = -pkin(2) - pkin(3);
t85 = t93 * t38 + (-m(5) * t47 + t95) * t37;
t84 = t93 * t37 - t38 * t95;
t63 = t45 * t69;
t17 = qJD(3) + t63;
t81 = t43 * pkin(1);
t30 = qJ(3) + t81;
t82 = t16 * t17 + t6 * t30;
t44 = sin(qJ(1));
t80 = t44 * pkin(1);
t79 = t45 * pkin(1);
t46 = cos(qJ(1));
t39 = t46 * pkin(1);
t75 = t40 * mrSges(4,1);
t74 = t40 * mrSges(5,1);
t71 = t38 * pkin(2) + t37 * qJ(3);
t65 = t39 + t71;
t64 = t43 * t69;
t33 = -pkin(2) - t79;
t22 = t38 * qJ(3);
t59 = -t37 * pkin(2) + t22;
t56 = t41 * t61;
t13 = -t43 * t60 + t45 * t67;
t53 = -g(1) * t37 + g(2) * t38;
t50 = qJDD(3) - t13;
t5 = t47 * t40 + t50;
t7 = -t40 * pkin(2) + t50;
t48 = t13 * mrSges(3,1) - t7 * mrSges(4,1) - t5 * mrSges(5,1) - t14 * mrSges(3,2) - t91 * t6 + (Ifges(4,2) + Ifges(3,3) + Ifges(5,3)) * t40;
t31 = t38 * pkin(3);
t20 = -pkin(3) + t33;
t15 = -t41 * pkin(2) + t52;
t9 = t47 * t41 + t52;
t1 = [m(3) * (t13 * t45 + t14 * t43) * pkin(1) + t48 + m(4) * (t15 * t64 + t7 * t33 + t82) + m(5) * (t5 * t20 + t9 * t64 + t82) - t20 * t74 - t33 * t75 + t40 * mrSges(3,1) * t79 + Ifges(2,3) * qJDD(1) + (-t40 * t81 - t41 * t63) * mrSges(3,2) - t64 * t92 + (-m(3) * t39 - m(4) * t65 - m(5) * (t31 + t65) - t46 * mrSges(2,1) + t44 * mrSges(2,2) + t84) * g(2) + (m(3) * t80 - m(4) * (t59 - t80) - m(5) * (t22 - t80) + t44 * mrSges(2,1) + t46 * mrSges(2,2) + t85) * g(1) - t91 * (t17 * t41 + t30 * t40); mrSges(3,2) * t56 + pkin(2) * t75 - t47 * t74 + t48 + (t5 * t47 - t9 * t62 + t87) * m(5) + (-t7 * pkin(2) - t15 * t62 + t87) * m(4) + t62 * t92 + (-m(4) * t71 - m(5) * (t31 + t71) + t84) * g(2) + (-m(4) * t59 - m(5) * t22 + t85) * g(1) - t91 * (-t56 + t89); t90 * t40 + (t5 + t53) * m(5) + (t53 + t7) * m(4) + (t91 * t41 + (-m(4) - m(5)) * t16) * t41; (g(3) + qJDD(4)) * m(5);];
tau  = t1;
