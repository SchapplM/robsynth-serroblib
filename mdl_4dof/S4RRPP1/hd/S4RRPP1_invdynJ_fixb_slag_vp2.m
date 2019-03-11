% Calculate vector of inverse dynamics joint torques for
% S4RRPP1
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:08
% EndTime: 2019-03-08 18:33:09
% DurationCPUTime: 0.39s
% Computational Cost: add. (425->114), mult. (698->139), div. (0->0), fcn. (318->10), ass. (0->60)
t74 = mrSges(4,1) + mrSges(5,1);
t81 = mrSges(4,2) - mrSges(5,3);
t78 = m(4) + m(5);
t57 = cos(pkin(6));
t60 = cos(qJ(2));
t73 = pkin(1) * qJD(1);
t56 = sin(pkin(6));
t58 = sin(qJ(2));
t76 = t56 * t58;
t18 = (t57 * t60 - t76) * t73;
t80 = qJD(4) - t18;
t79 = m(3) * pkin(1);
t55 = qJ(1) + qJ(2);
t51 = cos(t55);
t43 = pkin(2) * t51;
t77 = t60 * pkin(1);
t75 = t57 * t58;
t72 = pkin(1) * qJD(2);
t67 = t58 * t72;
t25 = -qJD(1) * t67 + qJDD(1) * t77;
t53 = qJDD(1) + qJDD(2);
t13 = t53 * pkin(2) + t25;
t71 = qJD(2) * t60;
t26 = (qJD(1) * t71 + qJDD(1) * t58) * pkin(1);
t6 = t56 * t13 + t57 * t26;
t44 = pkin(2) + t77;
t21 = pkin(1) * t75 + t56 * t44;
t49 = pkin(6) + t55;
t40 = sin(t49);
t41 = cos(t49);
t69 = t41 * pkin(3) + t40 * qJ(4) + t43;
t68 = t58 * t73;
t66 = mrSges(3,1) * t58 + mrSges(3,2) * t60;
t5 = t57 * t13 - t56 * t26;
t54 = qJD(1) + qJD(2);
t27 = t54 * pkin(2) + t60 * t73;
t10 = t56 * t27 + t57 * t68;
t19 = t57 * pkin(1) * t71 - t56 * t67;
t20 = -pkin(1) * t76 + t57 * t44;
t65 = pkin(1) * (t56 * t60 + t75);
t9 = t57 * t27 - t56 * t68;
t50 = sin(t55);
t64 = -t51 * mrSges(3,1) + t50 * mrSges(3,2) + t81 * t40 - t74 * t41;
t2 = t53 * qJ(4) + t54 * qJD(4) + t6;
t3 = -t53 * pkin(3) + qJDD(4) - t5;
t63 = t25 * mrSges(3,1) + t5 * mrSges(4,1) - t3 * mrSges(5,1) - t26 * mrSges(3,2) - t6 * mrSges(4,2) + t2 * mrSges(5,3) + (Ifges(5,2) + Ifges(3,3) + Ifges(4,3)) * t53;
t62 = t51 * mrSges(3,2) + (t78 * pkin(2) + mrSges(3,1)) * t50 + (m(5) * pkin(3) + t74) * t40 + (-m(5) * qJ(4) + t81) * t41;
t61 = cos(qJ(1));
t59 = sin(qJ(1));
t52 = t61 * pkin(1);
t42 = -t57 * pkin(2) - pkin(3);
t38 = t56 * pkin(2) + qJ(4);
t17 = qJD(2) * t65;
t16 = qJD(1) * t65;
t15 = -pkin(3) - t20;
t14 = qJ(4) + t21;
t12 = qJD(4) + t19;
t8 = t54 * qJ(4) + t10;
t7 = -t54 * pkin(3) + qJD(4) - t9;
t1 = [(t20 * mrSges(4,1) - t15 * mrSges(5,1) - t21 * mrSges(4,2) + t14 * mrSges(5,3) + (mrSges(3,1) * t60 - mrSges(3,2) * t58) * pkin(1)) * t53 + (-t19 * mrSges(4,2) + t12 * mrSges(5,3) - t74 * t17 - t66 * t72) * t54 + (-m(5) * (t52 + t69) - m(4) * (t43 + t52) + t59 * mrSges(2,2) + (-mrSges(2,1) - t79) * t61 + t64) * g(2) + (t25 * t60 + t26 * t58) * t79 + m(5) * (t8 * t12 + t2 * t14 + t3 * t15 + t7 * t17) + m(4) * (t10 * t19 - t9 * t17 + t5 * t20 + t6 * t21) + (t61 * mrSges(2,2) + (mrSges(2,1) + (m(3) + t78) * pkin(1)) * t59 + t62) * g(1) + t63 + Ifges(2,3) * qJDD(1); (t18 * mrSges(4,2) + t80 * mrSges(5,3) + t74 * t16 + t66 * t73) * t54 + (-t42 * mrSges(5,1) + t38 * mrSges(5,3) + (mrSges(4,1) * t57 - mrSges(4,2) * t56) * pkin(2)) * t53 + t64 * g(2) + t62 * g(1) + t63 + (-t69 * g(2) - t7 * t16 + t2 * t38 + t3 * t42 + t80 * t8) * m(5) + ((t5 * t57 + t56 * t6) * pkin(2) - t10 * t18 + t9 * t16 - t43 * g(2)) * m(4); t78 * (-g(3) + qJDD(3)); -t54 ^ 2 * mrSges(5,3) - t53 * mrSges(5,1) + (-g(1) * t40 + g(2) * t41 - t8 * t54 + t3) * m(5);];
tau  = t1;
