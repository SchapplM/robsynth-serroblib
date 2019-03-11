% Calculate vector of inverse dynamics joint torques for
% S4RPPP1
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
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPP1_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:19
% EndTime: 2019-03-08 18:26:21
% DurationCPUTime: 1.34s
% Computational Cost: add. (471->203), mult. (1384->254), div. (0->0), fcn. (924->6), ass. (0->95)
t121 = mrSges(4,1) + mrSges(3,3);
t120 = mrSges(5,2) + mrSges(4,3);
t118 = m(4) + m(5);
t119 = t118 * qJ(3) - mrSges(3,2) + t120;
t67 = cos(pkin(6));
t117 = pkin(1) * t67;
t65 = sin(pkin(6));
t66 = sin(pkin(4));
t116 = t65 * t66;
t115 = t66 * t67;
t68 = cos(pkin(4));
t114 = t68 * t65;
t70 = cos(qJ(1));
t113 = t68 * t70;
t69 = sin(qJ(1));
t112 = t69 * t67;
t111 = -mrSges(3,1) - mrSges(5,3);
t44 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t66;
t96 = qJDD(1) * t68;
t88 = pkin(1) * t96;
t14 = t67 * t44 + t65 * t88;
t109 = ((-mrSges(3,1) + mrSges(4,2)) * t68 + t121 * t116) * qJD(1);
t108 = (t120 * t68 + (mrSges(4,1) + mrSges(5,1)) * t115) * qJD(1);
t102 = qJD(1) * t66;
t86 = qJ(2) * t102;
t101 = qJD(1) * t68;
t89 = pkin(1) * t101;
t24 = t65 * t89 + t67 * t86;
t104 = qJ(2) * t66;
t35 = pkin(1) * t114 + t67 * t104;
t97 = qJDD(1) * t66;
t84 = t67 * t97;
t107 = mrSges(5,1) * t84 + mrSges(5,2) * t96;
t85 = t65 * t97;
t106 = mrSges(4,1) * t85 + mrSges(4,2) * t96;
t105 = t70 * pkin(1) + t69 * t104;
t103 = qJ(3) * t68;
t100 = qJD(2) * t66;
t99 = qJD(4) * t68;
t45 = t65 * t86;
t98 = qJD(3) + t45;
t27 = t65 * t44;
t95 = qJDD(3) + t27;
t94 = m(3) + t118;
t93 = qJD(1) * qJD(3);
t92 = 0.2e1 * t68;
t87 = -pkin(2) - t117;
t56 = -pkin(1) * t97 + qJDD(2);
t83 = -m(3) * t56 - mrSges(3,2) * t85;
t82 = -t69 * pkin(1) + t70 * t104;
t81 = -qJ(3) * t65 - pkin(1);
t79 = t87 * t68;
t77 = m(5) * qJ(4) - mrSges(4,2) - t111;
t76 = -qJD(3) * t65 - qJD(4) * t67;
t75 = -pkin(2) * t67 + t81;
t74 = pkin(3) * t115 + t103;
t20 = t75 * t66;
t73 = (-pkin(2) - qJ(4)) * t67 + t81;
t30 = -t67 * t113 + t65 * t69;
t32 = t68 * t112 + t65 * t70;
t72 = -g(1) * t32 - g(2) * t30 + g(3) * t115;
t12 = t73 * t66;
t71 = pkin(3) * t116 + (-qJ(4) + t87) * t68;
t57 = t65 * t104;
t49 = mrSges(4,2) * t84;
t47 = mrSges(5,1) * t85;
t43 = qJD(3) * t68 + t67 * t100;
t42 = t65 * t100 - t99;
t38 = (mrSges(5,1) * t116 - t68 * mrSges(5,3)) * qJD(1);
t37 = (-t68 * mrSges(3,2) + mrSges(3,3) * t115) * qJD(1);
t34 = t68 * t117 - t57;
t33 = -t69 * t114 + t67 * t70;
t31 = t65 * t113 + t112;
t29 = t76 * t66;
t23 = (-t65 * mrSges(5,2) - mrSges(5,3) * t67) * t102;
t22 = (mrSges(4,2) * t67 - t65 * mrSges(4,3)) * t102;
t21 = t67 * t89 - t45;
t19 = t57 + t79;
t18 = -t103 - t35;
t17 = -qJ(3) * t101 - t24;
t16 = qJD(1) * t20 + qJD(2);
t15 = qJD(1) * t79 + t98;
t13 = t67 * t88 - t27;
t11 = t74 + t35;
t10 = t57 + t71;
t9 = qJDD(1) * t79 + t95;
t8 = qJD(1) * t12 + qJD(2);
t7 = t74 * qJD(1) + qJD(4) + t24;
t6 = qJDD(2) + (t75 * qJDD(1) - t65 * t93) * t66;
t5 = (-qJ(3) * qJDD(1) - t93) * t68 - t14;
t4 = t71 * qJD(1) + t98;
t3 = t74 * qJDD(1) + t68 * t93 + qJDD(4) + t14;
t2 = qJDD(2) + (t76 * qJD(1) + t73 * qJDD(1)) * t66;
t1 = -qJD(1) * t99 + t71 * qJDD(1) + t95;
t25 = [t10 * t47 + t20 * t49 + t11 * t107 + t19 * t106 + t42 * t38 + t29 * t23 + t108 * t43 + m(3) * (t13 * t34 + t14 * t35) + m(4) * (-t17 * t43 + t18 * t5 + t19 * t9 + t20 * t6) + m(5) * (t1 * t10 + t11 * t3 + t12 * t2 + t29 * t8 + t4 * t42 + t43 * t7) + (-m(3) * t105 - t70 * mrSges(2,1) + t69 * mrSges(2,2) - t77 * t33 - t119 * t32 - t118 * (t33 * pkin(2) + t105)) * g(2) + (-m(3) * t82 + t69 * mrSges(2,1) + t70 * mrSges(2,2) + t77 * t31 + t119 * t30 + t118 * (t31 * pkin(2) - t82)) * g(1) + (t13 * mrSges(3,1) - t14 * mrSges(3,2) + t9 * mrSges(4,2) + t3 * mrSges(5,2) - t5 * mrSges(4,3) - t1 * mrSges(5,3)) * t68 + (t83 * pkin(1) + (-t56 * mrSges(3,1) - t5 * mrSges(4,1) + t3 * mrSges(5,1) + t6 * mrSges(4,2) + t14 * mrSges(3,3) - t2 * mrSges(5,3) + (m(3) * t24 + t37) * qJD(2)) * t67 + (t9 * mrSges(4,1) + t1 * mrSges(5,1) + t56 * mrSges(3,2) - t2 * mrSges(5,2) - t13 * mrSges(3,3) - t6 * mrSges(4,3) + (-m(4) * t16 - t22) * qJD(3) + (-m(3) * t21 + m(4) * t15 + t109) * qJD(2)) * t65 + (g(1) * t70 + g(2) * t69) * (-m(5) * pkin(3) - mrSges(5,1) - t121)) * t66 + (Ifges(2,3) + (t34 * mrSges(3,1) - t35 * mrSges(3,2) - t18 * mrSges(4,3) - t10 * mrSges(5,3) + (Ifges(5,1) + Ifges(4,1) + Ifges(3,3)) * t68) * t68 + ((-t12 * mrSges(5,2) - t34 * mrSges(3,3) - t20 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2) + Ifges(3,1)) * t116 + (-Ifges(4,4) + Ifges(3,5) + Ifges(5,5)) * t92) * t65 + (-t18 * mrSges(4,1) + t35 * mrSges(3,3) - t12 * mrSges(5,3) + (pkin(1) * mrSges(3,1) + (Ifges(5,2) + Ifges(4,3) + Ifges(3,2)) * t67) * t66 + 0.2e1 * (Ifges(3,4) + Ifges(4,6) - Ifges(5,6)) * t116 + (-Ifges(5,4) - Ifges(4,5) + Ifges(3,6)) * t92) * t67) * t66) * qJDD(1); t49 - t94 * t68 * g(3) + m(4) * t6 + m(5) * t2 + ((t111 * t67 - t120 * t65) * qJDD(1) + ((-t37 - t108) * t67 + (-t38 - t109) * t65 - m(5) * (t4 * t65 + t67 * t7) - m(4) * (t15 * t65 - t17 * t67) - m(3) * (-t21 * t65 + t24 * t67)) * qJD(1) + (-g(1) * t69 + g(2) * t70) * t94) * t66 - t83; -mrSges(5,3) * t96 + t47 + (t1 + t72) * m(5) + (t72 + t9) * m(4) + (-t108 * t68 + (t22 + t23) * t116 - m(4) * (-t16 * t116 - t17 * t68) - m(5) * (-t8 * t116 + t68 * t7)) * qJD(1) + t106; (t23 * t115 + t38 * t68) * qJD(1) + (t3 - g(1) * t33 - g(2) * t31 - g(3) * t116 - (-t8 * t115 - t4 * t68) * qJD(1)) * m(5) + t107;];
tau  = t25;
