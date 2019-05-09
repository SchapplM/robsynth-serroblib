% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:57:04
% EndTime: 2019-05-05 15:57:06
% DurationCPUTime: 0.97s
% Computational Cost: add. (11606->153), mult. (21883->185), div. (0->0), fcn. (12547->8), ass. (0->80)
t82 = sin(qJ(1));
t86 = cos(qJ(1));
t112 = t82 * g(1) - t86 * g(2);
t88 = qJD(1) ^ 2;
t48 = -qJDD(1) * pkin(1) - t88 * qJ(2) + qJDD(2) - t112;
t118 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t48;
t100 = -t86 * g(1) - t82 * g(2);
t117 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t100;
t116 = -m(3) - m(4);
t115 = mrSges(3,2) - mrSges(4,3);
t114 = -mrSges(4,2) - mrSges(3,3);
t109 = qJD(1) * qJD(4);
t81 = sin(qJ(4));
t103 = t81 * t109;
t85 = cos(qJ(4));
t67 = t85 * t109;
t61 = -t81 * qJDD(1) - t67;
t62 = t85 * qJDD(1) - t103;
t92 = -t88 * pkin(7) - t118;
t25 = (-t62 + t103) * pkin(8) + (-t61 + t67) * pkin(4) + t92;
t90 = qJDD(3) + (-pkin(1) - qJ(3)) * t88 + t117;
t39 = -qJDD(1) * pkin(7) + t90;
t105 = -t85 * g(3) + t81 * t39;
t60 = (pkin(4) * t81 - pkin(8) * t85) * qJD(1);
t69 = t81 * qJD(1);
t87 = qJD(4) ^ 2;
t28 = -t87 * pkin(4) + qJDD(4) * pkin(8) - t60 * t69 + t105;
t80 = sin(qJ(5));
t84 = cos(qJ(5));
t113 = t80 * t25 + t84 * t28;
t111 = qJD(1) * t85;
t66 = t69 + qJD(5);
t110 = -m(2) + t116;
t106 = mrSges(2,1) - t115;
t56 = qJDD(5) - t61;
t59 = (mrSges(5,1) * t81 + mrSges(5,2) * t85) * qJD(1);
t63 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t69;
t99 = t81 * g(3) + t85 * t39;
t27 = -qJDD(4) * pkin(4) - t87 * pkin(8) + t60 * t111 - t99;
t58 = t80 * qJD(4) + t84 * t111;
t31 = -t58 * qJD(5) + t84 * qJDD(4) - t80 * t62;
t57 = t84 * qJD(4) - t80 * t111;
t32 = t57 * qJD(5) + t80 * qJDD(4) + t84 * t62;
t44 = -t66 * mrSges(6,2) + t57 * mrSges(6,3);
t45 = t66 * mrSges(6,1) - t58 * mrSges(6,3);
t79 = sin(qJ(6));
t83 = cos(qJ(6));
t34 = t79 * t57 + t83 * t58;
t19 = -t34 * qJD(6) + t83 * t31 - t79 * t32;
t33 = t83 * t57 - t79 * t58;
t20 = t33 * qJD(6) + t79 * t31 + t83 * t32;
t65 = qJD(6) + t66;
t29 = -t65 * mrSges(7,2) + t33 * mrSges(7,3);
t30 = t65 * mrSges(7,1) - t34 * mrSges(7,3);
t46 = t66 * pkin(5) - t58 * pkin(9);
t55 = t57 ^ 2;
t93 = t19 * mrSges(7,1) + t33 * t29 - m(7) * (-t31 * pkin(5) - t55 * pkin(9) + t58 * t46 + t27) - t20 * mrSges(7,2) - t34 * t30;
t89 = m(6) * t27 - t31 * mrSges(6,1) + t32 * mrSges(6,2) - t57 * t44 + t58 * t45 - t93;
t11 = m(5) * t99 + qJDD(4) * mrSges(5,1) - t62 * mrSges(5,3) + qJD(4) * t63 - t59 * t111 - t89;
t64 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t111;
t102 = t84 * t25 - t80 * t28;
t14 = (t57 * t66 - t32) * pkin(9) + (t57 * t58 + t56) * pkin(5) + t102;
t15 = -t55 * pkin(5) + t31 * pkin(9) - t66 * t46 + t113;
t24 = -t33 * mrSges(7,1) + t34 * mrSges(7,2);
t50 = qJDD(6) + t56;
t12 = m(7) * (t83 * t14 - t79 * t15) - t20 * mrSges(7,3) + t50 * mrSges(7,1) - t34 * t24 + t65 * t29;
t13 = m(7) * (t79 * t14 + t83 * t15) + t19 * mrSges(7,3) - t50 * mrSges(7,2) + t33 * t24 - t65 * t30;
t36 = -t57 * mrSges(6,1) + t58 * mrSges(6,2);
t8 = m(6) * t102 + t56 * mrSges(6,1) - t32 * mrSges(6,3) + t83 * t12 + t79 * t13 - t58 * t36 + t66 * t44;
t9 = m(6) * t113 - t56 * mrSges(6,2) + t31 * mrSges(6,3) - t79 * t12 + t83 * t13 + t57 * t36 - t66 * t45;
t5 = m(5) * t105 - qJDD(4) * mrSges(5,2) + t61 * mrSges(5,3) - qJD(4) * t64 - t59 * t69 - t80 * t8 + t84 * t9;
t104 = -t81 * t11 + t85 * t5;
t101 = m(4) * t90 + qJDD(1) * mrSges(4,2) + t85 * t11 + t81 * t5;
t97 = m(3) * (t88 * pkin(1) - t117) - t101;
t95 = m(5) * t92 - t61 * mrSges(5,1) + t62 * mrSges(5,2) + t64 * t111 + t63 * t69 + t84 * t8 + t80 * t9;
t94 = m(4) * t118 - t95;
t91 = m(3) * t48 + t94;
t2 = m(2) * t112 + (-mrSges(2,2) - t114) * t88 + t106 * qJDD(1) - t91;
t1 = m(2) * t100 + (-mrSges(2,2) + mrSges(3,3)) * qJDD(1) - t106 * t88 - t97;
t3 = [-m(1) * g(1) + t86 * t1 - t82 * t2, t1, t116 * g(3) + t104, -m(4) * g(3) + t104, t5, t9, t13; -m(1) * g(2) + t82 * t1 + t86 * t2, t2, -qJDD(1) * mrSges(3,3) - t115 * t88 + t97, -t88 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t94, t11, t8, t12; (-m(1) + t110) * g(3) + t104, t110 * g(3) + t104, t115 * qJDD(1) + t114 * t88 + t91, -t88 * mrSges(4,3) + t101, t95, t89, -t93;];
f_new  = t3;
