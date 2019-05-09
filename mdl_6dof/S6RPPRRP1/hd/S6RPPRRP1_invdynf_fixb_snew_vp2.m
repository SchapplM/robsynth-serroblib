% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-05-05 14:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:43:38
% EndTime: 2019-05-05 14:43:40
% DurationCPUTime: 1.59s
% Computational Cost: add. (19425->164), mult. (43026->203), div. (0->0), fcn. (29476->10), ass. (0->85)
t92 = qJD(1) ^ 2;
t83 = cos(pkin(10));
t79 = t83 ^ 2;
t81 = sin(pkin(10));
t116 = t81 ^ 2 + t79;
t125 = t116 * mrSges(4,3);
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t98 = t81 * t86 - t83 * t89;
t64 = t98 * qJD(1);
t99 = t81 * t89 + t83 * t86;
t65 = t99 * qJD(1);
t113 = t65 * qJD(4);
t56 = -t98 * qJDD(1) - t113;
t115 = pkin(7) * qJDD(1);
t112 = qJD(1) * qJD(3);
t80 = -g(3) + qJDD(2);
t117 = -0.2e1 * t81 * t112 + t83 * t80;
t123 = pkin(3) * t92;
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t106 = t87 * g(1) - t90 * g(2);
t70 = qJDD(1) * pkin(1) + t106;
t102 = -t90 * g(1) - t87 * g(2);
t71 = -t92 * pkin(1) + t102;
t82 = sin(pkin(9));
t84 = cos(pkin(9));
t118 = t82 * t70 + t84 * t71;
t49 = -t92 * pkin(2) + qJDD(1) * qJ(3) + t118;
t32 = (t83 * t123 - t115 - t49) * t81 + t117;
t107 = t81 * t80 + (0.2e1 * t112 + t49) * t83;
t35 = t83 * t115 - t79 * t123 + t107;
t104 = t89 * t32 - t86 * t35;
t55 = t64 * pkin(4) - t65 * pkin(8);
t91 = qJD(4) ^ 2;
t20 = -qJDD(4) * pkin(4) - t91 * pkin(8) + t65 * t55 - t104;
t114 = t64 * qJD(4);
t57 = t99 * qJDD(1) - t114;
t85 = sin(qJ(5));
t88 = cos(qJ(5));
t60 = t85 * qJD(4) + t88 * t65;
t33 = -t60 * qJD(5) + t88 * qJDD(4) - t85 * t57;
t59 = t88 * qJD(4) - t85 * t65;
t34 = t59 * qJD(5) + t85 * qJDD(4) + t88 * t57;
t63 = qJD(5) + t64;
t44 = t63 * pkin(5) - t60 * qJ(6);
t45 = t63 * mrSges(7,1) - t60 * mrSges(7,3);
t58 = t59 ^ 2;
t108 = m(7) * (-t33 * pkin(5) - t58 * qJ(6) + t60 * t44 + qJDD(6) + t20) + t34 * mrSges(7,2) + t60 * t45;
t42 = -t63 * mrSges(7,2) + t59 * mrSges(7,3);
t43 = -t63 * mrSges(6,2) + t59 * mrSges(6,3);
t46 = t63 * mrSges(6,1) - t60 * mrSges(6,3);
t124 = m(6) * t20 + t34 * mrSges(6,2) - (t43 + t42) * t59 - (mrSges(6,1) + mrSges(7,1)) * t33 + t60 * t46 + t108;
t120 = t86 * t32 + t89 * t35;
t21 = -t91 * pkin(4) + qJDD(4) * pkin(8) - t64 * t55 + t120;
t103 = t84 * t70 - t82 * t71;
t101 = qJDD(3) - t103;
t93 = (-pkin(3) * t83 - pkin(2)) * qJDD(1) + (-t116 * pkin(7) - qJ(3)) * t92 + t101;
t24 = (-t57 + t114) * pkin(8) + (-t56 + t113) * pkin(4) + t93;
t121 = t88 * t21 + t85 * t24;
t52 = t64 * mrSges(5,1) + t65 * mrSges(5,2);
t61 = -qJD(4) * mrSges(5,2) - t64 * mrSges(5,3);
t14 = m(5) * t104 + qJDD(4) * mrSges(5,1) - t57 * mrSges(5,3) + qJD(4) * t61 - t65 * t52 - t124;
t105 = -t85 * t21 + t88 * t24;
t54 = qJDD(5) - t56;
t110 = m(7) * (-0.2e1 * qJD(6) * t60 + (t59 * t63 - t34) * qJ(6) + (t59 * t60 + t54) * pkin(5) + t105) + t63 * t42 + t54 * mrSges(7,1);
t39 = -t59 * mrSges(7,1) + t60 * mrSges(7,2);
t40 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t11 = m(6) * t105 + t54 * mrSges(6,1) + t63 * t43 + (-t40 - t39) * t60 + (-mrSges(6,3) - mrSges(7,3)) * t34 + t110;
t109 = m(7) * (-t58 * pkin(5) + t33 * qJ(6) + 0.2e1 * qJD(6) * t59 - t63 * t44 + t121) + t59 * t39 + t33 * mrSges(7,3);
t13 = m(6) * t121 + t33 * mrSges(6,3) + t59 * t40 + (-t46 - t45) * t63 + (-mrSges(6,2) - mrSges(7,2)) * t54 + t109;
t62 = qJD(4) * mrSges(5,1) - t65 * mrSges(5,3);
t9 = m(5) * t120 - qJDD(4) * mrSges(5,2) + t56 * mrSges(5,3) - qJD(4) * t62 - t85 * t11 + t88 * t13 - t64 * t52;
t100 = -t83 * mrSges(4,1) + t81 * mrSges(4,2);
t97 = qJDD(1) * mrSges(4,3) + t92 * t100;
t6 = m(4) * t117 + t86 * t9 + t89 * t14 + (-m(4) * t49 - t97) * t81;
t7 = m(4) * t107 - t86 * t14 + t97 * t83 + t89 * t9;
t111 = m(3) * t80 + t83 * t6 + t81 * t7;
t96 = -m(5) * t93 + t56 * mrSges(5,1) - t57 * mrSges(5,2) - t88 * t11 - t85 * t13 - t64 * t61 - t65 * t62;
t94 = m(4) * (-qJDD(1) * pkin(2) - t92 * qJ(3) + t101) - t96;
t8 = m(3) * t103 + (-mrSges(3,2) + t125) * t92 + (mrSges(3,1) - t100) * qJDD(1) - t94;
t3 = m(3) * t118 - t92 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t81 * t6 + t83 * t7;
t2 = m(2) * t102 - t92 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t84 * t3 - t82 * t8;
t1 = m(2) * t106 + qJDD(1) * mrSges(2,1) - t92 * mrSges(2,2) + t82 * t3 + t84 * t8;
t4 = [-m(1) * g(1) - t87 * t1 + t90 * t2, t2, t3, t7, t9, t13, -t54 * mrSges(7,2) - t63 * t45 + t109; -m(1) * g(2) + t90 * t1 + t87 * t2, t1, t8, t6, t14, t11, -t34 * mrSges(7,3) - t60 * t39 + t110; (-m(1) - m(2)) * g(3) + t111, -m(2) * g(3) + t111, t111, t100 * qJDD(1) - t92 * t125 + t94, -t96, t124, -t33 * mrSges(7,1) - t59 * t42 + t108;];
f_new  = t4;
