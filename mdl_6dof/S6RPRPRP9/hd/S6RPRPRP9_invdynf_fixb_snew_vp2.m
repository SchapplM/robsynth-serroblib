% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 18:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRP9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:04:49
% EndTime: 2019-05-05 18:04:53
% DurationCPUTime: 1.40s
% Computational Cost: add. (16004->173), mult. (33239->214), div. (0->0), fcn. (20594->8), ass. (0->85)
t85 = sin(qJ(1));
t87 = cos(qJ(1));
t102 = -t87 * g(1) - t85 * g(2);
t124 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t102;
t123 = -m(2) - m(3);
t122 = (-pkin(1) - pkin(7));
t120 = cos(qJ(5));
t86 = cos(qJ(3));
t114 = qJD(1) * t86;
t81 = sin(pkin(9));
t82 = cos(pkin(9));
t66 = t82 * qJD(3) - t81 * t114;
t67 = t81 * qJD(3) + t82 * t114;
t83 = sin(qJ(5));
t43 = -t120 * t66 + t83 * t67;
t44 = t120 * t67 + t83 * t66;
t28 = t43 * pkin(5) - t44 * qJ(6);
t112 = qJD(1) * qJD(3);
t104 = t86 * t112;
t84 = sin(qJ(3));
t72 = t84 * qJDD(1) + t104;
t69 = qJDD(5) + t72;
t113 = t84 * qJD(1);
t77 = qJD(5) + t113;
t76 = t77 ^ 2;
t105 = t84 * t112;
t73 = t86 * qJDD(1) - t105;
t89 = qJD(1) ^ 2;
t93 = (t122 * t89) - t124;
t33 = (-t73 + t105) * qJ(4) + (t72 + t104) * pkin(3) + t93;
t108 = t85 * g(1) - t87 * g(2);
t97 = -t89 * qJ(2) + qJDD(2) - t108;
t55 = t122 * qJDD(1) + t97;
t107 = -t86 * g(3) + t84 * t55;
t70 = (pkin(3) * t84 - qJ(4) * t86) * qJD(1);
t88 = qJD(3) ^ 2;
t37 = -t88 * pkin(3) + qJDD(3) * qJ(4) - t70 * t113 + t107;
t103 = -0.2e1 * qJD(4) * t67 + t82 * t33 - t81 * t37;
t53 = t81 * qJDD(3) + t82 * t73;
t17 = (t66 * t113 - t53) * pkin(8) + (t66 * t67 + t72) * pkin(4) + t103;
t109 = 0.2e1 * qJD(4) * t66 + t81 * t33 + t82 * t37;
t52 = t82 * qJDD(3) - t81 * t73;
t54 = pkin(4) * t113 - t67 * pkin(8);
t65 = t66 ^ 2;
t19 = -t65 * pkin(4) + t52 * pkin(8) - t54 * t113 + t109;
t99 = t120 * t17 - t83 * t19;
t121 = m(7) * (-t69 * pkin(5) - t76 * qJ(6) + t44 * t28 + qJDD(6) - t99);
t119 = (mrSges(2,1) - mrSges(3,2));
t118 = -mrSges(2,2) + mrSges(3,3);
t117 = -mrSges(6,3) - mrSges(7,2);
t116 = t120 * t19 + t83 * t17;
t29 = t43 * mrSges(7,1) - t44 * mrSges(7,3);
t115 = -t43 * mrSges(6,1) - t44 * mrSges(6,2) - t29;
t41 = -t77 * mrSges(7,1) + t44 * mrSges(7,2);
t110 = m(7) * (-t76 * pkin(5) + t69 * qJ(6) + 0.2e1 * qJD(6) * t77 - t43 * t28 + t116) + t77 * t41 + t69 * mrSges(7,3);
t101 = t84 * g(3) + t86 * t55;
t71 = (mrSges(4,1) * t84 + mrSges(4,2) * t86) * qJD(1);
t74 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t113;
t35 = -qJDD(3) * pkin(3) - t88 * qJ(4) + t70 * t114 + qJDD(4) - t101;
t50 = -mrSges(5,2) * t113 + t66 * mrSges(5,3);
t51 = mrSges(5,1) * t113 - t67 * mrSges(5,3);
t24 = t44 * qJD(5) - t120 * t52 + t83 * t53;
t25 = -t43 * qJD(5) + t120 * t53 + t83 * t52;
t38 = -t77 * mrSges(6,2) - t43 * mrSges(6,3);
t40 = t77 * mrSges(6,1) - t44 * mrSges(6,3);
t91 = -t52 * pkin(4) - t65 * pkin(8) + t67 * t54 + t35;
t39 = -t43 * mrSges(7,2) + t77 * mrSges(7,3);
t96 = t25 * mrSges(7,3) + t44 * t41 - m(7) * (-0.2e1 * qJD(6) * t44 + (t43 * t77 - t25) * qJ(6) + (t44 * t77 + t24) * pkin(5) + t91) - t24 * mrSges(7,1) - t43 * t39;
t92 = m(6) * t91 + t24 * mrSges(6,1) + t25 * mrSges(6,2) + t43 * t38 + t44 * t40 - t96;
t90 = m(5) * t35 - t52 * mrSges(5,1) + t53 * mrSges(5,2) - t66 * t50 + t67 * t51 + t92;
t11 = m(4) * t101 + qJDD(3) * mrSges(4,1) - t73 * mrSges(4,3) + qJD(3) * t74 - t71 * t114 - t90;
t10 = m(6) * t99 - t121 + (t38 + t39) * t77 + (mrSges(6,1) + mrSges(7,1)) * t69 + t115 * t44 + t117 * t25;
t45 = -t66 * mrSges(5,1) + t67 * mrSges(5,2);
t9 = m(6) * t116 - t69 * mrSges(6,2) + t115 * t43 + t117 * t24 - t77 * t40 + t110;
t7 = m(5) * t103 + t72 * mrSges(5,1) - t53 * mrSges(5,3) + t120 * t10 + t50 * t113 - t67 * t45 + t83 * t9;
t75 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t114;
t8 = m(5) * t109 - t72 * mrSges(5,2) + t52 * mrSges(5,3) - t83 * t10 - t51 * t113 + t120 * t9 + t66 * t45;
t4 = m(4) * t107 - qJDD(3) * mrSges(4,2) - t72 * mrSges(4,3) - qJD(3) * t75 - t71 * t113 - t81 * t7 + t82 * t8;
t106 = -t84 * t11 + t86 * t4;
t98 = -m(3) * (-qJDD(1) * pkin(1) + t97) - t86 * t11 - t84 * t4;
t95 = m(4) * t93 + t72 * mrSges(4,1) + t73 * mrSges(4,2) + t74 * t113 + t75 * t114 + t82 * t7 + t81 * t8;
t94 = -m(3) * (t89 * pkin(1) + t124) + t95;
t2 = m(2) * t102 + t118 * qJDD(1) - (t119 * t89) + t94;
t1 = m(2) * t108 + t119 * qJDD(1) + t118 * t89 + t98;
t3 = [-m(1) * g(1) - t85 * t1 + t87 * t2, t2, -m(3) * g(3) + t106, t4, t8, t9, -t24 * mrSges(7,2) - t43 * t29 + t110; -m(1) * g(2) + t87 * t1 + t85 * t2, t1, -(t89 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t94, t11, t7, t10, -t96; (-m(1) + t123) * g(3) + t106, t123 * g(3) + t106, qJDD(1) * mrSges(3,2) - t89 * mrSges(3,3) - t98, t95, t90, t92, -t69 * mrSges(7,1) + t25 * mrSges(7,2) + t44 * t29 - t77 * t39 + t121;];
f_new  = t3;
