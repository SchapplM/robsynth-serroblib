% Calculate vector of cutting forces with Newton-Euler
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-05-05 16:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:32:44
% EndTime: 2019-05-05 16:32:48
% DurationCPUTime: 1.55s
% Computational Cost: add. (16167->179), mult. (35105->225), div. (0->0), fcn. (21888->10), ass. (0->91)
t92 = cos(qJ(3));
t116 = qJD(1) * t92;
t89 = sin(qJ(3));
t117 = qJD(1) * t89;
t118 = cos(pkin(10));
t85 = sin(pkin(10));
t63 = -t118 * t116 + t85 * t117;
t64 = (t118 * t89 + t85 * t92) * qJD(1);
t42 = t63 * pkin(4) - t64 * qJ(5);
t134 = (2 * qJD(4)) + t42;
t90 = sin(qJ(1));
t93 = cos(qJ(1));
t112 = t90 * g(1) - t93 * g(2);
t70 = qJDD(1) * pkin(1) + t112;
t108 = -t93 * g(1) - t90 * g(2);
t95 = qJD(1) ^ 2;
t72 = -t95 * pkin(1) + t108;
t86 = sin(pkin(9));
t87 = cos(pkin(9));
t109 = t87 * t70 - t86 * t72;
t105 = -qJDD(1) * pkin(2) - t109;
t114 = qJD(1) * qJD(3);
t111 = t92 * t114;
t73 = t89 * qJDD(1) + t111;
t74 = t92 * qJDD(1) - t89 * t114;
t76 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t117;
t77 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t116;
t115 = qJD(3) * t63;
t119 = t86 * t70 + t87 * t72;
t39 = -t95 * pkin(2) + qJDD(1) * pkin(7) + t119;
t84 = -g(3) + qJDD(2);
t110 = -t89 * t39 + t92 * t84;
t25 = (-t73 + t111) * qJ(4) + (t89 * t92 * t95 + qJDD(3)) * pkin(3) + t110;
t122 = t92 * t39 + t89 * t84;
t75 = qJD(3) * pkin(3) - qJ(4) * t117;
t83 = t92 ^ 2;
t26 = -t83 * t95 * pkin(3) + t74 * qJ(4) - qJD(3) * t75 + t122;
t106 = t118 * t25 - t85 * t26;
t94 = qJD(3) ^ 2;
t19 = -qJDD(3) * pkin(4) - t94 * qJ(5) + t134 * t64 + qJDD(5) - t106;
t49 = t118 * t73 + t85 * t74;
t14 = (t63 * t64 - qJDD(3)) * pkin(8) + (t49 + t115) * pkin(5) + t19;
t48 = -t118 * t74 + t85 * t73;
t56 = t64 * pkin(5) - qJD(3) * pkin(8);
t62 = t63 ^ 2;
t129 = -2 * qJD(5);
t98 = -t74 * pkin(3) + qJDD(4) + t75 * t117 + (-qJ(4) * t83 - pkin(7)) * t95 + t105;
t96 = (-t49 + t115) * qJ(5) + t98 + (pkin(4) * qJD(3) + t129) * t64;
t17 = -t62 * pkin(5) - t64 * t56 + (pkin(4) + pkin(8)) * t48 + t96;
t88 = sin(qJ(6));
t91 = cos(qJ(6));
t50 = -t88 * qJD(3) + t91 * t63;
t32 = t50 * qJD(6) + t91 * qJDD(3) + t88 * t48;
t51 = t91 * qJD(3) + t88 * t63;
t33 = -mrSges(7,1) * t50 + mrSges(7,2) * t51;
t61 = qJD(6) + t64;
t35 = -t61 * mrSges(7,2) + t50 * mrSges(7,3);
t47 = qJDD(6) + t49;
t12 = m(7) * (t91 * t14 - t88 * t17) - t32 * mrSges(7,3) + t47 * mrSges(7,1) - t51 * t33 + t61 * t35;
t31 = -t51 * qJD(6) - t88 * qJDD(3) + t91 * t48;
t36 = t61 * mrSges(7,1) - t51 * mrSges(7,3);
t13 = m(7) * (t88 * t14 + t91 * t17) + t31 * mrSges(7,3) - t47 * mrSges(7,2) + t50 * t33 - t61 * t36;
t55 = t64 * mrSges(6,1) + qJD(3) * mrSges(6,2);
t104 = t88 * t12 - t91 * t13 - m(6) * (t48 * pkin(4) + t96) + t64 * t55 + t49 * mrSges(6,3);
t54 = t63 * mrSges(6,1) - qJD(3) * mrSges(6,3);
t120 = -qJD(3) * mrSges(5,2) - t63 * mrSges(5,3) - t54;
t125 = mrSges(5,1) - mrSges(6,2);
t53 = qJD(3) * mrSges(5,1) - t64 * mrSges(5,3);
t97 = m(5) * t98 + t49 * mrSges(5,2) + t120 * t63 + t125 * t48 + t64 * t53 - t104;
t133 = (t89 * t76 - t92 * t77) * qJD(1) - t74 * mrSges(4,1) + t73 * mrSges(4,2) + m(4) * (-t95 * pkin(7) + t105) + t97;
t131 = -2 * qJD(4);
t124 = -mrSges(5,3) - mrSges(6,1);
t123 = t118 * t26 + t85 * t25;
t44 = -t63 * mrSges(6,2) - t64 * mrSges(6,3);
t121 = -t63 * mrSges(5,1) - t64 * mrSges(5,2) - t44;
t59 = t63 * t131;
t102 = t94 * pkin(4) - qJDD(3) * qJ(5) - t123;
t101 = -t31 * mrSges(7,1) - t50 * t35 + m(7) * (-t48 * pkin(5) - t62 * pkin(8) - t63 * t42 + t59 + ((2 * qJD(5)) + t56) * qJD(3) - t102) + t51 * t36 + t32 * mrSges(7,2);
t99 = -m(6) * (qJD(3) * t129 + t134 * t63 + t102) + t101;
t10 = m(5) * (t59 + t123) + t121 * t63 + t124 * t48 + (-mrSges(5,2) + mrSges(6,3)) * qJDD(3) + (-t53 + t55) * qJD(3) + t99;
t71 = (-mrSges(4,1) * t92 + mrSges(4,2) * t89) * qJD(1);
t103 = -m(6) * t19 - t91 * t12 - t88 * t13;
t9 = m(5) * t106 + (m(5) * t131 + t121) * t64 + t124 * t49 + t125 * qJDD(3) + t120 * qJD(3) + t103;
t6 = m(4) * t110 + qJDD(3) * mrSges(4,1) - t73 * mrSges(4,3) + qJD(3) * t77 + t85 * t10 - t71 * t117 + t118 * t9;
t7 = m(4) * t122 - qJDD(3) * mrSges(4,2) + t74 * mrSges(4,3) - qJD(3) * t76 + t118 * t10 + t71 * t116 - t85 * t9;
t113 = m(3) * t84 + t92 * t6 + t89 * t7;
t8 = m(3) * t109 + qJDD(1) * mrSges(3,1) - t95 * mrSges(3,2) - t133;
t3 = m(3) * t119 - t95 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t89 * t6 + t92 * t7;
t2 = m(2) * t108 - t95 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t87 * t3 - t86 * t8;
t1 = m(2) * t112 + qJDD(1) * mrSges(2,1) - t95 * mrSges(2,2) + t86 * t3 + t87 * t8;
t4 = [-m(1) * g(1) - t90 * t1 + t93 * t2, t2, t3, t7, t10, -t48 * mrSges(6,2) - t63 * t54 - t104, t13; -m(1) * g(2) + t93 * t1 + t90 * t2, t1, t8, t6, t9, t48 * mrSges(6,1) - qJDD(3) * mrSges(6,3) - qJD(3) * t55 + t63 * t44 - t99, t12; (-m(1) - m(2)) * g(3) + t113, -m(2) * g(3) + t113, t113, t133, t97, t49 * mrSges(6,1) + qJDD(3) * mrSges(6,2) + qJD(3) * t54 + t64 * t44 - t103, t101;];
f_new  = t4;
