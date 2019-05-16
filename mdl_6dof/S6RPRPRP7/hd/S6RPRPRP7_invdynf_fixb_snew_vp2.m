% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRP7
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
% Datum: 2019-05-05 17:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:56:46
% EndTime: 2019-05-05 17:56:50
% DurationCPUTime: 1.38s
% Computational Cost: add. (14268->176), mult. (30976->217), div. (0->0), fcn. (19841->8), ass. (0->83)
t131 = -2 * qJD(4);
t89 = sin(qJ(1));
t92 = cos(qJ(1));
t105 = -t92 * g(1) - t89 * g(2);
t102 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t105;
t117 = qJD(1) * qJD(3);
t88 = sin(qJ(3));
t108 = t88 * t117;
t110 = t89 * g(1) - t92 * g(2);
t94 = qJD(1) ^ 2;
t100 = -t94 * qJ(2) + qJDD(2) - t110;
t127 = -pkin(1) - pkin(7);
t61 = t127 * qJDD(1) + t100;
t91 = cos(qJ(3));
t121 = t88 * g(3) + t91 * t61;
t75 = t91 * qJDD(1) - t108;
t32 = (-t75 - t108) * qJ(4) + (-t88 * t91 * t94 + qJDD(3)) * pkin(3) + t121;
t109 = -t91 * g(3) + t88 * t61;
t74 = -t88 * qJDD(1) - t91 * t117;
t119 = qJD(1) * t91;
t77 = qJD(3) * pkin(3) - qJ(4) * t119;
t84 = t88 ^ 2;
t33 = -t84 * t94 * pkin(3) + t74 * qJ(4) - qJD(3) * t77 + t109;
t120 = qJD(1) * t88;
t85 = sin(pkin(9));
t86 = cos(pkin(9));
t68 = t86 * t119 - t85 * t120;
t130 = t68 * t131 + t86 * t32 - t85 * t33;
t67 = (t85 * t91 + t86 * t88) * qJD(1);
t46 = t67 * pkin(4) - t68 * pkin(8);
t93 = qJD(3) ^ 2;
t17 = -qJDD(3) * pkin(4) - t93 * pkin(8) + t68 * t46 - t130;
t51 = t85 * t74 + t86 * t75;
t87 = sin(qJ(5));
t90 = cos(qJ(5));
t54 = t87 * qJD(3) + t90 * t68;
t26 = -t54 * qJD(5) + t90 * qJDD(3) - t87 * t51;
t53 = t90 * qJD(3) - t87 * t68;
t27 = t53 * qJD(5) + t87 * qJDD(3) + t90 * t51;
t65 = qJD(5) + t67;
t41 = t65 * pkin(5) - t54 * qJ(6);
t42 = t65 * mrSges(7,1) - t54 * mrSges(7,3);
t52 = t53 ^ 2;
t113 = m(7) * (-t26 * pkin(5) - t52 * qJ(6) + t54 * t41 + qJDD(6) + t17) + t54 * t42 + t27 * mrSges(7,2);
t39 = -t65 * mrSges(7,2) + t53 * mrSges(7,3);
t40 = -t65 * mrSges(6,2) + t53 * mrSges(6,3);
t43 = t65 * mrSges(6,1) - t54 * mrSges(6,3);
t129 = -m(6) * t17 - t27 * mrSges(6,2) + (t40 + t39) * t53 + (mrSges(6,1) + mrSges(7,1)) * t26 - t54 * t43 - t113;
t128 = -m(2) - m(3);
t126 = (mrSges(2,1) - mrSges(3,2));
t124 = -mrSges(2,2) + mrSges(3,3);
t112 = t67 * t131 + t85 * t32 + t86 * t33;
t18 = -t93 * pkin(4) + qJDD(3) * pkin(8) - t67 * t46 + t112;
t50 = t86 * t74 - t85 * t75;
t96 = -t74 * pkin(3) + qJDD(4) + t77 * t119 + (-qJ(4) * t84 + t127) * t94 + t102;
t21 = (qJD(3) * t67 - t51) * pkin(8) + (qJD(3) * t68 - t50) * pkin(4) + t96;
t123 = t90 * t18 + t87 * t21;
t107 = -t87 * t18 + t90 * t21;
t49 = qJDD(5) - t50;
t115 = m(7) * (-0.2e1 * qJD(6) * t54 + (t53 * t65 - t27) * qJ(6) + (t53 * t54 + t49) * pkin(5) + t107) + t65 * t39 + t49 * mrSges(7,1);
t36 = -t53 * mrSges(7,1) + t54 * mrSges(7,2);
t114 = m(7) * (-t52 * pkin(5) + t26 * qJ(6) + 0.2e1 * qJD(6) * t53 - t65 * t41 + t123) + t53 * t36 + t26 * mrSges(7,3);
t45 = t67 * mrSges(5,1) + t68 * mrSges(5,2);
t59 = -qJD(3) * mrSges(5,2) - t67 * mrSges(5,3);
t11 = m(5) * t130 + qJDD(3) * mrSges(5,1) - t51 * mrSges(5,3) + qJD(3) * t59 - t68 * t45 + t129;
t37 = -t53 * mrSges(6,1) + t54 * mrSges(6,2);
t10 = m(6) * t123 + t26 * mrSges(6,3) + t53 * t37 + (-t43 - t42) * t65 + (-mrSges(6,2) - mrSges(7,2)) * t49 + t114;
t60 = qJD(3) * mrSges(5,1) - t68 * mrSges(5,3);
t8 = m(6) * t107 + t49 * mrSges(6,1) + t65 * t40 + (-t37 - t36) * t54 + (-mrSges(6,3) - mrSges(7,3)) * t27 + t115;
t6 = m(5) * t112 - qJDD(3) * mrSges(5,2) + t50 * mrSges(5,3) - qJD(3) * t60 + t90 * t10 - t67 * t45 - t87 * t8;
t73 = (mrSges(4,1) * t88 + mrSges(4,2) * t91) * qJD(1);
t76 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t120;
t3 = m(4) * t121 + qJDD(3) * mrSges(4,1) - t75 * mrSges(4,3) + qJD(3) * t76 + t86 * t11 - t73 * t119 + t85 * t6;
t78 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t119;
t4 = m(4) * t109 - qJDD(3) * mrSges(4,2) + t74 * mrSges(4,3) - qJD(3) * t78 - t85 * t11 - t73 * t120 + t86 * t6;
t111 = -t88 * t3 + t91 * t4;
t101 = -m(3) * (-qJDD(1) * pkin(1) + t100) - t91 * t3 - t88 * t4;
t99 = m(5) * t96 - t50 * mrSges(5,1) + t51 * mrSges(5,2) + t87 * t10 + t67 * t59 + t68 * t60 + t90 * t8;
t97 = -t74 * mrSges(4,1) + m(4) * (t127 * t94 + t102) + t76 * t120 + t78 * t119 + t75 * mrSges(4,2) + t99;
t95 = -m(3) * (t94 * pkin(1) - t102) + t97;
t5 = m(2) * t105 + t124 * qJDD(1) - (t126 * t94) + t95;
t1 = m(2) * t110 + t126 * qJDD(1) + t124 * t94 + t101;
t2 = [-m(1) * g(1) - t89 * t1 + t92 * t5, t5, -m(3) * g(3) + t111, t4, t6, t10, -t49 * mrSges(7,2) - t65 * t42 + t114; -m(1) * g(2) + t92 * t1 + t89 * t5, t1, -(t94 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t95, t3, t11, t8, -t27 * mrSges(7,3) - t54 * t36 + t115; (-m(1) + t128) * g(3) + t111, t128 * g(3) + t111, qJDD(1) * mrSges(3,2) - t94 * mrSges(3,3) - t101, t97, t99, -t129, -t26 * mrSges(7,1) - t53 * t39 + t113;];
f_new  = t2;
