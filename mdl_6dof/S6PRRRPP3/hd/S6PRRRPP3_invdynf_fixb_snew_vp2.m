% Calculate vector of cutting forces with Newton-Euler
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-05-05 06:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRPP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:49:44
% EndTime: 2019-05-05 06:49:47
% DurationCPUTime: 1.25s
% Computational Cost: add. (14026->177), mult. (26436->208), div. (0->0), fcn. (17130->10), ass. (0->89)
t128 = cos(qJ(4));
t88 = sin(qJ(3));
t114 = qJD(2) * t88;
t87 = sin(qJ(4));
t66 = t87 * qJD(3) + t114 * t128;
t112 = qJD(2) * qJD(3);
t90 = cos(qJ(3));
t105 = t90 * t112;
t69 = t88 * qJDD(2) + t105;
t40 = t66 * qJD(4) - qJDD(3) * t128 + t87 * t69;
t65 = -qJD(3) * t128 + t114 * t87;
t41 = -t65 * qJD(4) + t87 * qJDD(3) + t128 * t69;
t113 = t90 * qJD(2);
t77 = -qJD(4) + t113;
t54 = t66 * mrSges(6,1) - t77 * mrSges(6,2);
t127 = t65 * t77;
t131 = -2 * qJD(5);
t84 = sin(pkin(6));
t89 = sin(qJ(2));
t125 = t84 * t89;
t83 = sin(pkin(10));
t85 = cos(pkin(10));
t71 = t83 * g(1) - t85 * g(2);
t86 = cos(pkin(6));
t126 = t71 * t86;
t72 = -t85 * g(1) - t83 * g(2);
t82 = -g(3) + qJDD(1);
t91 = cos(qJ(2));
t108 = t82 * t125 + t89 * t126 + t91 * t72;
t93 = qJD(2) ^ 2;
t32 = -t93 * pkin(2) + qJDD(2) * pkin(8) + t108;
t56 = -t84 * t71 + t86 * t82;
t104 = -t88 * t32 + t90 * t56;
t68 = (-pkin(3) * t90 - pkin(9) * t88) * qJD(2);
t92 = qJD(3) ^ 2;
t25 = -qJDD(3) * pkin(3) - t92 * pkin(9) + t68 * t114 - t104;
t94 = (-t41 - t127) * qJ(5) + t25 + (-t77 * pkin(4) + t131) * t66;
t135 = m(6) * (t40 * pkin(4) + t94) - t41 * mrSges(6,3) - t66 * t54;
t134 = (t82 * t84 + t126) * t91 - t89 * t72;
t130 = 2 * qJD(6);
t119 = t90 * t32 + t88 * t56;
t26 = -t92 * pkin(3) + qJDD(3) * pkin(9) + t113 * t68 + t119;
t106 = t88 * t112;
t31 = -qJDD(2) * pkin(2) - t93 * pkin(8) - t134;
t70 = t90 * qJDD(2) - t106;
t28 = (-t69 - t105) * pkin(9) + (-t70 + t106) * pkin(3) + t31;
t103 = t128 * t28 - t87 * t26;
t44 = t65 * pkin(4) - t66 * qJ(5);
t62 = qJDD(4) - t70;
t76 = t77 ^ 2;
t20 = -t62 * pkin(4) - t76 * qJ(5) + t66 * t44 + qJDD(5) - t103;
t43 = -t66 * mrSges(7,2) + t65 * mrSges(7,3);
t111 = m(7) * (t77 * t130 + (t65 * t66 - t62) * qJ(6) + (t41 - t127) * pkin(5) + t20) + t66 * t43 + t41 * mrSges(7,1);
t100 = m(6) * t20 + t111;
t52 = t65 * mrSges(6,1) + t77 * mrSges(6,3);
t53 = -t65 * mrSges(7,1) - t77 * mrSges(7,2);
t115 = -t52 + t53;
t46 = -t65 * mrSges(6,2) - t66 * mrSges(6,3);
t118 = -t65 * mrSges(5,1) - t66 * mrSges(5,2) - t46;
t121 = -mrSges(5,3) - mrSges(6,1);
t122 = mrSges(6,2) - mrSges(7,3);
t48 = t77 * mrSges(5,2) - t65 * mrSges(5,3);
t11 = m(5) * t103 + t118 * t66 + t121 * t41 + (-t48 - t115) * t77 + (mrSges(5,1) - t122) * t62 - t100;
t120 = t128 * t26 + t87 * t28;
t49 = -t77 * mrSges(5,1) - t66 * mrSges(5,3);
t50 = t66 * pkin(5) + t77 * qJ(6);
t51 = t66 * mrSges(7,1) + t77 * mrSges(7,3);
t61 = t65 ^ 2;
t96 = -t76 * pkin(4) + t62 * qJ(5) - t65 * t44 + t120;
t109 = t77 * t51 - m(7) * (-t40 * pkin(5) - t61 * qJ(6) + qJDD(6) + (t131 - t50) * t77 + t96) - t62 * mrSges(7,2);
t99 = m(6) * (0.2e1 * qJD(5) * t77 - t96) + t109;
t12 = m(5) * t120 + (t49 - t54) * t77 + (-mrSges(5,2) + mrSges(6,3)) * t62 + (-t43 + t118) * t65 + (-mrSges(7,1) + t121) * t40 - t99;
t73 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t114;
t74 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t113;
t133 = m(4) * t31 - t70 * mrSges(4,1) + t69 * mrSges(4,2) + (t73 * t88 - t74 * t90) * qJD(2) + t11 * t128 + t87 * t12;
t110 = m(7) * (-t61 * pkin(5) + t65 * t130 - t66 * t50 + (pkin(4) + qJ(6)) * t40 + t94) + t40 * mrSges(7,3) + t65 * t53;
t132 = m(5) * t25 + (t49 - t51) * t66 + (t48 - t52) * t65 + (mrSges(5,2) - mrSges(7,2)) * t41 + (mrSges(5,1) - mrSges(6,2)) * t40 + t110 + t135;
t8 = m(3) * t134 + qJDD(2) * mrSges(3,1) - t93 * mrSges(3,2) - t133;
t129 = t8 * t91;
t67 = (-mrSges(4,1) * t90 + mrSges(4,2) * t88) * qJD(2);
t10 = m(4) * t104 + qJDD(3) * mrSges(4,1) - t69 * mrSges(4,3) + qJD(3) * t74 - t114 * t67 - t132;
t9 = m(4) * t119 - qJDD(3) * mrSges(4,2) + t70 * mrSges(4,3) - qJD(3) * t73 - t87 * t11 + t113 * t67 + t12 * t128;
t4 = m(3) * t108 - t93 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t88 * t10 + t90 * t9;
t6 = m(3) * t56 + t90 * t10 + t88 * t9;
t107 = m(2) * t82 + t4 * t125 + t84 * t129 + t86 * t6;
t98 = -t41 * mrSges(7,2) - t66 * t51 + t110;
t2 = m(2) * t72 + t91 * t4 - t89 * t8;
t1 = m(2) * t71 - t84 * t6 + (t4 * t89 + t129) * t86;
t3 = [-m(1) * g(1) - t83 * t1 + t85 * t2, t2, t4, t9, t12, -t40 * mrSges(6,2) - t65 * t52 + t135 + t98, t98; -m(1) * g(2) + t85 * t1 + t83 * t2, t1, t8, t10, t11, -t62 * mrSges(6,3) + t77 * t54 + (t43 + t46) * t65 + (mrSges(6,1) + mrSges(7,1)) * t40 + t99, -t62 * mrSges(7,3) + t77 * t53 + t111; -m(1) * g(3) + t107, t107, t6, t133, t132, t41 * mrSges(6,1) + t115 * t77 + t122 * t62 + t66 * t46 + t100, -t40 * mrSges(7,1) - t65 * t43 - t109;];
f_new  = t3;
