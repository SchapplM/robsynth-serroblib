% Calculate vector of cutting forces with Newton-Euler
% S6PRRRPR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-05-05 07:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:27:33
% EndTime: 2019-05-05 07:27:39
% DurationCPUTime: 2.22s
% Computational Cost: add. (26673->176), mult. (53042->227), div. (0->0), fcn. (36724->12), ass. (0->95)
t100 = qJD(2) ^ 2;
t89 = sin(pkin(11));
t91 = cos(pkin(11));
t76 = t89 * g(1) - t91 * g(2);
t92 = cos(pkin(6));
t128 = t76 * t92;
t77 = -t91 * g(1) - t89 * g(2);
t88 = -g(3) + qJDD(1);
t90 = sin(pkin(6));
t96 = sin(qJ(2));
t99 = cos(qJ(2));
t137 = (t88 * t90 + t128) * t99 - t96 * t77;
t104 = -qJDD(2) * pkin(2) - t137;
t95 = sin(qJ(3));
t120 = qJD(2) * t95;
t118 = qJD(2) * qJD(3);
t98 = cos(qJ(3));
t75 = t98 * qJDD(2) - t95 * t118;
t81 = qJD(3) * pkin(3) - pkin(9) * t120;
t87 = t98 ^ 2;
t102 = -t75 * pkin(3) + t81 * t120 + (-pkin(9) * t87 - pkin(8)) * t100 + t104;
t119 = qJD(2) * t98;
t133 = cos(qJ(4));
t94 = sin(qJ(4));
t67 = -t133 * t119 + t94 * t120;
t86 = qJD(3) + qJD(4);
t131 = t67 * t86;
t135 = -2 * qJD(5);
t115 = t98 * t118;
t74 = t95 * qJDD(2) + t115;
t42 = -t67 * qJD(4) + t133 * t74 + t94 * t75;
t68 = (t133 * t95 + t94 * t98) * qJD(2);
t101 = (-t42 + t131) * qJ(5) + t102 + (t86 * pkin(4) + t135) * t68;
t127 = t90 * t96;
t117 = t88 * t127 + t96 * t128 + t99 * t77;
t46 = -t100 * pkin(2) + qJDD(2) * pkin(8) + t117;
t61 = -t90 * t76 + t92 * t88;
t114 = -t95 * t46 + t98 * t61;
t27 = (-t74 + t115) * pkin(9) + (t100 * t95 * t98 + qJDD(3)) * pkin(3) + t114;
t123 = t98 * t46 + t95 * t61;
t28 = -t87 * t100 * pkin(3) + t75 * pkin(9) - qJD(3) * t81 + t123;
t113 = t133 * t27 - t94 * t28;
t49 = t67 * pkin(4) - t68 * qJ(5);
t84 = t86 ^ 2;
t85 = qJDD(3) + qJDD(4);
t21 = -t85 * pkin(4) - t84 * qJ(5) + t68 * t49 + qJDD(5) - t113;
t16 = (t67 * t68 - t85) * pkin(10) + (t42 + t131) * pkin(5) + t21;
t41 = t68 * qJD(4) - t133 * t75 + t94 * t74;
t60 = t68 * pkin(5) - t86 * pkin(10);
t66 = t67 ^ 2;
t19 = -t66 * pkin(5) - t68 * t60 + (pkin(4) + pkin(10)) * t41 + t101;
t93 = sin(qJ(6));
t97 = cos(qJ(6));
t52 = t97 * t67 - t93 * t86;
t31 = t52 * qJD(6) + t93 * t41 + t97 * t85;
t53 = t93 * t67 + t97 * t86;
t35 = -t52 * mrSges(7,1) + t53 * mrSges(7,2);
t39 = qJDD(6) + t42;
t64 = qJD(6) + t68;
t43 = -t64 * mrSges(7,2) + t52 * mrSges(7,3);
t14 = m(7) * (t97 * t16 - t93 * t19) - t31 * mrSges(7,3) + t39 * mrSges(7,1) - t53 * t35 + t64 * t43;
t30 = -t53 * qJD(6) + t97 * t41 - t93 * t85;
t44 = t64 * mrSges(7,1) - t53 * mrSges(7,3);
t15 = m(7) * (t93 * t16 + t97 * t19) + t30 * mrSges(7,3) - t39 * mrSges(7,2) + t52 * t35 - t64 * t44;
t59 = t68 * mrSges(6,1) + t86 * mrSges(6,2);
t110 = t93 * t14 - t97 * t15 - m(6) * (t41 * pkin(4) + t101) + t42 * mrSges(6,3) + t68 * t59;
t58 = t67 * mrSges(6,1) - t86 * mrSges(6,3);
t121 = -t86 * mrSges(5,2) - t67 * mrSges(5,3) - t58;
t126 = mrSges(5,1) - mrSges(6,2);
t57 = t86 * mrSges(5,1) - t68 * mrSges(5,3);
t103 = m(5) * t102 + t42 * mrSges(5,2) + t121 * t67 + t126 * t41 + t68 * t57 - t110;
t78 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t120;
t79 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t119;
t138 = t103 + (t95 * t78 - t98 * t79) * qJD(2) - t75 * mrSges(4,1) + t74 * mrSges(4,2) + m(4) * (-t100 * pkin(8) + t104);
t10 = m(3) * t137 + qJDD(2) * mrSges(3,1) - t100 * mrSges(3,2) - t138;
t132 = t10 * t99;
t125 = -mrSges(5,3) - mrSges(6,1);
t124 = t133 * t28 + t94 * t27;
t51 = -t67 * mrSges(6,2) - t68 * mrSges(6,3);
t122 = -t67 * mrSges(5,1) - t68 * mrSges(5,2) - t51;
t109 = -m(6) * t21 - t97 * t14 - t93 * t15;
t11 = m(5) * t113 + t121 * t86 + t122 * t68 + t125 * t42 + t126 * t85 + t109;
t106 = -t84 * pkin(4) + t85 * qJ(5) - t67 * t49 + t124;
t108 = -t30 * mrSges(7,1) - t52 * t43 + m(7) * (-t41 * pkin(5) - t66 * pkin(10) + ((2 * qJD(5)) + t60) * t86 + t106) + t31 * mrSges(7,2) + t53 * t44;
t105 = -m(6) * (t86 * t135 - t106) + t108;
t12 = m(5) * t124 + (-t57 + t59) * t86 + (-mrSges(5,2) + mrSges(6,3)) * t85 + t122 * t67 + t125 * t41 + t105;
t73 = (-mrSges(4,1) * t98 + mrSges(4,2) * t95) * qJD(2);
t7 = m(4) * t114 + qJDD(3) * mrSges(4,1) - t74 * mrSges(4,3) + qJD(3) * t79 + t133 * t11 + t94 * t12 - t73 * t120;
t8 = m(4) * t123 - qJDD(3) * mrSges(4,2) + t75 * mrSges(4,3) - qJD(3) * t78 - t94 * t11 + t73 * t119 + t133 * t12;
t4 = m(3) * t117 - t100 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t95 * t7 + t98 * t8;
t6 = m(3) * t61 + t98 * t7 + t95 * t8;
t116 = m(2) * t88 + t4 * t127 + t90 * t132 + t92 * t6;
t2 = m(2) * t77 - t96 * t10 + t99 * t4;
t1 = m(2) * t76 - t90 * t6 + (t4 * t96 + t132) * t92;
t3 = [-m(1) * g(1) - t89 * t1 + t91 * t2, t2, t4, t8, t12, -t41 * mrSges(6,2) - t67 * t58 - t110, t15; -m(1) * g(2) + t91 * t1 + t89 * t2, t1, t10, t7, t11, t41 * mrSges(6,1) - t85 * mrSges(6,3) + t67 * t51 - t86 * t59 - t105, t14; -m(1) * g(3) + t116, t116, t6, t138, t103, t42 * mrSges(6,1) + t85 * mrSges(6,2) + t68 * t51 + t86 * t58 - t109, t108;];
f_new  = t3;
