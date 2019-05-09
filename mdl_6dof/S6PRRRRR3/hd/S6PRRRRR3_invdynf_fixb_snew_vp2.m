% Calculate vector of cutting forces with Newton-Euler
% S6PRRRRR3
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 11:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:01:13
% EndTime: 2019-05-05 11:01:24
% DurationCPUTime: 6.05s
% Computational Cost: add. (90898->178), mult. (178985->238), div. (0->0), fcn. (131473->14), ass. (0->102)
t100 = sin(qJ(2));
t105 = cos(qJ(2));
t92 = sin(pkin(12));
t94 = cos(pkin(12));
t81 = t92 * g(1) - t94 * g(2);
t95 = cos(pkin(6));
t128 = t81 * t95;
t82 = -t94 * g(1) - t92 * g(2);
t91 = -g(3) + qJDD(1);
t93 = sin(pkin(6));
t131 = -t100 * t82 + (t91 * t93 + t128) * t105;
t102 = cos(qJ(5));
t103 = cos(qJ(4));
t106 = qJD(3) ^ 2;
t104 = cos(qJ(3));
t122 = t104 * qJD(2);
t107 = qJD(2) ^ 2;
t124 = t100 * t93;
t120 = t100 * t128 + t105 * t82 + t91 * t124;
t51 = -t107 * pkin(2) + qJDD(2) * pkin(8) + t120;
t66 = -t93 * t81 + t95 * t91;
t99 = sin(qJ(3));
t125 = t104 * t51 + t99 * t66;
t78 = (-pkin(3) * t104 - pkin(9) * t99) * qJD(2);
t38 = -t106 * pkin(3) + qJDD(3) * pkin(9) + t78 * t122 + t125;
t121 = qJD(2) * qJD(3);
t116 = t104 * t121;
t50 = -qJDD(2) * pkin(2) - t107 * pkin(8) - t131;
t79 = t99 * qJDD(2) + t116;
t89 = t99 * t121;
t80 = t104 * qJDD(2) - t89;
t43 = (-t79 - t116) * pkin(9) + (-t80 + t89) * pkin(3) + t50;
t98 = sin(qJ(4));
t117 = t103 * t43 - t98 * t38;
t101 = cos(qJ(6));
t123 = qJD(2) * t99;
t75 = t103 * qJD(3) - t98 * t123;
t57 = t75 * qJD(4) + t98 * qJDD(3) + t103 * t79;
t72 = qJDD(4) - t80;
t76 = t98 * qJD(3) + t103 * t123;
t88 = qJD(4) - t122;
t26 = (t75 * t88 - t57) * pkin(10) + (t75 * t76 + t72) * pkin(4) + t117;
t126 = t103 * t38 + t98 * t43;
t56 = -t76 * qJD(4) + t103 * qJDD(3) - t98 * t79;
t65 = t88 * pkin(4) - t76 * pkin(10);
t71 = t75 ^ 2;
t28 = -t71 * pkin(4) + t56 * pkin(10) - t88 * t65 + t126;
t97 = sin(qJ(5));
t118 = t102 * t26 - t97 * t28;
t59 = t102 * t75 - t97 * t76;
t35 = t59 * qJD(5) + t102 * t57 + t97 * t56;
t60 = t102 * t76 + t97 * t75;
t70 = qJDD(5) + t72;
t87 = qJD(5) + t88;
t17 = (t59 * t87 - t35) * pkin(11) + (t59 * t60 + t70) * pkin(5) + t118;
t127 = t102 * t28 + t97 * t26;
t34 = -t60 * qJD(5) + t102 * t56 - t97 * t57;
t54 = t87 * pkin(5) - t60 * pkin(11);
t58 = t59 ^ 2;
t18 = -t58 * pkin(5) + t34 * pkin(11) - t87 * t54 + t127;
t96 = sin(qJ(6));
t45 = t101 * t59 - t96 * t60;
t23 = t45 * qJD(6) + t101 * t35 + t96 * t34;
t46 = t101 * t60 + t96 * t59;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t83 = qJD(6) + t87;
t41 = -t83 * mrSges(7,2) + t45 * mrSges(7,3);
t68 = qJDD(6) + t70;
t15 = m(7) * (t101 * t17 - t96 * t18) - t23 * mrSges(7,3) + t68 * mrSges(7,1) - t46 * t32 + t83 * t41;
t22 = -t46 * qJD(6) + t101 * t34 - t96 * t35;
t42 = t83 * mrSges(7,1) - t46 * mrSges(7,3);
t16 = m(7) * (t101 * t18 + t96 * t17) + t22 * mrSges(7,3) - t68 * mrSges(7,2) + t45 * t32 - t83 * t42;
t47 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t52 = -t87 * mrSges(6,2) + t59 * mrSges(6,3);
t12 = m(6) * t118 + t70 * mrSges(6,1) - t35 * mrSges(6,3) + t101 * t15 + t96 * t16 - t60 * t47 + t87 * t52;
t53 = t87 * mrSges(6,1) - t60 * mrSges(6,3);
t13 = m(6) * t127 - t70 * mrSges(6,2) + t34 * mrSges(6,3) + t101 * t16 - t96 * t15 + t59 * t47 - t87 * t53;
t61 = -t75 * mrSges(5,1) + t76 * mrSges(5,2);
t63 = -t88 * mrSges(5,2) + t75 * mrSges(5,3);
t10 = m(5) * t117 + t72 * mrSges(5,1) - t57 * mrSges(5,3) + t102 * t12 + t97 * t13 - t76 * t61 + t88 * t63;
t64 = t88 * mrSges(5,1) - t76 * mrSges(5,3);
t11 = m(5) * t126 - t72 * mrSges(5,2) + t56 * mrSges(5,3) + t102 * t13 - t97 * t12 + t75 * t61 - t88 * t64;
t84 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t123;
t85 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t122;
t130 = m(4) * t50 - t80 * mrSges(4,1) + t79 * mrSges(4,2) - (t104 * t85 - t99 * t84) * qJD(2) + t103 * t10 + t98 * t11;
t8 = m(3) * t131 + qJDD(2) * mrSges(3,1) - t107 * mrSges(3,2) - t130;
t129 = t105 * t8;
t115 = t104 * t66 - t99 * t51;
t37 = -qJDD(3) * pkin(3) - t106 * pkin(9) + t78 * t123 - t115;
t110 = -t56 * pkin(4) - t71 * pkin(10) + t76 * t65 + t37;
t112 = t22 * mrSges(7,1) + t45 * t41 - m(7) * (-t34 * pkin(5) - t58 * pkin(11) + t60 * t54 + t110) - t23 * mrSges(7,2) - t46 * t42;
t109 = -m(6) * t110 + t34 * mrSges(6,1) - t35 * mrSges(6,2) + t59 * t52 - t60 * t53 + t112;
t108 = m(5) * t37 - t56 * mrSges(5,1) + t57 * mrSges(5,2) - t75 * t63 + t76 * t64 - t109;
t77 = (-mrSges(4,1) * t104 + mrSges(4,2) * t99) * qJD(2);
t14 = m(4) * t115 + qJDD(3) * mrSges(4,1) - t79 * mrSges(4,3) + qJD(3) * t85 - t77 * t123 - t108;
t9 = m(4) * t125 - qJDD(3) * mrSges(4,2) + t80 * mrSges(4,3) - qJD(3) * t84 - t98 * t10 + t103 * t11 + t77 * t122;
t4 = m(3) * t120 - t107 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t104 * t9 - t99 * t14;
t6 = m(3) * t66 + t104 * t14 + t99 * t9;
t119 = m(2) * t91 + t4 * t124 + t93 * t129 + t95 * t6;
t2 = m(2) * t82 - t100 * t8 + t105 * t4;
t1 = m(2) * t81 - t93 * t6 + (t100 * t4 + t129) * t95;
t3 = [-m(1) * g(1) - t92 * t1 + t94 * t2, t2, t4, t9, t11, t13, t16; -m(1) * g(2) + t94 * t1 + t92 * t2, t1, t8, t14, t10, t12, t15; -m(1) * g(3) + t119, t119, t6, t130, t108, -t109, -t112;];
f_new  = t3;
