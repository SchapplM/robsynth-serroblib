% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 05:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:56:41
% EndTime: 2019-05-08 04:56:58
% DurationCPUTime: 4.94s
% Computational Cost: add. (67288->201), mult. (137610->254), div. (0->0), fcn. (99728->10), ass. (0->99)
t109 = qJD(2) ^ 2;
t102 = sin(qJ(2));
t132 = qJD(1) * t102;
t107 = cos(qJ(2));
t110 = qJD(1) ^ 2;
t103 = sin(qJ(1));
t108 = cos(qJ(1));
t118 = -t108 * g(1) - t103 * g(2);
t81 = -t110 * pkin(1) + qJDD(1) * pkin(7) + t118;
t133 = -t107 * g(3) - t102 * t81;
t88 = (-pkin(2) * t107 - pkin(8) * t102) * qJD(1);
t60 = -qJDD(2) * pkin(2) - t109 * pkin(8) + t88 * t132 - t133;
t101 = sin(qJ(3));
t106 = cos(qJ(3));
t86 = t101 * qJD(2) + t106 * t132;
t130 = qJD(1) * qJD(2);
t121 = t107 * t130;
t89 = t102 * qJDD(1) + t121;
t66 = -t86 * qJD(3) + t106 * qJDD(2) - t101 * t89;
t131 = t107 * qJD(1);
t95 = qJD(3) - t131;
t74 = t95 * pkin(3) - t86 * pkin(9);
t85 = t106 * qJD(2) - t101 * t132;
t83 = t85 ^ 2;
t113 = -t66 * pkin(3) - t83 * pkin(9) + t86 * t74 + t60;
t100 = sin(qJ(4));
t105 = cos(qJ(4));
t67 = t85 * qJD(3) + t101 * qJDD(2) + t106 * t89;
t70 = t100 * t85 + t105 * t86;
t42 = -t70 * qJD(4) - t100 * t67 + t105 * t66;
t94 = qJD(4) + t95;
t64 = t94 * pkin(4) - t70 * pkin(10);
t69 = -t100 * t86 + t105 * t85;
t68 = t69 ^ 2;
t112 = -t42 * pkin(4) - t68 * pkin(10) + t70 * t64 + t113;
t104 = cos(qJ(5));
t43 = t69 * qJD(4) + t100 * t66 + t105 * t67;
t99 = sin(qJ(5));
t54 = t104 * t70 + t99 * t69;
t27 = -t54 * qJD(5) + t104 * t42 - t99 * t43;
t53 = t104 * t69 - t99 * t70;
t28 = t53 * qJD(5) + t104 * t43 + t99 * t42;
t91 = qJD(5) + t94;
t47 = t91 * pkin(5) - t54 * qJ(6);
t48 = t91 * mrSges(7,1) - t54 * mrSges(7,3);
t52 = t53 ^ 2;
t127 = m(7) * (-t27 * pkin(5) - t52 * qJ(6) + t54 * t47 + qJDD(6) + t112) + t28 * mrSges(7,2) + t54 * t48;
t45 = -t91 * mrSges(7,2) + t53 * mrSges(7,3);
t46 = -t91 * mrSges(6,2) + t53 * mrSges(6,3);
t49 = t91 * mrSges(6,1) - t54 * mrSges(6,3);
t114 = m(6) * t112 + t28 * mrSges(6,2) + t54 * t49 + t127 - (t45 + t46) * t53 - (mrSges(6,1) + mrSges(7,1)) * t27;
t62 = -t94 * mrSges(5,2) + t69 * mrSges(5,3);
t63 = t94 * mrSges(5,1) - t70 * mrSges(5,3);
t146 = m(5) * t113 - t42 * mrSges(5,1) + t43 * mrSges(5,2) - t69 * t62 + t70 * t63 + t114;
t72 = -t95 * mrSges(4,2) + t85 * mrSges(4,3);
t73 = t95 * mrSges(4,1) - t86 * mrSges(4,3);
t145 = m(4) * t60 - t66 * mrSges(4,1) + t67 * mrSges(4,2) - t85 * t72 + t86 * t73 + t146;
t124 = t103 * g(1) - t108 * g(2);
t80 = -qJDD(1) * pkin(1) - t110 * pkin(7) - t124;
t96 = t102 * t130;
t90 = t107 * qJDD(1) - t96;
t58 = (-t89 - t121) * pkin(8) + (-t90 + t96) * pkin(2) + t80;
t123 = -t102 * g(3) + t107 * t81;
t61 = -t109 * pkin(2) + qJDD(2) * pkin(8) + t88 * t131 + t123;
t119 = -t101 * t61 + t106 * t58;
t84 = qJDD(3) - t90;
t32 = (t85 * t95 - t67) * pkin(9) + (t85 * t86 + t84) * pkin(3) + t119;
t134 = t101 * t58 + t106 * t61;
t34 = -t83 * pkin(3) + t66 * pkin(9) - t95 * t74 + t134;
t120 = -t100 * t34 + t105 * t32;
t82 = qJDD(4) + t84;
t19 = (t69 * t94 - t43) * pkin(10) + (t69 * t70 + t82) * pkin(4) + t120;
t136 = t100 * t32 + t105 * t34;
t21 = -t68 * pkin(4) + t42 * pkin(10) - t94 * t64 + t136;
t137 = t104 * t21 + t99 * t19;
t37 = -t53 * mrSges(7,1) + t54 * mrSges(7,2);
t128 = m(7) * (-t52 * pkin(5) + t27 * qJ(6) + 0.2e1 * qJD(6) * t53 - t91 * t47 + t137) + t27 * mrSges(7,3) + t53 * t37;
t38 = -t53 * mrSges(6,1) + t54 * mrSges(6,2);
t79 = qJDD(5) + t82;
t10 = m(6) * t137 + t27 * mrSges(6,3) + t53 * t38 + (-t49 - t48) * t91 + (-mrSges(6,2) - mrSges(7,2)) * t79 + t128;
t55 = -t69 * mrSges(5,1) + t70 * mrSges(5,2);
t122 = t104 * t19 - t99 * t21;
t129 = m(7) * (-0.2e1 * qJD(6) * t54 + (t53 * t91 - t28) * qJ(6) + (t53 * t54 + t79) * pkin(5) + t122) + t91 * t45 + t79 * mrSges(7,1);
t9 = m(6) * t122 + t79 * mrSges(6,1) + t91 * t46 + (-t38 - t37) * t54 + (-mrSges(6,3) - mrSges(7,3)) * t28 + t129;
t7 = m(5) * t120 + t82 * mrSges(5,1) - t43 * mrSges(5,3) + t99 * t10 + t104 * t9 - t70 * t55 + t94 * t62;
t71 = -t85 * mrSges(4,1) + t86 * mrSges(4,2);
t8 = m(5) * t136 - t82 * mrSges(5,2) + t42 * mrSges(5,3) + t104 * t10 + t69 * t55 - t94 * t63 - t99 * t9;
t5 = m(4) * t119 + t84 * mrSges(4,1) - t67 * mrSges(4,3) + t100 * t8 + t105 * t7 - t86 * t71 + t95 * t72;
t6 = m(4) * t134 - t84 * mrSges(4,2) + t66 * mrSges(4,3) - t100 * t7 + t105 * t8 + t85 * t71 - t95 * t73;
t92 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t132;
t93 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t131;
t143 = m(3) * t80 - t90 * mrSges(3,1) + t89 * mrSges(3,2) + t101 * t6 + t106 * t5 + (t102 * t92 - t107 * t93) * qJD(1);
t87 = (-mrSges(3,1) * t107 + mrSges(3,2) * t102) * qJD(1);
t12 = m(3) * t133 + qJDD(2) * mrSges(3,1) - t89 * mrSges(3,3) + qJD(2) * t93 - t87 * t132 - t145;
t4 = m(3) * t123 - qJDD(2) * mrSges(3,2) + t90 * mrSges(3,3) - qJD(2) * t92 - t101 * t5 + t106 * t6 + t87 * t131;
t141 = t102 * t4 + t107 * t12;
t2 = m(2) * t124 + qJDD(1) * mrSges(2,1) - t110 * mrSges(2,2) - t143;
t1 = m(2) * t118 - t110 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t102 * t12 + t107 * t4;
t3 = [-m(1) * g(1) + t108 * t1 - t103 * t2, t1, t4, t6, t8, t10, -t79 * mrSges(7,2) - t91 * t48 + t128; -m(1) * g(2) + t103 * t1 + t108 * t2, t2, t12, t5, t7, t9, -t28 * mrSges(7,3) - t54 * t37 + t129; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t143, t145, t146, t114, -t27 * mrSges(7,1) - t53 * t45 + t127;];
f_new  = t3;
