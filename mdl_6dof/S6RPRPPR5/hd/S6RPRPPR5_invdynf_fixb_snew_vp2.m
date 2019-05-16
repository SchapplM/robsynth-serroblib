% Calculate vector of cutting forces with Newton-Euler
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-05-05 17:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:57:43
% EndTime: 2019-05-05 16:57:48
% DurationCPUTime: 2.47s
% Computational Cost: add. (27602->190), mult. (68561->239), div. (0->0), fcn. (49414->10), ass. (0->98)
t149 = -2 * qJD(4);
t105 = qJD(1) ^ 2;
t101 = sin(qJ(1));
t103 = cos(qJ(1));
t125 = t101 * g(1) - t103 * g(2);
t121 = qJDD(2) - t125;
t98 = cos(pkin(9));
t93 = t98 ^ 2;
t96 = sin(pkin(9));
t135 = t96 ^ 2 + t93;
t108 = (-pkin(2) * t98 - pkin(1)) * qJDD(1) + (-t135 * pkin(7) - qJ(2)) * t105 + t121;
t102 = cos(qJ(6));
t100 = sin(qJ(3));
t142 = cos(qJ(3));
t115 = t100 * t98 + t142 * t96;
t83 = t115 * qJD(1);
t131 = t83 * qJD(3);
t145 = t100 * t96 - t98 * t142;
t82 = t145 * qJD(1);
t132 = t82 * qJD(3);
t64 = t115 * qJDD(1) - t132;
t106 = pkin(3) * t131 + t83 * t149 + (-t64 + t132) * qJ(4) + t108;
t63 = t145 * qJDD(1) + t131;
t73 = t83 * pkin(4) - qJD(3) * qJ(5);
t81 = t82 ^ 2;
t19 = -t81 * pkin(4) - t83 * t73 + (pkin(3) + qJ(5)) * t63 + t106;
t104 = qJD(3) ^ 2;
t130 = qJD(1) * qJD(2);
t126 = -t98 * g(3) - 0.2e1 * t96 * t130;
t133 = pkin(7) * qJDD(1);
t141 = pkin(2) * t105;
t119 = -t103 * g(1) - t101 * g(2);
t84 = -t105 * pkin(1) + qJDD(1) * qJ(2) + t119;
t45 = (t98 * t141 - t133 - t84) * t96 + t126;
t124 = -t96 * g(3) + (0.2e1 * t130 + t84) * t98;
t52 = t98 * t133 - t93 * t141 + t124;
t122 = -t100 * t52 + t142 * t45;
t56 = t82 * pkin(3) - t83 * qJ(4);
t29 = -qJDD(3) * pkin(3) - t104 * qJ(4) + t83 * t56 + qJDD(4) - t122;
t22 = (t82 * t83 - qJDD(3)) * qJ(5) + (t64 + t132) * pkin(4) + t29;
t95 = sin(pkin(10));
t97 = cos(pkin(10));
t69 = t97 * qJD(3) + t95 * t82;
t123 = -0.2e1 * qJD(5) * t69 - t95 * t19 + t97 * t22;
t51 = t97 * qJDD(3) + t95 * t63;
t68 = -t95 * qJD(3) + t97 * t82;
t14 = (t68 * t83 - t51) * pkin(8) + (t68 * t69 + t64) * pkin(5) + t123;
t128 = 0.2e1 * qJD(5) * t68 + t97 * t19 + t95 * t22;
t49 = t83 * pkin(5) - t69 * pkin(8);
t50 = -t95 * qJDD(3) + t97 * t63;
t67 = t68 ^ 2;
t15 = -t67 * pkin(5) + t50 * pkin(8) - t83 * t49 + t128;
t99 = sin(qJ(6));
t37 = t102 * t68 - t99 * t69;
t31 = t37 * qJD(6) + t102 * t51 + t99 * t50;
t38 = t102 * t69 + t99 * t68;
t33 = -t37 * mrSges(7,1) + t38 * mrSges(7,2);
t79 = qJD(6) + t83;
t34 = -t79 * mrSges(7,2) + t37 * mrSges(7,3);
t61 = qJDD(6) + t64;
t12 = m(7) * (t102 * t14 - t99 * t15) - t31 * mrSges(7,3) + t61 * mrSges(7,1) - t38 * t33 + t79 * t34;
t30 = -t38 * qJD(6) + t102 * t50 - t99 * t51;
t35 = t79 * mrSges(7,1) - t38 * mrSges(7,3);
t13 = m(7) * (t102 * t15 + t99 * t14) + t30 * mrSges(7,3) - t61 * mrSges(7,2) + t37 * t33 - t79 * t35;
t39 = -t68 * mrSges(6,1) + t69 * mrSges(6,2);
t48 = t83 * mrSges(6,1) - t69 * mrSges(6,3);
t10 = m(6) * t128 - t64 * mrSges(6,2) + t50 * mrSges(6,3) + t102 * t13 - t99 * t12 + t68 * t39 - t83 * t48;
t75 = t83 * mrSges(5,1) + qJD(3) * mrSges(5,2);
t47 = -t83 * mrSges(6,2) + t68 * mrSges(6,3);
t9 = m(6) * t123 + t64 * mrSges(6,1) - t51 * mrSges(6,3) + t102 * t12 + t99 * t13 - t69 * t39 + t83 * t47;
t117 = -t97 * t10 + t95 * t9 - m(5) * (t63 * pkin(3) + t106) + t83 * t75 + t64 * mrSges(5,3);
t74 = t82 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t136 = -qJD(3) * mrSges(4,2) - t82 * mrSges(4,3) - t74;
t140 = mrSges(4,1) - mrSges(5,2);
t72 = qJD(3) * mrSges(4,1) - t83 * mrSges(4,3);
t109 = m(4) * t108 + t64 * mrSges(4,2) + t136 * t82 + t140 * t63 + t83 * t72 - t117;
t148 = -t109 - m(3) * (-qJDD(1) * pkin(1) - t105 * qJ(2) + t121);
t147 = t135 * mrSges(3,3);
t138 = t100 * t45 + t142 * t52;
t146 = t104 * pkin(3) - qJDD(3) * qJ(4) + qJD(3) * t149 + t82 * t56 - t138;
t111 = -t63 * pkin(4) - t81 * qJ(5) + qJD(3) * t73 + qJDD(5) - t146;
t113 = -t30 * mrSges(7,1) - t37 * t34 + m(7) * (-t50 * pkin(5) - t67 * pkin(8) + t69 * t49 + t111) + t31 * mrSges(7,2) + t38 * t35;
t110 = m(6) * t111 - t50 * mrSges(6,1) + t51 * mrSges(6,2) - t68 * t47 + t69 * t48 + t113;
t107 = -m(5) * t146 + t110;
t58 = -t82 * mrSges(5,2) - t83 * mrSges(5,3);
t137 = -t82 * mrSges(4,1) - t83 * mrSges(4,2) - t58;
t139 = -mrSges(4,3) - mrSges(5,1);
t11 = m(4) * t138 + t107 + (-t72 + t75) * qJD(3) + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) + t139 * t63 + t137 * t82;
t120 = -t98 * mrSges(3,1) + t96 * mrSges(3,2);
t118 = qJDD(1) * mrSges(3,3) + t105 * t120;
t114 = -m(5) * t29 - t95 * t10 - t97 * t9;
t7 = m(4) * t122 + t136 * qJD(3) + t140 * qJDD(3) + t137 * t83 + t139 * t64 + t114;
t4 = m(3) * t126 + t100 * t11 + t142 * t7 + (-m(3) * t84 - t118) * t96;
t5 = m(3) * t124 - t100 * t7 + t142 * t11 + t118 * t98;
t144 = t98 * t4 + t96 * t5;
t6 = m(2) * t125 + (-mrSges(2,2) + t147) * t105 + (mrSges(2,1) - t120) * qJDD(1) + t148;
t1 = m(2) * t119 - t105 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t96 * t4 + t98 * t5;
t2 = [-m(1) * g(1) + t103 * t1 - t101 * t6, t1, t5, t11, -t63 * mrSges(5,2) - t82 * t74 - t117, t10, t13; -m(1) * g(2) + t101 * t1 + t103 * t6, t6, t4, t7, t63 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t75 + t82 * t58 - t107, t9, t12; (-m(1) - m(2)) * g(3) + t144, -m(2) * g(3) + t144, t120 * qJDD(1) - t105 * t147 - t148, t109, t64 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t74 + t83 * t58 - t114, t110, t113;];
f_new  = t2;
