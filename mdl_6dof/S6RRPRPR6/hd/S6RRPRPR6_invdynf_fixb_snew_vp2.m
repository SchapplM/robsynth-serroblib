% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 14:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:20:25
% EndTime: 2019-05-06 14:20:37
% DurationCPUTime: 4.50s
% Computational Cost: add. (53851->211), mult. (141600->280), div. (0->0), fcn. (109997->12), ass. (0->110)
t104 = sin(pkin(11));
t106 = cos(pkin(11));
t107 = cos(pkin(6));
t100 = t107 * qJD(1) + qJD(2);
t105 = sin(pkin(6));
t110 = sin(qJ(2));
t113 = cos(qJ(2));
t115 = qJD(1) ^ 2;
t111 = sin(qJ(1));
t114 = cos(qJ(1));
t131 = t111 * g(1) - t114 * g(2);
t150 = pkin(8) * t105;
t90 = qJDD(1) * pkin(1) + t115 * t150 + t131;
t143 = t107 * t90;
t127 = -t114 * g(1) - t111 * g(2);
t91 = -t115 * pkin(1) + qJDD(1) * t150 + t127;
t129 = -t110 * t91 + t113 * t143;
t141 = t105 ^ 2 * t115;
t136 = qJD(1) * qJD(2);
t93 = (qJDD(1) * t110 + t113 * t136) * t105;
t99 = t107 * qJDD(1) + qJDD(2);
t41 = t99 * pkin(2) - t93 * qJ(3) + (pkin(2) * t110 * t141 + (qJ(3) * qJD(1) * t100 - g(3)) * t105) * t113 + t129;
t140 = t105 * t110;
t125 = -g(3) * t140 + t110 * t143 + t113 * t91;
t134 = t113 ^ 2 * t141;
t138 = qJD(1) * t105;
t133 = t110 * t138;
t87 = t100 * pkin(2) - qJ(3) * t133;
t94 = (qJDD(1) * t113 - t110 * t136) * t105;
t44 = -pkin(2) * t134 + t94 * qJ(3) - t100 * t87 + t125;
t85 = (t104 * t113 + t106 * t110) * t138;
t154 = -0.2e1 * qJD(3) * t85 - t104 * t44 + t106 * t41;
t108 = sin(qJ(6));
t112 = cos(qJ(6));
t109 = sin(qJ(4));
t151 = cos(qJ(4));
t73 = -t151 * t100 + t109 * t85;
t132 = t113 * t138;
t84 = -t104 * t133 + t106 * t132;
t83 = qJD(4) - t84;
t149 = t73 * t83;
t152 = -2 * qJD(5);
t65 = -t84 * pkin(3) - t85 * pkin(9);
t98 = t100 ^ 2;
t28 = -t99 * pkin(3) - t98 * pkin(9) + t85 * t65 - t154;
t69 = t104 * t94 + t106 * t93;
t49 = -t73 * qJD(4) + t109 * t99 + t151 * t69;
t74 = t109 * t100 + t151 * t85;
t116 = (-t49 + t149) * qJ(5) + t28 + (t83 * pkin(4) + t152) * t74;
t135 = 0.2e1 * qJD(3) * t84 + t104 * t41 + t106 * t44;
t29 = -t98 * pkin(3) + t99 * pkin(9) + t84 * t65 + t135;
t126 = -t107 * g(3) - t105 * t90;
t117 = -t94 * pkin(2) - qJ(3) * t134 + t87 * t133 + qJDD(3) + t126;
t68 = -t104 * t93 + t106 * t94;
t31 = (-t100 * t84 - t69) * pkin(9) + (t100 * t85 - t68) * pkin(3) + t117;
t128 = -t109 * t29 + t151 * t31;
t51 = t73 * pkin(4) - t74 * qJ(5);
t67 = qJDD(4) - t68;
t82 = t83 ^ 2;
t23 = -t67 * pkin(4) - t82 * qJ(5) + t74 * t51 + qJDD(5) - t128;
t18 = (t73 * t74 - t67) * pkin(10) + (t49 + t149) * pkin(5) + t23;
t48 = t74 * qJD(4) + t109 * t69 - t151 * t99;
t61 = t74 * pkin(5) - t83 * pkin(10);
t72 = t73 ^ 2;
t21 = -t72 * pkin(5) - t74 * t61 + (pkin(4) + pkin(10)) * t48 + t116;
t55 = -t108 * t83 + t112 * t73;
t34 = t55 * qJD(6) + t108 * t48 + t112 * t67;
t56 = t108 * t73 + t112 * t83;
t36 = -mrSges(7,1) * t55 + mrSges(7,2) * t56;
t71 = qJD(6) + t74;
t42 = -t71 * mrSges(7,2) + t55 * mrSges(7,3);
t47 = qJDD(6) + t49;
t16 = m(7) * (-t108 * t21 + t112 * t18) - t34 * mrSges(7,3) + t47 * mrSges(7,1) - t56 * t36 + t71 * t42;
t33 = -t56 * qJD(6) - t108 * t67 + t112 * t48;
t43 = t71 * mrSges(7,1) - t56 * mrSges(7,3);
t17 = m(7) * (t108 * t18 + t112 * t21) + t33 * mrSges(7,3) - t47 * mrSges(7,2) + t55 * t36 - t71 * t43;
t58 = t74 * mrSges(6,1) + t83 * mrSges(6,2);
t124 = t108 * t16 - t112 * t17 - m(6) * (t48 * pkin(4) + t116) + t74 * t58 + t49 * mrSges(6,3);
t57 = t73 * mrSges(6,1) - t83 * mrSges(6,3);
t144 = -t83 * mrSges(5,2) - t73 * mrSges(5,3) - t57;
t148 = mrSges(5,1) - mrSges(6,2);
t60 = t83 * mrSges(5,1) - t74 * mrSges(5,3);
t153 = m(5) * t28 + t49 * mrSges(5,2) + t144 * t73 + t148 * t48 + t74 * t60 - t124;
t147 = -mrSges(5,3) - mrSges(6,1);
t146 = t109 * t31 + t151 * t29;
t53 = -t73 * mrSges(6,2) - t74 * mrSges(6,3);
t145 = -t73 * mrSges(5,1) - t74 * mrSges(5,2) - t53;
t139 = t105 * t113;
t64 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t75 = -t100 * mrSges(4,2) + t84 * mrSges(4,3);
t10 = m(4) * t154 + t99 * mrSges(4,1) - t69 * mrSges(4,3) + t100 * t75 - t85 * t64 - t153;
t123 = -m(6) * t23 - t108 * t17 - t112 * t16;
t12 = m(5) * t128 + t144 * t83 + t145 * t74 + t147 * t49 + t148 * t67 + t123;
t121 = -t82 * pkin(4) + t67 * qJ(5) - t73 * t51 + t146;
t122 = -t33 * mrSges(7,1) - t55 * t42 + m(7) * (-t48 * pkin(5) - t72 * pkin(10) + ((2 * qJD(5)) + t61) * t83 + t121) + t34 * mrSges(7,2) + t56 * t43;
t119 = -m(6) * (t83 * t152 - t121) + t122;
t14 = m(5) * t146 + (-t60 + t58) * t83 + t145 * t73 + (-mrSges(5,2) + mrSges(6,3)) * t67 + t147 * t48 + t119;
t76 = t100 * mrSges(4,1) - t85 * mrSges(4,3);
t7 = m(4) * t135 - t99 * mrSges(4,2) + t68 * mrSges(4,3) - t100 * t76 - t109 * t12 + t151 * t14 + t84 * t64;
t89 = -t100 * mrSges(3,2) + mrSges(3,3) * t132;
t92 = (-mrSges(3,1) * t113 + mrSges(3,2) * t110) * t138;
t5 = m(3) * (-g(3) * t139 + t129) - t93 * mrSges(3,3) + t99 * mrSges(3,1) - t92 * t133 + t100 * t89 + t104 * t7 + t106 * t10;
t88 = t100 * mrSges(3,1) - mrSges(3,3) * t133;
t6 = m(3) * t125 - t99 * mrSges(3,2) + t94 * mrSges(3,3) - t104 * t10 - t100 * t88 + t106 * t7 + t92 * t132;
t120 = m(4) * t117 - t68 * mrSges(4,1) + t69 * mrSges(4,2) + t109 * t14 + t151 * t12 - t84 * t75 + t85 * t76;
t9 = m(3) * t126 + t93 * mrSges(3,2) - t94 * mrSges(3,1) + (t110 * t88 - t113 * t89) * t138 + t120;
t137 = t107 * t9 + t5 * t139 + t6 * t140;
t2 = m(2) * t127 - t115 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t110 * t5 + t113 * t6;
t1 = m(2) * t131 + qJDD(1) * mrSges(2,1) - t115 * mrSges(2,2) - t105 * t9 + (t110 * t6 + t113 * t5) * t107;
t3 = [-m(1) * g(1) - t111 * t1 + t114 * t2, t2, t6, t7, t14, -t48 * mrSges(6,2) - t73 * t57 - t124, t17; -m(1) * g(2) + t114 * t1 + t111 * t2, t1, t5, t10, t12, t48 * mrSges(6,1) - t67 * mrSges(6,3) + t73 * t53 - t83 * t58 - t119, t16; (-m(1) - m(2)) * g(3) + t137, -m(2) * g(3) + t137, t9, t120, t153, t49 * mrSges(6,1) + t67 * mrSges(6,2) + t74 * t53 + t83 * t57 - t123, t122;];
f_new  = t3;
