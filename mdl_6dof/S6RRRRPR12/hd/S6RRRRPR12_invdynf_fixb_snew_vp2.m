% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR12
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-08 00:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR12_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR12_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR12_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 00:09:51
% EndTime: 2019-05-08 00:10:52
% DurationCPUTime: 21.84s
% Computational Cost: add. (409853->229), mult. (1016046->326), div. (0->0), fcn. (862321->16), ass. (0->130)
t105 = sin(pkin(7));
t112 = sin(qJ(3));
t117 = cos(qJ(3));
t108 = cos(pkin(7));
t109 = cos(pkin(6));
t101 = t109 * qJDD(1) + qJDD(2);
t102 = t109 * qJD(1) + qJD(2);
t106 = sin(pkin(6));
t118 = cos(qJ(2));
t113 = sin(qJ(2));
t120 = qJD(1) ^ 2;
t114 = sin(qJ(1));
t119 = cos(qJ(1));
t135 = t114 * g(1) - t119 * g(2);
t158 = pkin(9) * t106;
t96 = qJDD(1) * pkin(1) + t120 * t158 + t135;
t151 = t109 * t96;
t131 = -t119 * g(1) - t114 * g(2);
t97 = -t120 * pkin(1) + qJDD(1) * t158 + t131;
t132 = -t113 * t97 + t118 * t151;
t144 = qJD(1) * t113;
t156 = pkin(10) * t108;
t143 = qJD(1) * t118;
t136 = t106 * t143;
t150 = t102 * t105;
t89 = (t108 * t136 + t150) * pkin(10);
t145 = qJD(1) * t106;
t157 = pkin(10) * t105;
t93 = (-pkin(2) * t118 - t113 * t157) * t145;
t141 = qJD(1) * qJD(2);
t99 = (qJDD(1) * t113 + t118 * t141) * t106;
t59 = -t99 * t156 + t101 * pkin(2) + t102 * t89 + (-g(3) * t118 - t93 * t144) * t106 + t132;
t152 = t108 * t59;
t100 = (qJDD(1) * t118 - t113 * t141) * t106;
t126 = t100 * t108 + t101 * t105;
t153 = t113 * t151 + t118 * t97;
t137 = t106 * t144;
t92 = t102 * pkin(2) - t137 * t156;
t60 = -t102 * t92 + (-g(3) * t113 + t93 * t143) * t106 + t126 * pkin(10) + t153;
t155 = t109 * g(3);
t66 = -t99 * t157 - t100 * pkin(2) - t155 + (-t96 + (t113 * t92 - t118 * t89) * qJD(1)) * t106;
t160 = -t112 * t60 + (t105 * t66 + t152) * t117;
t146 = t108 * t118;
t149 = t105 * t112;
t84 = t102 * t149 + (t112 * t146 + t113 * t117) * t145;
t72 = -t84 * qJD(3) - t112 * t99 + t126 * t117;
t83 = (-t112 * t113 + t117 * t146) * t145 + t117 * t150;
t159 = 2 * qJD(5);
t111 = sin(qJ(4));
t116 = cos(qJ(4));
t139 = t112 * t152 + t117 * t60 + t66 * t149;
t75 = -t83 * pkin(3) - t84 * pkin(11);
t85 = -t105 * t100 + t108 * t101 + qJDD(3);
t90 = t108 * t102 - t105 * t136 + qJD(3);
t88 = t90 ^ 2;
t33 = -t88 * pkin(3) + t85 * pkin(11) + t83 * t75 + t139;
t134 = -t105 * t59 + t108 * t66;
t73 = t83 * qJD(3) + t126 * t112 + t117 * t99;
t36 = (-t83 * t90 - t73) * pkin(11) + (t84 * t90 - t72) * pkin(3) + t134;
t154 = t111 * t36 + t116 * t33;
t148 = t106 * t113;
t147 = t106 * t118;
t104 = sin(pkin(13));
t107 = cos(pkin(13));
t110 = sin(qJ(6));
t115 = cos(qJ(6));
t133 = -t111 * t33 + t116 * t36;
t77 = -t111 * t84 + t116 * t90;
t48 = t77 * qJD(4) + t111 * t85 + t116 * t73;
t71 = qJDD(4) - t72;
t78 = t111 * t90 + t116 * t84;
t82 = qJD(4) - t83;
t24 = (t77 * t82 - t48) * qJ(5) + (t77 * t78 + t71) * pkin(4) + t133;
t47 = -t78 * qJD(4) - t111 * t73 + t116 * t85;
t69 = t82 * pkin(4) - t78 * qJ(5);
t76 = t77 ^ 2;
t26 = -t76 * pkin(4) + t47 * qJ(5) - t82 * t69 + t154;
t61 = -t104 * t78 + t107 * t77;
t140 = t104 * t24 + t107 * t26 + t61 * t159;
t62 = t104 * t77 + t107 * t78;
t46 = -t61 * pkin(5) - t62 * pkin(12);
t81 = t82 ^ 2;
t21 = -t81 * pkin(5) + t71 * pkin(12) + t61 * t46 + t140;
t32 = -t85 * pkin(3) - t88 * pkin(11) + t84 * t75 - t160;
t122 = -t47 * pkin(4) - t76 * qJ(5) + t78 * t69 + qJDD(5) + t32;
t39 = -t104 * t48 + t107 * t47;
t40 = t104 * t47 + t107 * t48;
t22 = (-t61 * t82 - t40) * pkin(12) + (t62 * t82 - t39) * pkin(5) + t122;
t49 = -t110 * t62 + t115 * t82;
t30 = t49 * qJD(6) + t110 * t71 + t115 * t40;
t38 = qJDD(6) - t39;
t50 = t110 * t82 + t115 * t62;
t41 = -t49 * mrSges(7,1) + t50 * mrSges(7,2);
t58 = qJD(6) - t61;
t42 = -t58 * mrSges(7,2) + t49 * mrSges(7,3);
t18 = m(7) * (-t110 * t21 + t115 * t22) - t30 * mrSges(7,3) + t38 * mrSges(7,1) - t50 * t41 + t58 * t42;
t29 = -t50 * qJD(6) - t110 * t40 + t115 * t71;
t43 = t58 * mrSges(7,1) - t50 * mrSges(7,3);
t19 = m(7) * (t110 * t22 + t115 * t21) + t29 * mrSges(7,3) - t38 * mrSges(7,2) + t49 * t41 - t58 * t43;
t45 = -t61 * mrSges(6,1) + t62 * mrSges(6,2);
t52 = t82 * mrSges(6,1) - t62 * mrSges(6,3);
t13 = m(6) * t140 - t71 * mrSges(6,2) + t39 * mrSges(6,3) - t110 * t18 + t115 * t19 + t61 * t45 - t82 * t52;
t129 = t104 * t26 - t107 * t24;
t123 = m(7) * (-t71 * pkin(5) - t81 * pkin(12) + (t159 + t46) * t62 + t129) - t29 * mrSges(7,1) + t30 * mrSges(7,2) - t49 * t42 + t50 * t43;
t51 = -t82 * mrSges(6,2) + t61 * mrSges(6,3);
t15 = m(6) * (-0.2e1 * qJD(5) * t62 - t129) - t40 * mrSges(6,3) + t71 * mrSges(6,1) - t62 * t45 + t82 * t51 - t123;
t63 = -t77 * mrSges(5,1) + t78 * mrSges(5,2);
t68 = -t82 * mrSges(5,2) + t77 * mrSges(5,3);
t11 = m(5) * t133 + t71 * mrSges(5,1) - t48 * mrSges(5,3) + t104 * t13 + t107 * t15 - t78 * t63 + t82 * t68;
t70 = t82 * mrSges(5,1) - t78 * mrSges(5,3);
t12 = m(5) * t154 - t71 * mrSges(5,2) + t47 * mrSges(5,3) - t104 * t15 + t107 * t13 + t77 * t63 - t82 * t70;
t79 = -t90 * mrSges(4,2) + t83 * mrSges(4,3);
t80 = t90 * mrSges(4,1) - t84 * mrSges(4,3);
t10 = m(4) * t134 - t72 * mrSges(4,1) + t73 * mrSges(4,2) + t116 * t11 + t111 * t12 - t83 * t79 + t84 * t80;
t124 = -m(6) * t122 + t39 * mrSges(6,1) - t40 * mrSges(6,2) - t110 * t19 - t115 * t18 + t61 * t51 - t62 * t52;
t121 = m(5) * t32 - t47 * mrSges(5,1) + t48 * mrSges(5,2) - t77 * t68 + t78 * t70 - t124;
t74 = -t83 * mrSges(4,1) + t84 * mrSges(4,2);
t14 = m(4) * t160 + t85 * mrSges(4,1) - t73 * mrSges(4,3) - t84 * t74 + t90 * t79 - t121;
t9 = m(4) * t139 - t85 * mrSges(4,2) + t72 * mrSges(4,3) - t111 * t11 + t116 * t12 + t83 * t74 - t90 * t80;
t130 = t112 * t9 + t117 * t14;
t95 = -t102 * mrSges(3,2) + mrSges(3,3) * t136;
t98 = (-mrSges(3,1) * t118 + mrSges(3,2) * t113) * t145;
t4 = m(3) * (-g(3) * t147 + t132) - t99 * mrSges(3,3) + t101 * mrSges(3,1) - t98 * t137 + t102 * t95 - t105 * t10 + t130 * t108;
t94 = t102 * mrSges(3,1) - mrSges(3,3) * t137;
t6 = m(3) * (-t106 * t96 - t155) + t99 * mrSges(3,2) - t100 * mrSges(3,1) + t108 * t10 + t130 * t105 + (t113 * t94 - t118 * t95) * t145;
t8 = m(3) * (-g(3) * t148 + t153) + t100 * mrSges(3,3) - t101 * mrSges(3,2) + t98 * t136 - t102 * t94 + t117 * t9 - t112 * t14;
t142 = t109 * t6 + t4 * t147 + t8 * t148;
t2 = m(2) * t131 - t120 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t113 * t4 + t118 * t8;
t1 = m(2) * t135 + qJDD(1) * mrSges(2,1) - t120 * mrSges(2,2) - t106 * t6 + (t113 * t8 + t118 * t4) * t109;
t3 = [-m(1) * g(1) - t114 * t1 + t119 * t2, t2, t8, t9, t12, t13, t19; -m(1) * g(2) + t119 * t1 + t114 * t2, t1, t4, t14, t11, t15, t18; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t6, t10, t121, -t124, t123;];
f_new  = t3;
