% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 02:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRP12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP12_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP12_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP12_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:15:07
% EndTime: 2019-05-06 02:15:27
% DurationCPUTime: 9.67s
% Computational Cost: add. (150726->213), mult. (467845->297), div. (0->0), fcn. (397274->14), ass. (0->119)
t105 = sin(qJ(3));
t108 = cos(qJ(3));
t100 = cos(pkin(12));
t102 = cos(pkin(6));
t101 = cos(pkin(7));
t99 = sin(pkin(6));
t143 = t101 * t99;
t98 = sin(pkin(7));
t118 = t100 * t143 + t102 * t98;
t97 = sin(pkin(12));
t153 = t97 * t99;
t114 = -t105 * t153 + t108 * t118;
t78 = t114 * qJD(1);
t112 = t105 * t118 + t108 * t153;
t79 = t112 * qJD(1);
t67 = -t79 * qJD(3) + qJDD(1) * t114;
t156 = pkin(9) * t97;
t122 = -pkin(2) * t100 - t156 * t98;
t140 = qJD(1) * t99;
t141 = pkin(9) * qJDD(1);
t116 = qJD(1) * t122 * t140 + t101 * t141;
t132 = qJD(2) * t140;
t110 = qJD(1) ^ 2;
t106 = sin(qJ(1));
t109 = cos(qJ(1));
t131 = t106 * g(1) - t109 * g(2);
t146 = qJ(2) * t99;
t89 = qJDD(1) * pkin(1) + t110 * t146 + t131;
t142 = t102 * t89;
t145 = t100 * t99;
t123 = -g(3) * t145 + t100 * t142 - 0.2e1 * t97 * t132;
t86 = t118 * qJD(1) * pkin(9);
t127 = -t109 * g(1) - t106 * g(2);
t90 = -t110 * pkin(1) + qJDD(1) * t146 + t127;
t51 = (pkin(2) * qJDD(1) + qJD(1) * t86) * t102 + (-t116 * t99 - t90) * t97 + t123;
t128 = -t102 * g(3) + qJDD(2);
t91 = (pkin(2) * t102 - t143 * t156) * qJD(1);
t60 = (-t89 + t122 * qJDD(1) + (-t100 * t86 + t91 * t97) * qJD(1)) * t99 + t128;
t160 = t101 * t51 + t60 * t98;
t134 = t97 * t142 + (0.2e1 * t132 + t90) * t100;
t52 = (-qJD(1) * t91 + t141 * t98) * t102 + (-g(3) * t97 + t100 * t116) * t99 + t134;
t159 = -t105 * t52 + t108 * t160;
t104 = sin(qJ(4));
t107 = cos(qJ(4));
t135 = t105 * t160 + t108 * t52;
t66 = -t78 * pkin(3) - t79 * pkin(10);
t119 = t101 * t102 - t145 * t98;
t87 = qJD(1) * t119 + qJD(3);
t83 = t87 ^ 2;
t84 = qJDD(1) * t119 + qJDD(3);
t28 = -t83 * pkin(3) + t84 * pkin(10) + t78 * t66 + t135;
t130 = t101 * t60 - t98 * t51;
t68 = t78 * qJD(3) + qJDD(1) * t112;
t30 = (-t78 * t87 - t68) * pkin(10) + (t79 * t87 - t67) * pkin(3) + t130;
t129 = -t104 * t28 + t107 * t30;
t71 = -t104 * t79 + t107 * t87;
t72 = t104 * t87 + t107 * t79;
t54 = -t71 * pkin(4) - t72 * pkin(11);
t64 = qJDD(4) - t67;
t77 = qJD(4) - t78;
t76 = t77 ^ 2;
t21 = -t64 * pkin(4) - t76 * pkin(11) + t72 * t54 - t129;
t103 = sin(qJ(5));
t155 = cos(qJ(5));
t47 = t71 * qJD(4) + t104 * t84 + t107 * t68;
t59 = t103 * t77 + t155 * t72;
t32 = t59 * qJD(5) + t103 * t47 - t155 * t64;
t58 = t103 * t72 - t155 * t77;
t33 = -t58 * qJD(5) + t103 * t64 + t155 * t47;
t70 = qJD(5) - t71;
t40 = -t58 * mrSges(7,2) + t70 * mrSges(7,3);
t136 = m(7) * (-0.2e1 * qJD(6) * t59 + (t58 * t70 - t33) * qJ(6) + (t59 * t70 + t32) * pkin(5) + t21) + t32 * mrSges(7,1) + t58 * t40;
t41 = -t70 * mrSges(6,2) - t58 * mrSges(6,3);
t42 = t70 * mrSges(6,1) - t59 * mrSges(6,3);
t43 = -t70 * mrSges(7,1) + t59 * mrSges(7,2);
t158 = m(6) * t21 + t32 * mrSges(6,1) + (t42 - t43) * t59 + (mrSges(6,2) - mrSges(7,3)) * t33 + t58 * t41 + t136;
t149 = t104 * t30 + t107 * t28;
t22 = -t76 * pkin(4) + t64 * pkin(11) + t71 * t54 + t149;
t27 = -t84 * pkin(3) - t83 * pkin(10) + t79 * t66 - t159;
t46 = -t72 * qJD(4) - t104 * t68 + t107 * t84;
t24 = (-t71 * t77 - t47) * pkin(11) + (t72 * t77 - t46) * pkin(4) + t27;
t117 = -t103 * t22 + t155 * t24;
t36 = t58 * pkin(5) - t59 * qJ(6);
t45 = qJDD(5) - t46;
t69 = t70 ^ 2;
t157 = m(7) * (-t45 * pkin(5) - t69 * qJ(6) + t59 * t36 + qJDD(6) - t117);
t151 = -mrSges(6,3) - mrSges(7,2);
t150 = t103 * t24 + t155 * t22;
t37 = t58 * mrSges(7,1) - t59 * mrSges(7,3);
t148 = -t58 * mrSges(6,1) - t59 * mrSges(6,2) - t37;
t137 = m(7) * (-t69 * pkin(5) + t45 * qJ(6) + 0.2e1 * qJD(6) * t70 - t58 * t36 + t150) + t70 * t43 + t45 * mrSges(7,3);
t14 = m(6) * t150 - t45 * mrSges(6,2) + t148 * t58 + t151 * t32 - t70 * t42 + t137;
t15 = m(6) * t117 - t157 + (t41 + t40) * t70 + t148 * t59 + (mrSges(6,1) + mrSges(7,1)) * t45 + t151 * t33;
t53 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t62 = t77 * mrSges(5,1) - t72 * mrSges(5,3);
t12 = m(5) * t149 - t64 * mrSges(5,2) + t46 * mrSges(5,3) - t103 * t15 + t14 * t155 + t71 * t53 - t77 * t62;
t61 = -t77 * mrSges(5,2) + t71 * mrSges(5,3);
t13 = m(5) * t129 + t64 * mrSges(5,1) - t47 * mrSges(5,3) - t72 * t53 + t77 * t61 - t158;
t73 = -t87 * mrSges(4,2) + t78 * mrSges(4,3);
t74 = t87 * mrSges(4,1) - t79 * mrSges(4,3);
t10 = m(4) * t130 - t67 * mrSges(4,1) + t68 * mrSges(4,2) + t104 * t12 + t107 * t13 - t78 * t73 + t79 * t74;
t121 = t102 * mrSges(3,1) - mrSges(3,3) * t153;
t111 = m(5) * t27 - t46 * mrSges(5,1) + t47 * mrSges(5,2) + t103 * t14 + t15 * t155 - t71 * t61 + t72 * t62;
t65 = -t78 * mrSges(4,1) + t79 * mrSges(4,2);
t11 = m(4) * t159 + t84 * mrSges(4,1) - t68 * mrSges(4,3) - t79 * t65 + t87 * t73 - t111;
t9 = m(4) * t135 - t84 * mrSges(4,2) + t67 * mrSges(4,3) - t104 * t13 + t107 * t12 + t78 * t65 - t87 * t74;
t124 = t105 * t9 + t108 * t11;
t126 = -mrSges(3,1) * t100 + mrSges(3,2) * t97;
t88 = t126 * t140;
t120 = -t102 * mrSges(3,2) + mrSges(3,3) * t145;
t93 = t120 * qJD(1);
t4 = m(3) * (-t97 * t90 + t123) - t98 * t10 + t124 * t101 + t121 * qJDD(1) + (t102 * t93 - t153 * t88) * qJD(1);
t92 = t121 * qJD(1);
t6 = m(3) * t128 + t101 * t10 + t124 * t98 + (-m(3) * t89 + t126 * qJDD(1) + (-t100 * t93 + t92 * t97) * qJD(1)) * t99;
t8 = m(3) * (-g(3) * t153 + t134) + t108 * t9 - t105 * t11 + t120 * qJDD(1) + (-t102 * t92 + t145 * t88) * qJD(1);
t138 = t102 * t6 + t4 * t145 + t8 * t153;
t2 = m(2) * t127 - t110 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t100 * t8 - t97 * t4;
t1 = m(2) * t131 + qJDD(1) * mrSges(2,1) - t110 * mrSges(2,2) - t99 * t6 + (t100 * t4 + t97 * t8) * t102;
t3 = [-m(1) * g(1) - t106 * t1 + t109 * t2, t2, t8, t9, t12, t14, -t32 * mrSges(7,2) - t58 * t37 + t137; -m(1) * g(2) + t109 * t1 + t106 * t2, t1, t4, t11, t13, t15, -t33 * mrSges(7,3) - t59 * t43 + t136; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t6, t10, t111, t158, -t45 * mrSges(7,1) + t33 * mrSges(7,2) + t59 * t37 - t70 * t40 + t157;];
f_new  = t3;
