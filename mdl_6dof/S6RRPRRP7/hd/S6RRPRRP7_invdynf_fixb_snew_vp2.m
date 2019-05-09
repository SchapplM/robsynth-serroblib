% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 18:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:12:34
% EndTime: 2019-05-06 18:12:38
% DurationCPUTime: 1.58s
% Computational Cost: add. (16758->200), mult. (34139->246), div. (0->0), fcn. (21126->8), ass. (0->92)
t106 = sin(qJ(2));
t109 = cos(qJ(2));
t131 = qJD(1) * qJD(2);
t127 = t109 * t131;
t112 = qJD(1) ^ 2;
t107 = sin(qJ(1));
t110 = cos(qJ(1));
t135 = t107 * g(1) - t110 * g(2);
t70 = -qJDD(1) * pkin(1) - t112 * pkin(7) - t135;
t79 = qJDD(1) * t106 + t127;
t126 = t106 * t131;
t80 = qJDD(1) * t109 - t126;
t120 = -t80 * pkin(2) + t70 + (-t127 - t79) * qJ(3);
t104 = sin(qJ(5));
t147 = cos(qJ(5));
t133 = qJD(1) * t106;
t134 = t109 ^ 2 * t112;
t150 = 2 * qJD(3);
t86 = -qJD(2) * pkin(3) - pkin(8) * t133;
t113 = -pkin(2) * t126 + pkin(3) * t80 - pkin(8) * t134 - t120 + (t150 + t86) * t133;
t105 = sin(qJ(4));
t108 = cos(qJ(4));
t69 = (-t105 * t109 + t106 * t108) * qJD(1);
t47 = -qJD(4) * t69 - t105 * t79 - t108 * t80;
t68 = (t105 * t106 + t108 * t109) * qJD(1);
t48 = -qJD(4) * t68 - t105 * t80 + t108 * t79;
t98 = -qJD(2) + qJD(4);
t19 = (t69 * t98 - t47) * pkin(4) + (t68 * t98 - t48) * pkin(9) + t113;
t111 = qJD(2) ^ 2;
t124 = -g(1) * t110 - g(2) * t107;
t71 = -pkin(1) * t112 + qJDD(1) * pkin(7) + t124;
t128 = -g(3) * t106 + t109 * t71;
t132 = qJD(1) * t109;
t76 = (-pkin(2) * t109 - qJ(3) * t106) * qJD(1);
t116 = -pkin(2) * t111 + qJDD(2) * qJ(3) + qJD(2) * t150 + t76 * t132 + t128;
t31 = -pkin(3) * t134 - pkin(8) * t80 + qJD(2) * t86 + t116;
t138 = -t109 * g(3) - t106 * t71;
t46 = -qJDD(2) * pkin(2) - qJ(3) * t111 + t76 * t133 + qJDD(3) - t138;
t32 = (-t79 + t127) * pkin(8) + (-t106 * t109 * t112 - qJDD(2)) * pkin(3) + t46;
t141 = t105 * t32 + t108 * t31;
t56 = pkin(4) * t68 - pkin(9) * t69;
t96 = t98 ^ 2;
t97 = -qJDD(2) + qJDD(4);
t22 = -pkin(4) * t96 + pkin(9) * t97 - t56 * t68 + t141;
t142 = t104 * t19 + t147 * t22;
t57 = t104 * t69 - t147 * t98;
t58 = t104 * t98 + t147 * t69;
t37 = pkin(5) * t57 - qJ(6) * t58;
t45 = qJDD(5) - t47;
t67 = qJD(5) + t68;
t52 = -mrSges(7,1) * t67 + mrSges(7,2) * t58;
t63 = t67 ^ 2;
t130 = m(7) * (-pkin(5) * t63 + qJ(6) * t45 + 0.2e1 * qJD(6) * t67 - t37 * t57 + t142) + t67 * t52 + t45 * mrSges(7,3);
t38 = mrSges(7,1) * t57 - mrSges(7,3) * t58;
t140 = -mrSges(6,1) * t57 - mrSges(6,2) * t58 - t38;
t143 = -mrSges(6,3) - mrSges(7,2);
t26 = qJD(5) * t58 + t104 * t48 - t147 * t97;
t51 = mrSges(6,1) * t67 - mrSges(6,3) * t58;
t11 = m(6) * t142 - t45 * mrSges(6,2) + t140 * t57 + t143 * t26 - t67 * t51 + t130;
t122 = -t104 * t22 + t147 * t19;
t148 = m(7) * (-t45 * pkin(5) - t63 * qJ(6) + t58 * t37 + qJDD(6) - t122);
t27 = -t57 * qJD(5) + t104 * t97 + t147 * t48;
t49 = -mrSges(7,2) * t57 + mrSges(7,3) * t67;
t50 = -mrSges(6,2) * t67 - mrSges(6,3) * t57;
t13 = m(6) * t122 - t148 + (t50 + t49) * t67 + t140 * t58 + (mrSges(6,1) + mrSges(7,1)) * t45 + t143 * t27;
t59 = -mrSges(5,2) * t98 - mrSges(5,3) * t68;
t60 = mrSges(5,1) * t98 - mrSges(5,3) * t69;
t121 = m(5) * t113 - t47 * mrSges(5,1) + t48 * mrSges(5,2) + t104 * t11 + t147 * t13 + t68 * t59 + t69 * t60;
t117 = m(4) * ((pkin(2) * qJD(2) - (2 * qJD(3))) * t133 + t120) - t80 * mrSges(4,1) - t121;
t85 = mrSges(4,2) * t132 + qJD(2) * mrSges(4,3);
t136 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t132 + t85;
t82 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t133;
t83 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t133;
t154 = ((t82 - t83) * t106 - t136 * t109) * qJD(1) + m(3) * t70 - t80 * mrSges(3,1) + (mrSges(3,2) - mrSges(4,3)) * t79 + t117;
t125 = -t105 * t31 + t108 * t32;
t21 = -pkin(4) * t97 - pkin(9) * t96 + t69 * t56 - t125;
t129 = m(7) * (-0.2e1 * qJD(6) * t58 + (t57 * t67 - t27) * qJ(6) + (t58 * t67 + t26) * pkin(5) + t21) + t26 * mrSges(7,1) + t57 * t49;
t151 = m(6) * t21 + t26 * mrSges(6,1) + (t51 - t52) * t58 + (mrSges(6,2) - mrSges(7,3)) * t27 + t57 * t50 + t129;
t77 = (-mrSges(4,1) * t109 - mrSges(4,3) * t106) * qJD(1);
t55 = mrSges(5,1) * t68 + mrSges(5,2) * t69;
t8 = m(5) * t141 - t97 * mrSges(5,2) + t47 * mrSges(5,3) - t104 * t13 + t147 * t11 - t68 * t55 - t98 * t60;
t9 = m(5) * t125 + t97 * mrSges(5,1) - t48 * mrSges(5,3) - t69 * t55 + t98 * t59 - t151;
t119 = m(4) * t116 + qJDD(2) * mrSges(4,3) + qJD(2) * t83 - t105 * t9 + t108 * t8 + t77 * t132;
t144 = mrSges(3,3) + mrSges(4,2);
t78 = (-mrSges(3,1) * t109 + mrSges(3,2) * t106) * qJD(1);
t4 = m(3) * t128 - qJDD(2) * mrSges(3,2) - qJD(2) * t82 + t78 * t132 + t144 * t80 + t119;
t118 = -m(4) * t46 - t105 * t8 - t108 * t9;
t5 = m(3) * t138 - t144 * t79 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t136 * qJD(2) + (-t77 - t78) * t133 + t118;
t149 = t106 * t4 + t109 * t5;
t6 = m(2) * t135 + qJDD(1) * mrSges(2,1) - t112 * mrSges(2,2) - t154;
t1 = m(2) * t124 - t112 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t5 + t109 * t4;
t2 = [-m(1) * g(1) + t1 * t110 - t107 * t6, t1, t4, t80 * mrSges(4,2) + t119, t8, t11, -t26 * mrSges(7,2) - t57 * t38 + t130; -m(1) * g(2) + t1 * t107 + t110 * t6, t6, t5, -t79 * mrSges(4,3) + (-t106 * t83 - t109 * t85) * qJD(1) + t117, t9, t13, -t27 * mrSges(7,3) - t58 * t52 + t129; (-m(1) - m(2)) * g(3) + t149, -m(2) * g(3) + t149, t154, -qJDD(2) * mrSges(4,1) + t79 * mrSges(4,2) - qJD(2) * t85 + t77 * t133 - t118, t121, t151, -t45 * mrSges(7,1) + t27 * mrSges(7,2) + t58 * t38 - t67 * t49 + t148;];
f_new  = t2;
