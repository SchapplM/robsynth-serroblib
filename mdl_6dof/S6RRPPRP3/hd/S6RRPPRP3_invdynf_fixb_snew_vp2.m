% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-05-06 09:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:14:24
% EndTime: 2019-05-06 09:14:27
% DurationCPUTime: 1.02s
% Computational Cost: add. (6203->203), mult. (12829->236), div. (0->0), fcn. (6318->6), ass. (0->89)
t114 = qJD(2) ^ 2;
t112 = cos(qJ(2));
t143 = qJD(1) * t112;
t161 = -pkin(2) - pkin(8);
t140 = qJD(1) * qJD(4);
t115 = qJD(1) ^ 2;
t145 = t112 ^ 2 * t115;
t109 = sin(qJ(2));
t110 = sin(qJ(1));
t113 = cos(qJ(1));
t128 = -t113 * g(1) - t110 * g(2);
t54 = -t115 * pkin(1) + qJDD(1) * pkin(7) + t128;
t134 = -t109 * g(3) + t112 * t54;
t71 = (-pkin(2) * t112 - qJ(3) * t109) * qJD(1);
t166 = qJDD(2) * qJ(3) + t71 * t143 + t134;
t141 = qJD(1) * qJD(2);
t132 = t109 * t141;
t77 = t112 * qJDD(1) - t132;
t142 = t109 * qJD(1);
t82 = -qJD(2) * pkin(3) - qJ(4) * t142;
t168 = pkin(3) * t145 + t77 * qJ(4) - qJD(2) * t82 + 0.2e1 * t112 * t140 - t166;
t75 = (pkin(4) * t109 + pkin(8) * t112) * qJD(1);
t139 = qJD(3) * qJD(2);
t96 = 0.2e1 * t139;
t117 = qJDD(2) * pkin(4) + t161 * t114 - t75 * t143 - t168 + t96;
t108 = sin(qJ(5));
t111 = cos(qJ(5));
t70 = t108 * qJD(2) + t111 * t143;
t38 = t70 * qJD(5) - t111 * qJDD(2) + t108 * t77;
t69 = -t111 * qJD(2) + t108 * t143;
t39 = t69 * qJD(5) - t108 * qJDD(2) - t111 * t77;
t91 = qJD(5) + t142;
t45 = t91 * pkin(5) + t70 * qJ(6);
t46 = t91 * mrSges(7,1) + t70 * mrSges(7,3);
t66 = t69 ^ 2;
t136 = t70 * t46 - m(7) * (-t38 * pkin(5) - t66 * qJ(6) - t70 * t45 + qJDD(6) + t117) - t39 * mrSges(7,2);
t43 = -t91 * mrSges(7,2) + t69 * mrSges(7,3);
t44 = -t91 * mrSges(6,2) + t69 * mrSges(6,3);
t47 = t91 * mrSges(6,1) + t70 * mrSges(6,3);
t119 = m(6) * t117 + t39 * mrSges(6,2) - (t43 + t44) * t69 - (mrSges(6,1) + mrSges(7,1)) * t38 - t70 * t47 - t136;
t157 = t114 * pkin(2);
t83 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t142;
t170 = -t119 - qJDD(2) * mrSges(5,1) + t77 * mrSges(5,3) - qJD(2) * t83 + m(5) * (-0.2e1 * t139 + t157 + t168);
t131 = t112 * t141;
t147 = t110 * g(1) - t113 * g(2);
t53 = -qJDD(1) * pkin(1) - t115 * pkin(7) - t147;
t76 = t109 * qJDD(1) + t131;
t125 = -t77 * pkin(2) + t53 + (-t131 - t76) * qJ(3);
t118 = -qJ(4) * t145 + qJDD(4) - t125 + (0.2e1 * qJD(3) + t82) * t142;
t160 = pkin(3) + pkin(8);
t18 = t118 + t160 * t77 + (pkin(4) * t112 + t161 * t109) * t141 + t76 * pkin(4);
t151 = -t112 * g(3) - t109 * t54;
t129 = t71 * t142 + qJDD(3) - t151;
t144 = t112 * t115;
t167 = -0.2e1 * t109 * t140 + (t131 - t76) * qJ(4);
t22 = (-pkin(4) - qJ(3)) * t114 + (-pkin(3) * t144 - qJD(1) * t75) * t109 + (-pkin(2) - t160) * qJDD(2) + t129 + t167;
t153 = t108 * t18 + t111 * t22;
t162 = 2 * qJD(6);
t41 = -t69 * mrSges(7,1) - t70 * mrSges(7,2);
t137 = m(7) * (-t66 * pkin(5) + t38 * qJ(6) + t69 * t162 - t91 * t45 + t153) + t69 * t41 + t38 * mrSges(7,3);
t42 = -t69 * mrSges(6,1) - t70 * mrSges(6,2);
t67 = qJDD(5) + t76;
t11 = m(6) * t153 + t38 * mrSges(6,3) + t69 * t42 + (-t47 - t46) * t91 + (-mrSges(6,2) - mrSges(7,2)) * t67 + t137;
t130 = -t108 * t22 + t111 * t18;
t138 = m(7) * (t70 * t162 + (t69 * t91 - t39) * qJ(6) + (-t69 * t70 + t67) * pkin(5) + t130) + t91 * t43 + t67 * mrSges(7,1);
t8 = m(6) * t130 + t67 * mrSges(6,1) + t91 * t44 + (t42 + t41) * t70 + (-mrSges(6,3) - mrSges(7,3)) * t39 + t138;
t86 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t143;
t126 = t108 * t11 + t111 * t8 + m(5) * (-pkin(2) * t132 + t77 * pkin(3) + t118) + t83 * t142 - t86 * t143 + t76 * mrSges(5,1) - t77 * mrSges(5,2);
t123 = m(4) * ((pkin(2) * qJD(2) - 0.2e1 * qJD(3)) * t142 + t125) - t77 * mrSges(4,1) - t126;
t88 = mrSges(4,2) * t143 + qJD(2) * mrSges(4,3);
t148 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t143 + t88;
t84 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t142;
t85 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t142;
t169 = (-t148 * t112 + (t84 - t85) * t109) * qJD(1) + m(3) * t53 - t77 * mrSges(3,1) + (mrSges(3,2) - mrSges(4,3)) * t76 + t123;
t31 = -qJDD(2) * pkin(2) - t114 * qJ(3) + t129;
t124 = t108 * t8 - t111 * t11 - qJDD(2) * mrSges(5,2) - m(5) * ((-t109 * t144 - qJDD(2)) * pkin(3) + t31 + t167) + t76 * mrSges(5,3) - qJD(2) * t86;
t122 = m(4) * t31 - t124;
t74 = (mrSges(5,1) * t109 - mrSges(5,2) * t112) * qJD(1);
t150 = (-mrSges(3,1) * t112 + mrSges(3,2) * t109) * qJD(1) - t74;
t154 = mrSges(3,3) + mrSges(4,2);
t72 = (-mrSges(4,1) * t112 - mrSges(4,3) * t109) * qJD(1);
t4 = m(3) * t151 - t154 * t76 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t148 * qJD(2) + (-t72 - t150) * t142 - t122;
t116 = qJDD(2) * mrSges(4,3) + m(4) * (t96 - t157 + t166) + t72 * t143 + qJD(2) * t85 - t170;
t6 = m(3) * t134 - qJDD(2) * mrSges(3,2) - qJD(2) * t84 + t150 * t143 + t154 * t77 + t116;
t159 = t109 * t6 + t112 * t4;
t135 = t74 * t143;
t2 = m(2) * t147 + qJDD(1) * mrSges(2,1) - t115 * mrSges(2,2) - t169;
t1 = m(2) * t128 - t115 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t109 * t4 + t112 * t6;
t3 = [-m(1) * g(1) + t113 * t1 - t110 * t2, t1, t6, t77 * mrSges(4,2) + t116 - t135, -t74 * t142 - t124, t11, -t67 * mrSges(7,2) - t91 * t46 + t137; -m(1) * g(2) + t110 * t1 + t113 * t2, t2, t4, -t76 * mrSges(4,3) + (-t109 * t85 - t112 * t88) * qJD(1) + t123, t135 + t170, t8, -t39 * mrSges(7,3) + t70 * t41 + t138; (-m(1) - m(2)) * g(3) + t159, -m(2) * g(3) + t159, t169, -qJDD(2) * mrSges(4,1) + t76 * mrSges(4,2) - qJD(2) * t88 + (t72 - t74) * t142 + t122, t126, t119, -t38 * mrSges(7,1) - t69 * t43 - t136;];
f_new  = t3;
