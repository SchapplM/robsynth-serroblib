% Calculate vector of cutting forces with Newton-Euler
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-05-07 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:30:48
% EndTime: 2019-05-07 04:30:52
% DurationCPUTime: 1.41s
% Computational Cost: add. (14096->207), mult. (29388->245), div. (0->0), fcn. (19042->8), ass. (0->97)
t105 = qJDD(2) + qJDD(3);
t106 = qJD(2) + qJD(3);
t112 = sin(qJ(3));
t152 = cos(qJ(3));
t113 = sin(qJ(2));
t116 = cos(qJ(2));
t118 = qJD(1) ^ 2;
t139 = qJD(1) * qJD(2);
t114 = sin(qJ(1));
t117 = cos(qJ(1));
t136 = -t117 * g(1) - t114 * g(2);
t86 = -t118 * pkin(1) + qJDD(1) * pkin(7) + t136;
t144 = t113 * t86;
t91 = t113 * qJDD(1) + t116 * t139;
t35 = qJDD(2) * pkin(2) - t91 * pkin(8) - t144 + (pkin(2) * t113 * t118 + pkin(8) * t139 - g(3)) * t116;
t138 = -t113 * g(3) + t116 * t86;
t143 = t116 ^ 2 * t118;
t92 = t116 * qJDD(1) - t113 * t139;
t141 = qJD(1) * t113;
t95 = qJD(2) * pkin(2) - pkin(8) * t141;
t36 = -pkin(2) * t143 + t92 * pkin(8) - qJD(2) * t95 + t138;
t148 = t112 * t35 + t152 * t36;
t158 = t106 ^ 2;
t140 = qJD(1) * t116;
t83 = t112 * t141 - t152 * t140;
t84 = (t112 * t116 + t152 * t113) * qJD(1);
t62 = t83 * pkin(3) - t84 * qJ(4);
t134 = pkin(3) * t158 - t105 * qJ(4) + t83 * t62 - t148;
t74 = -t106 * mrSges(5,1) + t84 * mrSges(5,2);
t156 = 2 * qJD(4);
t97 = t106 * t156;
t160 = t105 * mrSges(5,3) + m(5) * (-t134 + t97) + t106 * t74;
t142 = t114 * g(1) - t117 * g(2);
t85 = -qJDD(1) * pkin(1) - t118 * pkin(7) - t142;
t131 = -t92 * pkin(2) - pkin(8) * t143 + t95 * t141 + t85;
t145 = t106 * t83;
t51 = t84 * qJD(3) + t112 * t91 - t152 * t92;
t52 = -t83 * qJD(3) + t112 * t92 + t152 * t91;
t127 = t51 * pkin(3) + t131 + (t145 - t52) * qJ(4);
t111 = sin(qJ(6));
t115 = cos(qJ(6));
t72 = -t106 * pkin(4) - t84 * qJ(5);
t82 = t83 ^ 2;
t123 = -t82 * qJ(5) + qJDD(5) - t127 + (t156 + t72) * t84;
t154 = -pkin(4) - pkin(9);
t14 = (-pkin(5) * t83 + (-pkin(3) - pkin(9)) * t84) * t106 + t154 * t51 + t52 * pkin(5) + t123;
t137 = -t112 * t36 + t152 * t35;
t128 = -qJ(4) * t158 + t84 * t62 + qJDD(4) - t137;
t155 = -2 * qJD(5);
t121 = (-t52 - t145) * qJ(5) + t128 + (t83 * pkin(4) + t155) * t84;
t65 = t84 * pkin(5) - t83 * pkin(9);
t15 = -t158 * pkin(5) - t84 * t65 + (-pkin(3) + t154) * t105 + t121;
t68 = -t115 * t106 - t111 * t83;
t29 = t68 * qJD(6) - t111 * t105 + t115 * t51;
t69 = -t111 * t106 + t115 * t83;
t34 = -t68 * mrSges(7,1) + t69 * mrSges(7,2);
t48 = qJDD(6) + t52;
t81 = qJD(6) + t84;
t54 = -t81 * mrSges(7,2) + t68 * mrSges(7,3);
t12 = m(7) * (-t111 * t15 + t115 * t14) - t29 * mrSges(7,3) + t48 * mrSges(7,1) - t69 * t34 + t81 * t54;
t28 = -t69 * qJD(6) - t115 * t105 - t111 * t51;
t55 = t81 * mrSges(7,1) - t69 * mrSges(7,3);
t13 = m(7) * (t111 * t14 + t115 * t15) + t28 * mrSges(7,3) - t48 * mrSges(7,2) + t68 * t34 - t81 * t55;
t151 = pkin(3) * t106;
t70 = -t106 * mrSges(6,1) - t83 * mrSges(6,3);
t75 = t106 * mrSges(6,2) - t84 * mrSges(6,3);
t133 = t111 * t13 + t115 * t12 + m(6) * (-t51 * pkin(4) - t84 * t151 + t123) + t51 * mrSges(6,2) + t52 * mrSges(6,1) + t83 * t70 + t84 * t75;
t157 = -2 * qJD(4);
t76 = -t83 * mrSges(5,2) + t106 * mrSges(5,3);
t124 = t52 * mrSges(5,3) + t84 * t74 + t133 - m(5) * ((t157 + t151) * t84 + t127) - t51 * mrSges(5,1) - t83 * t76;
t71 = -t106 * mrSges(4,2) - t83 * mrSges(4,3);
t73 = t106 * mrSges(4,1) - t84 * mrSges(4,3);
t120 = m(4) * t131 + t51 * mrSges(4,1) + t52 * mrSges(4,2) + t83 * t71 + t84 * t73 - t124;
t93 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t141;
t94 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t140;
t159 = m(3) * t85 - t92 * mrSges(3,1) + t91 * mrSges(3,2) + (t113 * t93 - t116 * t94) * qJD(1) + t120;
t61 = t84 * mrSges(6,1) + t83 * mrSges(6,2);
t132 = t111 * t12 - t115 * t13 - m(6) * ((-pkin(3) - pkin(4)) * t105 + t121) + t84 * t61 + t52 * mrSges(6,3);
t129 = m(5) * (-t105 * pkin(3) + t128) - t132;
t146 = t70 - t76;
t63 = t83 * mrSges(5,1) - t84 * mrSges(5,3);
t147 = -t83 * mrSges(4,1) - t84 * mrSges(4,2) - t63;
t149 = -mrSges(4,3) - mrSges(5,2);
t150 = -mrSges(5,1) + mrSges(6,2);
t7 = m(4) * t137 + t147 * t84 + t149 * t52 + (t71 - t146) * t106 + (mrSges(4,1) - t150) * t105 - t129;
t126 = t82 * pkin(4) - t51 * qJ(5) + t134;
t130 = -t28 * mrSges(7,1) - t68 * t54 + m(7) * (t105 * pkin(5) - t158 * pkin(9) + t106 * t72 + t97 + ((2 * qJD(5)) + t65) * t83 - t126) + t29 * mrSges(7,2) + t69 * t55;
t125 = -m(6) * (t83 * t155 + (t157 - t72) * t106 + t126) + t130 + t83 * t61 + t51 * mrSges(6,3);
t8 = t147 * t83 + (-t73 + t75) * t106 + t149 * t51 + (-mrSges(4,2) + mrSges(6,1)) * t105 + m(4) * t148 + t125 + t160;
t90 = (-mrSges(3,1) * t116 + mrSges(3,2) * t113) * qJD(1);
t4 = m(3) * (-t116 * g(3) - t144) - t91 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t90 * t141 + qJD(2) * t94 + t112 * t8 + t152 * t7;
t5 = m(3) * t138 - qJDD(2) * mrSges(3,2) + t92 * mrSges(3,3) - qJD(2) * t93 - t112 * t7 + t90 * t140 + t152 * t8;
t153 = t113 * t5 + t116 * t4;
t122 = -t105 * mrSges(6,1) - t106 * t75 - t125;
t6 = m(2) * t142 + qJDD(1) * mrSges(2,1) - t118 * mrSges(2,2) - t159;
t1 = m(2) * t136 - t118 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t113 * t4 + t116 * t5;
t2 = [-m(1) * g(1) + t117 * t1 - t114 * t6, t1, t5, t8, -t51 * mrSges(5,2) - t83 * t63 - t122 + t160, t105 * mrSges(6,2) + t106 * t70 - t132, t13; -m(1) * g(2) + t114 * t1 + t117 * t6, t6, t4, t7, -t124, t122, t12; (-m(1) - m(2)) * g(3) + t153, -m(2) * g(3) + t153, t159, t120, t52 * mrSges(5,2) + t150 * t105 + t146 * t106 + t84 * t63 + t129, t133, t130;];
f_new  = t2;
