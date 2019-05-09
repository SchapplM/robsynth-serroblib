% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-05-06 11:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:02:50
% EndTime: 2019-05-06 11:02:55
% DurationCPUTime: 1.94s
% Computational Cost: add. (20345->224), mult. (46867->276), div. (0->0), fcn. (31801->10), ass. (0->112)
t108 = cos(pkin(6));
t100 = t108 * qJD(1) + qJD(2);
t107 = sin(pkin(6));
t115 = cos(qJ(2));
t148 = qJD(1) * t115;
t140 = t107 * t148;
t136 = t100 * t140;
t117 = qJD(1) ^ 2;
t112 = sin(qJ(1));
t116 = cos(qJ(1));
t138 = t112 * g(1) - t116 * g(2);
t77 = t117 * t107 * pkin(8) + qJDD(1) * pkin(1) + t138;
t153 = -t108 * g(3) - t107 * t77;
t111 = sin(qJ(2));
t84 = (qJD(2) * t148 + qJDD(1) * t111) * t107;
t149 = qJD(1) * t107;
t139 = t111 * t149;
t144 = qJDD(1) * t107;
t85 = -qJD(2) * t139 + t115 * t144;
t135 = -t85 * pkin(2) + t153 + (-t136 - t84) * qJ(3);
t163 = pkin(2) * t100;
t168 = m(4) * ((-(2 * qJD(3)) + t163) * t139 + t135) - t85 * mrSges(4,1);
t154 = t108 * t77;
t134 = -t116 * g(1) - t112 * g(2);
t78 = -t117 * pkin(1) + pkin(8) * t144 + t134;
t158 = t111 * t154 + t115 * t78;
t166 = t100 ^ 2;
t79 = (-pkin(2) * t115 - qJ(3) * t111) * t149;
t99 = t108 * qJDD(1) + qJDD(2);
t130 = pkin(2) * t166 - t99 * qJ(3) - t79 * t140 - t158;
t151 = t107 * t111;
t143 = g(3) * t151;
t73 = -t100 * mrSges(4,1) + mrSges(4,2) * t139;
t80 = (-mrSges(4,1) * t115 - mrSges(4,3) * t111) * t149;
t147 = qJD(3) * t100;
t89 = 0.2e1 * t147;
t167 = m(4) * (-t130 + t89 - t143) + t100 * t73 + t80 * t140 + t99 * mrSges(4,3);
t165 = pkin(3) + pkin(9);
t164 = t99 * pkin(2);
t162 = t111 * g(3);
t161 = -mrSges(4,1) + mrSges(5,2);
t160 = mrSges(3,3) + mrSges(4,2);
t110 = sin(qJ(5));
t114 = cos(qJ(5));
t152 = t107 ^ 2 * t117;
t142 = t115 ^ 2 * t152;
t70 = -t100 * pkin(3) - qJ(4) * t139;
t121 = -qJ(4) * t142 + qJDD(4) - t135 + ((2 * qJD(3)) + t70) * t139;
t20 = t84 * pkin(4) + t165 * t85 + (pkin(4) * t115 + (-pkin(2) - pkin(9)) * t111) * t100 * t149 + t121;
t150 = t107 * t115;
t131 = -g(3) * t150 - t111 * t78 + t115 * t154;
t126 = -qJ(3) * t166 + t79 * t139 + qJDD(3) - t131;
t145 = qJD(1) * qJD(4);
t137 = -0.2e1 * t107 * t145;
t120 = t111 * t137 + t126 + (t136 - t84) * qJ(4);
t141 = t115 * t152;
t83 = (pkin(4) * t111 + pkin(9) * t115) * t149;
t24 = -t166 * pkin(4) + (-pkin(3) * t141 - t83 * t149) * t111 + (-pkin(2) - t165) * t99 + t120;
t159 = t110 * t20 + t114 * t24;
t71 = t100 * mrSges(5,2) - mrSges(5,3) * t139;
t157 = -t100 * mrSges(3,1) + mrSges(3,3) * t139 + t71;
t76 = mrSges(4,2) * t140 + t100 * mrSges(4,3);
t156 = t100 * mrSges(3,2) - mrSges(3,3) * t140 - t76;
t82 = (mrSges(5,1) * t111 - mrSges(5,2) * t115) * t149;
t155 = (-mrSges(3,1) * t115 + mrSges(3,2) * t111) * t149 - t82;
t109 = sin(qJ(6));
t113 = cos(qJ(6));
t58 = -t114 * t100 + t110 * t140;
t59 = -t110 * t100 - t114 * t140;
t43 = -t58 * pkin(5) - t59 * pkin(10);
t68 = qJDD(5) + t84;
t90 = qJD(5) + t139;
t88 = t90 ^ 2;
t17 = -t88 * pkin(5) + t68 * pkin(10) + t58 * t43 + t159;
t123 = pkin(3) * t142 + t85 * qJ(4) - t100 * t70 + t130;
t118 = t99 * pkin(4) - t166 * pkin(9) + t89 + t115 * t137 + (-t83 * t148 - t162) * t107 - t123;
t40 = -t59 * qJD(5) + t110 * t85 - t114 * t99;
t41 = t58 * qJD(5) - t110 * t99 - t114 * t85;
t18 = (t59 * t90 - t40) * pkin(5) + (-t58 * t90 - t41) * pkin(10) + t118;
t44 = -t109 * t59 + t113 * t90;
t30 = t44 * qJD(6) + t109 * t68 + t113 * t41;
t45 = t109 * t90 + t113 * t59;
t34 = -t44 * mrSges(7,1) + t45 * mrSges(7,2);
t56 = qJD(6) - t58;
t35 = -t56 * mrSges(7,2) + t44 * mrSges(7,3);
t38 = qJDD(6) - t40;
t14 = m(7) * (-t109 * t17 + t113 * t18) - t30 * mrSges(7,3) + t38 * mrSges(7,1) - t45 * t34 + t56 * t35;
t29 = -t45 * qJD(6) - t109 * t41 + t113 * t68;
t36 = t56 * mrSges(7,1) - t45 * mrSges(7,3);
t15 = m(7) * (t109 * t18 + t113 * t17) + t29 * mrSges(7,3) - t38 * mrSges(7,2) + t44 * t34 - t56 * t36;
t42 = -t58 * mrSges(6,1) + t59 * mrSges(6,2);
t47 = t90 * mrSges(6,1) - t59 * mrSges(6,3);
t10 = m(6) * t159 - t68 * mrSges(6,2) + t40 * mrSges(6,3) - t109 * t14 + t113 * t15 + t58 * t42 - t90 * t47;
t133 = -t110 * t24 + t114 * t20;
t124 = m(7) * (-t68 * pkin(5) - t88 * pkin(10) + t59 * t43 - t133) - t29 * mrSges(7,1) + t30 * mrSges(7,2) - t44 * t35 + t45 * t36;
t46 = -t90 * mrSges(6,2) + t58 * mrSges(6,3);
t11 = m(6) * t133 + t68 * mrSges(6,1) - t41 * mrSges(6,3) - t59 * t42 + t90 * t46 - t124;
t74 = -t100 * mrSges(5,1) + mrSges(5,3) * t140;
t128 = -m(5) * (t85 * pkin(3) - t139 * t163 + t121) - t110 * t10 - t114 * t11 + t74 * t140 - t84 * mrSges(5,1);
t4 = m(3) * t153 + (-mrSges(3,1) + mrSges(5,2)) * t85 + (mrSges(3,2) - mrSges(4,3)) * t84 + (t156 * t115 + (-t73 - t157) * t111) * t149 + t128 + t168;
t132 = -t114 * t10 + t110 * t11 - m(5) * (-t164 + (-t111 * t141 - t99) * pkin(3) + t120) + t84 * mrSges(5,3);
t129 = m(4) * (t126 - t164) - t132;
t6 = m(3) * t131 - t160 * t84 + (mrSges(3,1) - t161) * t99 + (-t74 - t156) * t100 + (-t80 - t155) * t139 - t129;
t127 = m(6) * t118 - t40 * mrSges(6,1) + t41 * mrSges(6,2) + t109 * t15 + t113 * t14 - t58 * t46 + t59 * t47;
t122 = -m(5) * (-0.2e1 * t147 + (0.2e1 * t115 * t145 + t162) * t107 + t123) + t127 - t85 * mrSges(5,3);
t8 = t155 * t140 + t160 * t85 + (-mrSges(3,2) + mrSges(5,1)) * t99 + t157 * t100 + m(3) * (-t143 + t158) + t122 + t167;
t146 = t108 * t4 + t6 * t150 + t8 * t151;
t125 = t85 * mrSges(5,2) + t128;
t119 = -t99 * mrSges(5,1) - t100 * t71 + t82 * t140 - t122;
t2 = m(2) * t134 - t117 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t111 * t6 + t115 * t8;
t1 = m(2) * t138 + qJDD(1) * mrSges(2,1) - t117 * mrSges(2,2) - t107 * t4 + (t111 * t8 + t115 * t6) * t108;
t3 = [-m(1) * g(1) - t112 * t1 + t116 * t2, t2, t8, t85 * mrSges(4,2) - t119 + t167, t99 * mrSges(5,2) + t100 * t74 - t82 * t139 - t132, t10, t15; -m(1) * g(2) + t116 * t1 + t112 * t2, t1, t6, -t84 * mrSges(4,3) + (-t115 * t76 + (-t71 - t73) * t111) * t149 + t125 + t168, t119, t11, t14; (-m(1) - m(2)) * g(3) + t146, -m(2) * g(3) + t146, t4, t84 * mrSges(4,2) + t161 * t99 + (t74 - t76) * t100 + (t80 - t82) * t139 + t129, t71 * t139 - t125, t127, t124;];
f_new  = t3;
