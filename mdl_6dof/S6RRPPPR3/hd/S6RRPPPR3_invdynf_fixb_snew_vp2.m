% Calculate vector of cutting forces with Newton-Euler
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-05-06 08:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:33:32
% EndTime: 2019-05-06 08:33:35
% DurationCPUTime: 1.36s
% Computational Cost: add. (12404->205), mult. (27087->252), div. (0->0), fcn. (14414->8), ass. (0->97)
t115 = sin(qJ(2));
t118 = cos(qJ(2));
t146 = qJD(1) * qJD(2);
t140 = t118 * t146;
t121 = qJD(1) ^ 2;
t116 = sin(qJ(1));
t119 = cos(qJ(1));
t149 = t116 * g(1) - t119 * g(2);
t58 = -qJDD(1) * pkin(1) - t121 * pkin(7) - t149;
t80 = qJDD(1) * t115 + t140;
t138 = t115 * t146;
t81 = qJDD(1) * t118 - t138;
t133 = -t81 * pkin(2) + t58 + (-t140 - t80) * qJ(3);
t111 = sin(pkin(9));
t112 = cos(pkin(9));
t148 = qJD(1) * t115;
t151 = t118 ^ 2 * t121;
t86 = -qJD(2) * pkin(3) - qJ(4) * t148;
t126 = -qJ(4) * t151 + qJDD(4) - t133 + ((2 * qJD(3)) + t86) * t148;
t147 = qJD(1) * t118;
t114 = sin(qJ(6));
t117 = cos(qJ(6));
t163 = 2 * qJD(5);
t157 = pkin(3) + qJ(5);
t158 = -pkin(2) - qJ(5);
t20 = t126 + t157 * t81 + pkin(4) * t80 + (pkin(4) * t118 + t158 * t115) * t146;
t120 = qJD(2) ^ 2;
t135 = -g(1) * t119 - g(2) * t116;
t59 = -pkin(1) * t121 + qJDD(1) * pkin(7) + t135;
t152 = -t118 * g(3) - t115 * t59;
t75 = (-pkin(2) * t118 - qJ(3) * t115) * qJD(1);
t137 = t75 * t148 + qJDD(3) - t152;
t150 = t118 * t121;
t145 = qJD(1) * qJD(4);
t167 = -0.2e1 * t115 * t145 + (t140 - t80) * qJ(4);
t78 = (pkin(4) * t115 + qJ(5) * t118) * qJD(1);
t24 = (-pkin(4) - qJ(3)) * t120 + (-pkin(3) * t150 - qJD(1) * t78) * t115 + (-pkin(2) - t157) * qJDD(2) + t137 + t167;
t67 = qJD(2) * t111 + t112 * t147;
t136 = -t111 * t24 + t112 * t20 + t67 * t163;
t49 = -qJDD(2) * t111 - t112 * t81;
t66 = -qJD(2) * t112 + t111 * t147;
t14 = (t66 * t148 - t49) * pkin(8) + (-t66 * t67 + t80) * pkin(5) + t136;
t143 = t111 * t20 + t112 * t24 + t66 * t163;
t48 = -qJDD(2) * t112 + t111 * t81;
t50 = pkin(5) * t148 + t67 * pkin(8);
t63 = t66 ^ 2;
t15 = -t63 * pkin(5) + t48 * pkin(8) - t50 * t148 + t143;
t42 = t114 * t67 + t117 * t66;
t32 = t42 * qJD(6) + t114 * t48 + t117 * t49;
t43 = t114 * t66 - t117 * t67;
t34 = -mrSges(7,1) * t42 + mrSges(7,2) * t43;
t95 = qJD(6) + t148;
t39 = -mrSges(7,2) * t95 + t42 * mrSges(7,3);
t73 = qJDD(6) + t80;
t12 = m(7) * (-t114 * t15 + t117 * t14) - t32 * mrSges(7,3) + t73 * mrSges(7,1) - t43 * t34 + t95 * t39;
t31 = -t43 * qJD(6) - t114 * t49 + t117 * t48;
t40 = mrSges(7,1) * t95 - t43 * mrSges(7,3);
t13 = m(7) * (t114 * t14 + t117 * t15) + t31 * mrSges(7,3) - t73 * mrSges(7,2) + t42 * t34 - t95 * t40;
t44 = -mrSges(6,1) * t66 - mrSges(6,2) * t67;
t46 = -mrSges(6,2) * t148 + t66 * mrSges(6,3);
t8 = m(6) * t136 + t80 * mrSges(6,1) - t49 * mrSges(6,3) + t114 * t13 + t117 * t12 + t46 * t148 + t67 * t44;
t87 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t148;
t47 = mrSges(6,1) * t148 + t67 * mrSges(6,3);
t9 = m(6) * t143 - t80 * mrSges(6,2) + t48 * mrSges(6,3) - t114 * t12 + t117 * t13 - t47 * t148 + t66 * t44;
t90 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t147;
t134 = t111 * t9 + t112 * t8 + m(5) * (-pkin(2) * t138 + pkin(3) * t81 + t126) + t87 * t148 - t90 * t147 + t80 * mrSges(5,1) - t81 * mrSges(5,2);
t130 = m(4) * ((pkin(2) * qJD(2) - (2 * qJD(3))) * t148 + t133) - t81 * mrSges(4,1) - t134;
t92 = mrSges(4,2) * t147 + qJD(2) * mrSges(4,3);
t154 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t147 + t92;
t88 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t148;
t89 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t148;
t169 = ((t88 - t89) * t115 - t154 * t118) * qJD(1) + m(3) * t58 - t81 * mrSges(3,1) + (mrSges(3,2) - mrSges(4,3)) * t80 + t130;
t141 = -g(3) * t115 + t118 * t59;
t166 = qJDD(2) * qJ(3) + t75 * t147 + t141;
t168 = pkin(3) * t151 + t81 * qJ(4) - qJD(2) * t86 + 0.2e1 * t118 * t145 - t166;
t144 = qJD(3) * qJD(2);
t100 = 0.2e1 * t144;
t123 = qJDD(2) * pkin(4) + t158 * t120 - t78 * t147 + qJDD(5) + t100 - t168;
t131 = -t31 * mrSges(7,1) - t42 * t39 + m(7) * (-t48 * pkin(5) - t63 * pkin(8) - t67 * t50 + t123) + t32 * mrSges(7,2) + t43 * t40;
t125 = m(6) * t123 - t48 * mrSges(6,1) + t49 * mrSges(6,2) - t66 * t46 - t67 * t47 + t131;
t161 = pkin(2) * t120;
t124 = -m(5) * (-0.2e1 * t144 + t161 + t168) + qJDD(2) * mrSges(5,1) + t125 - t81 * mrSges(5,3) + qJD(2) * t87;
t76 = (-mrSges(4,1) * t118 - mrSges(4,3) * t115) * qJD(1);
t122 = qJDD(2) * mrSges(4,3) + t124 + m(4) * (t100 - t161 + t166) + t76 * t147 + qJD(2) * t89;
t79 = (mrSges(5,1) * t115 - mrSges(5,2) * t118) * qJD(1);
t156 = (-mrSges(3,1) * t118 + mrSges(3,2) * t115) * qJD(1) - t79;
t159 = mrSges(3,3) + mrSges(4,2);
t11 = m(3) * t141 - qJDD(2) * mrSges(3,2) - qJD(2) * t88 + t156 * t147 + t159 * t81 + t122;
t38 = -qJDD(2) * pkin(2) - t120 * qJ(3) + t137;
t132 = t111 * t8 - t112 * t9 - qJDD(2) * mrSges(5,2) - m(5) * ((-t115 * t150 - qJDD(2)) * pkin(3) + t38 + t167) + t80 * mrSges(5,3) - qJD(2) * t90;
t129 = m(4) * t38 - t132;
t4 = m(3) * t152 - t159 * t80 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t154 * qJD(2) + (-t76 - t156) * t148 - t129;
t162 = t115 * t11 + t118 * t4;
t142 = t79 * t147;
t2 = m(2) * t149 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t169;
t1 = m(2) * t135 - t121 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t118 * t11 - t115 * t4;
t3 = [-m(1) * g(1) + t1 * t119 - t116 * t2, t1, t11, t81 * mrSges(4,2) + t122 - t142, -t79 * t148 - t132, t9, t13; -m(1) * g(2) + t1 * t116 + t119 * t2, t2, t4, -t80 * mrSges(4,3) + (-t115 * t89 - t118 * t92) * qJD(1) + t130, -t124 + t142, t8, t12; (-m(1) - m(2)) * g(3) + t162, -m(2) * g(3) + t162, t169, -qJDD(2) * mrSges(4,1) + t80 * mrSges(4,2) - qJD(2) * t92 + (t76 - t79) * t148 + t129, t134, t125, t131;];
f_new  = t3;
