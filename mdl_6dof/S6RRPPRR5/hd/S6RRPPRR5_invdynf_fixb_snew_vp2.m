% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRR5
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
% Datum: 2019-05-06 10:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:38:12
% EndTime: 2019-05-06 10:38:19
% DurationCPUTime: 1.98s
% Computational Cost: add. (20394->220), mult. (47014->276), div. (0->0), fcn. (31801->10), ass. (0->105)
t116 = cos(pkin(6));
t107 = t116 * qJDD(1) + qJDD(2);
t108 = t116 * qJD(1) + qJD(2);
t115 = sin(pkin(6));
t123 = cos(qJ(2));
t152 = qJD(1) * t123;
t145 = t115 * t152;
t119 = sin(qJ(2));
t125 = qJD(1) ^ 2;
t120 = sin(qJ(1));
t124 = cos(qJ(1));
t144 = t120 * g(1) - t124 * g(2);
t81 = t125 * t115 * pkin(8) + qJDD(1) * pkin(1) + t144;
t158 = t116 * t81;
t140 = -t124 * g(1) - t120 * g(2);
t150 = qJDD(1) * t115;
t82 = -t125 * pkin(1) + pkin(8) * t150 + t140;
t162 = t119 * t158 + t123 * t82;
t168 = t108 ^ 2;
t153 = qJD(1) * t115;
t83 = (-pkin(2) * t123 - qJ(3) * t119) * t153;
t136 = -pkin(2) * t168 + t107 * qJ(3) + 0.2e1 * qJD(3) * t108 + t83 * t145 + t162;
t155 = t115 * t119;
t149 = g(3) * t155;
t146 = t119 * t153;
t77 = -t108 * mrSges(4,1) + mrSges(4,2) * t146;
t84 = (-mrSges(4,1) * t123 - mrSges(4,3) * t119) * t153;
t169 = m(4) * (t136 - t149) + t108 * t77 + t84 * t145 + t107 * mrSges(4,3);
t167 = -2 * qJD(4);
t166 = pkin(2) * t108;
t165 = -mrSges(5,2) - mrSges(4,3);
t164 = mrSges(3,3) + mrSges(4,2);
t118 = sin(qJ(5));
t122 = cos(qJ(5));
t142 = t108 * t145;
t157 = -t116 * g(3) - t115 * t81;
t88 = (qJD(2) * t152 + qJDD(1) * t119) * t115;
t89 = -qJD(2) * t146 + t123 * t150;
t141 = -t89 * pkin(2) + t157 + (-t142 - t88) * qJ(3);
t143 = 0.2e1 * t146;
t156 = t115 ^ 2 * t125;
t147 = t123 ^ 2 * t156;
t74 = -t108 * pkin(3) - qJ(4) * t146;
t129 = -qJ(4) * t147 + qJD(3) * t143 + t74 * t146 + qJDD(4) - t141;
t20 = -t88 * pkin(9) + (pkin(3) + pkin(4)) * t89 + (-pkin(9) * t123 + (-pkin(2) - pkin(4)) * t119) * t108 * t153 + t129;
t128 = -pkin(3) * t147 - t89 * qJ(4) + t108 * t74 + t145 * t167 + t136;
t87 = (pkin(4) * t123 - pkin(9) * t119) * t153;
t22 = -t168 * pkin(4) - t107 * pkin(9) + (-g(3) * t119 - t87 * t152) * t115 + t128;
t163 = t118 * t20 + t122 * t22;
t75 = -t108 * mrSges(5,1) - mrSges(5,3) * t146;
t161 = -t108 * mrSges(3,1) + mrSges(3,3) * t146 + t75;
t80 = mrSges(4,2) * t145 + t108 * mrSges(4,3);
t160 = -t108 * mrSges(3,2) + mrSges(3,3) * t145 + t80;
t85 = (mrSges(5,1) * t123 + mrSges(5,2) * t119) * t153;
t159 = -t85 + (-mrSges(3,1) * t123 + mrSges(3,2) * t119) * t153;
t154 = t115 * t123;
t117 = sin(qJ(6));
t121 = cos(qJ(6));
t62 = -t122 * t108 - t118 * t146;
t63 = -t118 * t108 + t122 * t146;
t45 = -t62 * pkin(5) - t63 * pkin(10);
t72 = qJDD(5) + t89;
t94 = qJD(5) + t145;
t92 = t94 ^ 2;
t17 = -t92 * pkin(5) + t72 * pkin(10) + t62 * t45 + t163;
t148 = -g(3) * t154 - t119 * t82 + t123 * t158;
t33 = -t107 * pkin(2) - qJ(3) * t168 + t83 * t146 + qJDD(3) - t148;
t133 = -t88 * qJ(4) + t33 + (-t119 * t123 * t156 - t107) * pkin(3);
t126 = t107 * pkin(4) - pkin(9) * t168 - qJ(4) * t142 + qJD(4) * t143 + t87 * t146 - t133;
t42 = -t63 * qJD(5) - t122 * t107 - t118 * t88;
t43 = t62 * qJD(5) - t118 * t107 + t122 * t88;
t18 = (-t62 * t94 - t43) * pkin(10) + (t63 * t94 - t42) * pkin(5) + t126;
t46 = -t117 * t63 + t121 * t94;
t30 = t46 * qJD(6) + t117 * t72 + t121 * t43;
t47 = t117 * t94 + t121 * t63;
t34 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t60 = qJD(6) - t62;
t35 = -t60 * mrSges(7,2) + t46 * mrSges(7,3);
t39 = qJDD(6) - t42;
t14 = m(7) * (-t117 * t17 + t121 * t18) - t30 * mrSges(7,3) + t39 * mrSges(7,1) - t47 * t34 + t60 * t35;
t29 = -t47 * qJD(6) - t117 * t43 + t121 * t72;
t36 = t60 * mrSges(7,1) - t47 * mrSges(7,3);
t15 = m(7) * (t117 * t18 + t121 * t17) + t29 * mrSges(7,3) - t39 * mrSges(7,2) + t46 * t34 - t60 * t36;
t44 = -t62 * mrSges(6,1) + t63 * mrSges(6,2);
t49 = t94 * mrSges(6,1) - t63 * mrSges(6,3);
t10 = m(6) * t163 - t72 * mrSges(6,2) + t42 * mrSges(6,3) - t117 * t14 + t121 * t15 + t62 * t44 - t94 * t49;
t139 = -t118 * t22 + t122 * t20;
t130 = m(7) * (-t72 * pkin(5) - t92 * pkin(10) + t63 * t45 - t139) - t29 * mrSges(7,1) + t30 * mrSges(7,2) - t46 * t35 + t47 * t36;
t48 = -t94 * mrSges(6,2) + t62 * mrSges(6,3);
t11 = m(6) * t139 + t72 * mrSges(6,1) - t43 * mrSges(6,3) - t63 * t44 + t94 * t48 - t130;
t138 = t122 * t10 - t118 * t11 + m(5) * (t128 - t149) - t89 * mrSges(5,3);
t4 = m(3) * (-t149 + t162) + t164 * t89 + t161 * t108 + (-mrSges(3,2) + mrSges(5,2)) * t107 + t159 * t145 + t138 + t169;
t135 = m(5) * (t89 * pkin(3) - t146 * t166 + t129) + t118 * t10 + t122 * t11 + t89 * mrSges(5,1);
t132 = -t135 + m(4) * ((-0.2e1 * qJD(3) + t166) * t146 + t141) - t89 * mrSges(4,1);
t78 = t108 * mrSges(5,2) - mrSges(5,3) * t145;
t6 = m(3) * t157 - t89 * mrSges(3,1) + (mrSges(3,2) + t165) * t88 + ((-t78 - t160) * t123 + (-t77 - t161) * t119) * t153 + t132;
t137 = m(6) * t126 - t42 * mrSges(6,1) + t43 * mrSges(6,2) + t117 * t15 + t121 * t14 - t62 * t48 + t63 * t49;
t134 = m(5) * ((qJ(4) * t108 * t123 + t119 * t167) * t153 + t133) - t107 * mrSges(5,1) - t108 * t78 - t137;
t131 = m(4) * t33 + t134;
t8 = -t131 + (mrSges(5,3) - t164) * t88 + t160 * t108 + (mrSges(3,1) + mrSges(4,1)) * t107 + m(3) * t148 + (-t84 - t159) * t146;
t151 = t116 * t6 + t8 * t154 + t4 * t155;
t127 = t107 * mrSges(5,2) + t108 * t75 - t85 * t145 + t138;
t2 = m(2) * t140 - t125 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t119 * t8 + t123 * t4;
t1 = m(2) * t144 + qJDD(1) * mrSges(2,1) - t125 * mrSges(2,2) - t115 * t6 + (t119 * t4 + t123 * t8) * t116;
t3 = [-m(1) * g(1) - t120 * t1 + t124 * t2, t2, t4, t89 * mrSges(4,2) + t127 + t169, t127, t10, t15; -m(1) * g(2) + t124 * t1 + t120 * t2, t1, t8, t165 * t88 + ((-t78 - t80) * t123 + (-t75 - t77) * t119) * t153 + t132, -t88 * mrSges(5,3) - t85 * t146 + t134, t11, t14; (-m(1) - m(2)) * g(3) + t151, -m(2) * g(3) + t151, t6, (mrSges(4,2) - mrSges(5,3)) * t88 + t131 - t107 * mrSges(4,1) - t108 * t80 + (t84 - t85) * t146, t88 * mrSges(5,2) + (t119 * t75 + t123 * t78) * t153 + t135, t137, t130;];
f_new  = t3;
