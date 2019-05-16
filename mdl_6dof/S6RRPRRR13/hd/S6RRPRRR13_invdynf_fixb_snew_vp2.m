% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-07 01:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR13_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR13_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR13_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 01:08:22
% EndTime: 2019-05-07 01:08:37
% DurationCPUTime: 5.17s
% Computational Cost: add. (65156->217), mult. (144512->283), div. (0->0), fcn. (107217->12), ass. (0->114)
t162 = -2 * qJD(3);
t111 = cos(pkin(6));
t105 = t111 * qJD(1) + qJD(2);
t115 = sin(qJ(2));
t110 = sin(pkin(6));
t146 = qJD(1) * t110;
t138 = t115 * t146;
t161 = (pkin(2) * t105 + t162) * t138;
t103 = t105 ^ 2;
t104 = t111 * qJDD(1) + qJDD(2);
t120 = cos(qJ(2));
t148 = t110 * t115;
t122 = qJD(1) ^ 2;
t116 = sin(qJ(1));
t121 = cos(qJ(1));
t137 = t116 * g(1) - t121 * g(2);
t85 = t122 * t110 * pkin(8) + qJDD(1) * pkin(1) + t137;
t151 = t111 * t85;
t134 = -t121 * g(1) - t116 * g(2);
t142 = qJDD(1) * t110;
t86 = -t122 * pkin(1) + pkin(8) * t142 + t134;
t132 = -g(3) * t148 + t115 * t151 + t120 * t86;
t145 = qJD(1) * t120;
t139 = t110 * t145;
t87 = (-pkin(2) * t120 - qJ(3) * t115) * t146;
t160 = t103 * pkin(2) - t104 * qJ(3) + t105 * t162 - t87 * t139 - t132;
t159 = -pkin(2) - pkin(9);
t158 = t111 * g(3);
t157 = mrSges(3,1) - mrSges(4,2);
t156 = mrSges(3,3) + mrSges(4,1);
t113 = sin(qJ(5));
t118 = cos(qJ(5));
t114 = sin(qJ(4));
t119 = cos(qJ(4));
t149 = t110 ^ 2 * t122;
t140 = t120 ^ 2 * t149;
t90 = pkin(3) * t138 - t105 * pkin(9);
t91 = (qJD(2) * t145 + qJDD(1) * t115) * t110;
t92 = -qJD(2) * t138 + t120 * t142;
t37 = -pkin(3) * t140 - t158 - t91 * qJ(3) + t159 * t92 + (-t85 + (-qJ(3) * t105 * t120 - t115 * t90) * qJD(1)) * t110 + t161;
t147 = t110 * t120;
t150 = g(3) * t147 + t115 * t86;
t130 = -t103 * qJ(3) + t87 * t138 + qJDD(3) + t150;
t39 = t91 * pkin(3) + t159 * t104 + (-pkin(3) * t105 * t146 - pkin(9) * t115 * t149 - t151) * t120 + t130;
t154 = t114 * t39 + t119 * t37;
t75 = -t114 * t105 - t119 * t139;
t76 = t119 * t105 - t114 * t139;
t61 = -t75 * pkin(4) - t76 * pkin(10);
t80 = qJDD(4) + t91;
t99 = qJD(4) + t138;
t96 = t99 ^ 2;
t24 = -t96 * pkin(4) + t80 * pkin(10) + t75 * t61 + t154;
t124 = t92 * pkin(3) - pkin(9) * t140 + t105 * t90 - t160;
t58 = -t76 * qJD(4) - t114 * t104 - t119 * t92;
t59 = t75 * qJD(4) + t119 * t104 - t114 * t92;
t28 = (-t75 * t99 - t59) * pkin(10) + (t76 * t99 - t58) * pkin(4) + t124;
t155 = t113 * t28 + t118 * t24;
t84 = mrSges(4,1) * t138 + t105 * mrSges(4,2);
t153 = t105 * mrSges(3,1) - mrSges(3,3) * t138 - t84;
t88 = (mrSges(4,2) * t120 - mrSges(4,3) * t115) * t146;
t152 = t88 + (-mrSges(3,1) * t120 + mrSges(3,2) * t115) * t146;
t112 = sin(qJ(6));
t117 = cos(qJ(6));
t136 = -t113 * t24 + t118 * t28;
t63 = -t113 * t76 + t118 * t99;
t41 = t63 * qJD(5) + t113 * t80 + t118 * t59;
t56 = qJDD(5) - t58;
t64 = t113 * t99 + t118 * t76;
t74 = qJD(5) - t75;
t18 = (t63 * t74 - t41) * pkin(11) + (t63 * t64 + t56) * pkin(5) + t136;
t40 = -t64 * qJD(5) - t113 * t59 + t118 * t80;
t53 = t74 * pkin(5) - t64 * pkin(11);
t62 = t63 ^ 2;
t19 = -t62 * pkin(5) + t40 * pkin(11) - t74 * t53 + t155;
t47 = -t112 * t64 + t117 * t63;
t30 = t47 * qJD(6) + t112 * t40 + t117 * t41;
t48 = t112 * t63 + t117 * t64;
t32 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t72 = qJD(6) + t74;
t42 = -t72 * mrSges(7,2) + t47 * mrSges(7,3);
t54 = qJDD(6) + t56;
t16 = m(7) * (-t112 * t19 + t117 * t18) - t30 * mrSges(7,3) + t54 * mrSges(7,1) - t48 * t32 + t72 * t42;
t29 = -t48 * qJD(6) - t112 * t41 + t117 * t40;
t43 = t72 * mrSges(7,1) - t48 * mrSges(7,3);
t17 = m(7) * (t112 * t18 + t117 * t19) + t29 * mrSges(7,3) - t54 * mrSges(7,2) + t47 * t32 - t72 * t43;
t49 = -t63 * mrSges(6,1) + t64 * mrSges(6,2);
t51 = -t74 * mrSges(6,2) + t63 * mrSges(6,3);
t13 = m(6) * t136 + t56 * mrSges(6,1) - t41 * mrSges(6,3) + t112 * t17 + t117 * t16 - t64 * t49 + t74 * t51;
t52 = t74 * mrSges(6,1) - t64 * mrSges(6,3);
t14 = m(6) * t155 - t56 * mrSges(6,2) + t40 * mrSges(6,3) - t112 * t16 + t117 * t17 + t63 * t49 - t74 * t52;
t60 = -t75 * mrSges(5,1) + t76 * mrSges(5,2);
t66 = t99 * mrSges(5,1) - t76 * mrSges(5,3);
t10 = m(5) * t154 - t80 * mrSges(5,2) + t58 * mrSges(5,3) - t113 * t13 + t118 * t14 + t75 * t60 - t99 * t66;
t133 = -t110 * t85 - t158;
t135 = -t114 * t37 + t119 * t39;
t23 = -t80 * pkin(4) - t96 * pkin(10) + t76 * t61 - t135;
t128 = t29 * mrSges(7,1) + t47 * t42 - m(7) * (-t40 * pkin(5) - t62 * pkin(11) + t64 * t53 + t23) - t30 * mrSges(7,2) - t48 * t43;
t123 = m(6) * t23 - t40 * mrSges(6,1) + t41 * mrSges(6,2) - t63 * t51 + t64 * t52 - t128;
t65 = -t99 * mrSges(5,2) + t75 * mrSges(5,3);
t15 = m(5) * t135 + t80 * mrSges(5,1) - t59 * mrSges(5,3) - t76 * t60 + t99 * t65 - t123;
t83 = -mrSges(4,1) * t139 - t105 * mrSges(4,3);
t131 = t119 * t10 - t114 * t15 + m(4) * (-t92 * pkin(2) + (-t105 * t139 - t91) * qJ(3) + t133 + t161) + t83 * t139 - t91 * mrSges(4,3);
t82 = -t105 * mrSges(3,2) + mrSges(3,3) * t139;
t5 = m(3) * t133 + t91 * mrSges(3,2) - t157 * t92 + (t153 * t115 - t120 * t82) * t146 + t131;
t141 = t120 * t151;
t129 = -m(4) * (-t104 * pkin(2) + t130 - t141) - t114 * t10 - t119 * t15;
t6 = m(3) * (t141 - t150) - t156 * t91 + (t82 - t83) * t105 + t157 * t104 - t152 * t138 + t129;
t127 = m(5) * t124 - t58 * mrSges(5,1) + t59 * mrSges(5,2) + t113 * t14 + t118 * t13 - t75 * t65 + t76 * t66;
t125 = -m(4) * t160 + t127;
t8 = m(3) * t132 + t156 * t92 - t153 * t105 + (-mrSges(3,2) + mrSges(4,3)) * t104 + t152 * t139 + t125;
t143 = t111 * t5 + t6 * t147 + t8 * t148;
t2 = m(2) * t134 - t122 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t115 * t6 + t120 * t8;
t1 = m(2) * t137 + qJDD(1) * mrSges(2,1) - t122 * mrSges(2,2) - t110 * t5 + (t115 * t8 + t120 * t6) * t111;
t3 = [-m(1) * g(1) - t116 * t1 + t121 * t2, t2, t8, t92 * mrSges(4,2) - t84 * t138 + t131, t10, t14, t17; -m(1) * g(2) + t121 * t1 + t116 * t2, t1, t6, -t92 * mrSges(4,1) - t104 * mrSges(4,3) - t105 * t84 - t88 * t139 - t125, t15, t13, t16; (-m(1) - m(2)) * g(3) + t143, -m(2) * g(3) + t143, t5, t91 * mrSges(4,1) + t104 * mrSges(4,2) + t105 * t83 + t88 * t138 - t129, t127, t123, -t128;];
f_new  = t3;
