% Calculate vector of cutting forces with Newton-Euler
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-05-05 03:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:57:36
% EndTime: 2019-05-05 02:57:39
% DurationCPUTime: 1.06s
% Computational Cost: add. (8788->180), mult. (17404->222), div. (0->0), fcn. (9579->10), ass. (0->92)
t104 = sin(qJ(3));
t107 = cos(qJ(3));
t130 = qJD(2) * qJD(3);
t124 = t107 * t130;
t110 = qJD(2) ^ 2;
t105 = sin(qJ(2));
t108 = cos(qJ(2));
t99 = sin(pkin(6));
t136 = t108 * t99;
t101 = cos(pkin(6));
t100 = cos(pkin(10));
t98 = sin(pkin(10));
t72 = t98 * g(1) - t100 * g(2);
t138 = t101 * t72;
t73 = -t100 * g(1) - t98 * g(2);
t96 = -g(3) + qJDD(1);
t128 = -t105 * t73 + t108 * t138 + t96 * t136;
t31 = -qJDD(2) * pkin(2) - t110 * pkin(8) - t128;
t67 = t104 * qJDD(2) + t124;
t125 = t104 * t130;
t68 = t107 * qJDD(2) - t125;
t119 = -t68 * pkin(3) + t31 + (-t124 - t67) * qJ(4);
t103 = sin(qJ(6));
t106 = cos(qJ(6));
t131 = t104 * qJD(2);
t135 = t110 * t107 ^ 2;
t148 = 2 * qJD(4);
t74 = -qJD(3) * pkin(4) - qJ(5) * t131;
t112 = -qJ(5) * t135 + qJDD(5) - t119 + (t148 + t74) * t131;
t132 = qJD(2) * t107;
t146 = pkin(4) + pkin(9);
t147 = -pkin(3) - pkin(9);
t16 = (pkin(5) * t107 + t104 * t147) * t130 + t146 * t68 + t112 + t67 * pkin(5);
t109 = qJD(3) ^ 2;
t137 = t105 * t99;
t127 = t105 * t138 + t108 * t73 + t96 * t137;
t32 = -t110 * pkin(2) + qJDD(2) * pkin(8) + t127;
t41 = t101 * t96 - t99 * t72;
t123 = -t104 * t32 + t107 * t41;
t62 = (-pkin(3) * t107 - qJ(4) * t104) * qJD(2);
t122 = t62 * t131 + qJDD(4) - t123;
t133 = t107 * t110;
t129 = qJD(2) * qJD(5);
t153 = -0.2e1 * t104 * t129 + (t124 - t67) * qJ(5);
t66 = (pkin(5) * t104 + pkin(9) * t107) * qJD(2);
t19 = (-pkin(5) - qJ(4)) * t109 + (-pkin(4) * t133 - qJD(2) * t66) * t104 + (-pkin(3) - t146) * qJDD(3) + t122 + t153;
t60 = -t106 * qJD(3) + t103 * t132;
t36 = t60 * qJD(6) - t103 * qJDD(3) - t106 * t68;
t61 = -t103 * qJD(3) - t106 * t132;
t37 = -t60 * mrSges(7,1) + t61 * mrSges(7,2);
t85 = qJD(6) + t131;
t39 = -t85 * mrSges(7,2) + t60 * mrSges(7,3);
t56 = qJDD(6) + t67;
t14 = m(7) * (-t103 * t19 + t106 * t16) - t36 * mrSges(7,3) + t56 * mrSges(7,1) - t61 * t37 + t85 * t39;
t35 = -t61 * qJD(6) - t106 * qJDD(3) + t103 * t68;
t40 = t85 * mrSges(7,1) - t61 * mrSges(7,3);
t15 = m(7) * (t103 * t16 + t106 * t19) + t35 * mrSges(7,3) - t56 * mrSges(7,2) + t60 * t37 - t85 * t40;
t75 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t131;
t78 = -qJD(3) * mrSges(6,1) + mrSges(6,3) * t132;
t121 = t103 * t15 + t106 * t14 + m(6) * (-pkin(3) * t125 + t68 * pkin(4) + t112) + t75 * t131 - t78 * t132 + t67 * mrSges(6,1) - t68 * mrSges(6,2);
t149 = -2 * qJD(4);
t117 = m(5) * ((pkin(3) * qJD(3) + t149) * t131 + t119) - t68 * mrSges(5,1) - t121;
t80 = mrSges(5,2) * t132 + qJD(3) * mrSges(5,3);
t139 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t132 + t80;
t76 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t131;
t77 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t131;
t155 = ((t76 - t77) * t104 - t107 * t139) * qJD(2) + m(4) * t31 - t68 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t67 + t117;
t142 = t104 * t41 + t107 * t32;
t154 = qJDD(3) * qJ(4) + t62 * t132 + t142;
t145 = t109 * pkin(3);
t63 = (-mrSges(5,1) * t107 - mrSges(5,3) * t104) * qJD(2);
t90 = qJD(3) * t148;
t152 = m(5) * (t90 - t145 + t154) + t63 * t132 + qJD(3) * t77 + qJDD(3) * mrSges(5,3);
t143 = mrSges(4,3) + mrSges(5,2);
t65 = (mrSges(6,1) * t104 - mrSges(6,2) * t107) * qJD(2);
t141 = (-mrSges(4,1) * t107 + mrSges(4,2) * t104) * qJD(2) - t65;
t115 = pkin(4) * t135 + t68 * qJ(5) - t154;
t118 = -t35 * mrSges(7,1) - t60 * t39 + m(7) * (qJDD(3) * pkin(5) + qJD(3) * t74 + t90 + t147 * t109 + (-0.2e1 * qJD(5) - t66) * t132 - t115) + t36 * mrSges(7,2) + t61 * t40;
t114 = -m(6) * (0.2e1 * t107 * t129 + t145 + (t149 - t74) * qJD(3) + t115) + t118 - t68 * mrSges(6,3) + qJD(3) * t75 + qJDD(3) * mrSges(6,1);
t10 = m(4) * t142 - qJDD(3) * mrSges(4,2) - qJD(3) * t76 + t132 * t141 + t143 * t68 + t114 + t152;
t26 = -qJDD(3) * pkin(3) - t109 * qJ(4) + t122;
t120 = t103 * t14 - t106 * t15 - m(6) * ((-t104 * t133 - qJDD(3)) * pkin(4) + t26 + t153) + t67 * mrSges(6,3) - qJD(3) * t78 - qJDD(3) * mrSges(6,2);
t116 = m(5) * t26 - t120;
t9 = m(4) * t123 - t143 * t67 + (mrSges(4,1) + mrSges(5,1)) * qJDD(3) + t139 * qJD(3) + (-t63 - t141) * t131 - t116;
t4 = m(3) * t127 - t110 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t107 * t10 - t104 * t9;
t6 = m(3) * t41 + t104 * t10 + t107 * t9;
t8 = m(3) * t128 + qJDD(2) * mrSges(3,1) - t110 * mrSges(3,2) - t155;
t126 = m(2) * t96 + t101 * t6 + t8 * t136 + t4 * t137;
t111 = t132 * t65 - t114;
t2 = m(2) * t73 - t105 * t8 + t108 * t4;
t1 = m(2) * t72 - t99 * t6 + (t105 * t4 + t108 * t8) * t101;
t3 = [-m(1) * g(1) - t98 * t1 + t100 * t2, t2, t4, t10, t68 * mrSges(5,2) - t111 + t152, -t131 * t65 - t120, t15; -m(1) * g(2) + t100 * t1 + t98 * t2, t1, t8, t9, -t67 * mrSges(5,3) + (-t104 * t77 - t107 * t80) * qJD(2) + t117, t111, t14; -m(1) * g(3) + t126, t126, t6, t155, -qJDD(3) * mrSges(5,1) + t67 * mrSges(5,2) - qJD(3) * t80 + (t63 - t65) * t131 + t116, t121, t118;];
f_new  = t3;
