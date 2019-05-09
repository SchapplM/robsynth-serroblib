% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-05-05 18:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:52:04
% EndTime: 2019-05-05 18:52:09
% DurationCPUTime: 2.42s
% Computational Cost: add. (28932->190), mult. (69940->240), div. (0->0), fcn. (51560->10), ass. (0->99)
t115 = qJD(1) ^ 2;
t105 = sin(pkin(10));
t98 = t105 ^ 2;
t106 = cos(pkin(10));
t99 = t106 ^ 2;
t152 = (t98 + t99) * mrSges(3,3);
t109 = sin(qJ(3));
t148 = cos(qJ(3));
t151 = t105 * t109 - t106 * t148;
t150 = 2 * qJD(4);
t130 = -mrSges(3,1) * t106 + mrSges(3,2) * t105;
t128 = mrSges(3,3) * qJDD(1) + t115 * t130;
t137 = qJD(1) * qJD(2);
t134 = -t106 * g(3) - 0.2e1 * t105 * t137;
t100 = -qJD(3) + qJD(5);
t107 = sin(qJ(6));
t111 = cos(qJ(6));
t108 = sin(qJ(5));
t112 = cos(qJ(5));
t82 = t151 * qJD(1);
t139 = t82 * qJD(3);
t114 = qJD(3) ^ 2;
t110 = sin(qJ(1));
t113 = cos(qJ(1));
t131 = -g(1) * t113 - g(2) * t110;
t84 = -pkin(1) * t115 + qJDD(1) * qJ(2) + t131;
t52 = (pkin(2) * t106 * t115 - pkin(7) * qJDD(1) - t84) * t105 + t134;
t133 = -g(3) * t105 + (0.2e1 * t137 + t84) * t106;
t136 = qJDD(1) * t106;
t141 = t115 * t99;
t53 = -pkin(2) * t141 + pkin(7) * t136 + t133;
t132 = -t109 * t53 + t148 * t52;
t125 = t148 * t105 + t106 * t109;
t83 = t125 * qJD(1);
t62 = pkin(3) * t82 - qJ(4) * t83;
t30 = -qJDD(3) * pkin(3) - t114 * qJ(4) + t83 * t62 + qJDD(4) - t132;
t70 = t125 * qJDD(1) - t139;
t22 = (-t70 - t139) * pkin(8) + (t82 * t83 - qJDD(3)) * pkin(4) + t30;
t145 = t109 * t52 + t148 * t53;
t124 = -pkin(3) * t114 + qJDD(3) * qJ(4) + qJD(3) * t150 - t82 * t62 + t145;
t138 = t83 * qJD(3);
t69 = qJDD(1) * t151 + t138;
t76 = -qJD(3) * pkin(4) - pkin(8) * t83;
t81 = t82 ^ 2;
t24 = -pkin(4) * t81 + pkin(8) * t69 + qJD(3) * t76 + t124;
t146 = t108 * t22 + t112 * t24;
t55 = -t108 * t83 + t112 * t82;
t56 = t108 * t82 + t112 * t83;
t42 = -pkin(5) * t55 - pkin(9) * t56;
t96 = t100 ^ 2;
t97 = -qJDD(3) + qJDD(5);
t17 = -pkin(5) * t96 + pkin(9) * t97 + t42 * t55 + t146;
t143 = t110 * g(1) - t113 * g(2);
t80 = -qJDD(1) * pkin(1) - t115 * qJ(2) + qJDD(2) - t143;
t122 = -pkin(2) * t136 + t80 + (-t115 * t98 - t141) * pkin(7);
t121 = t69 * pkin(3) + t122 + (t139 - t70) * qJ(4);
t118 = -pkin(3) * t138 - t69 * pkin(4) - t81 * pkin(8) - t121 + (t150 + t76) * t83;
t34 = -qJD(5) * t56 - t108 * t70 + t112 * t69;
t35 = qJD(5) * t55 + t108 * t69 + t112 * t70;
t18 = t118 + (t100 * t56 - t34) * pkin(5) + (-t100 * t55 - t35) * pkin(9);
t45 = t100 * t111 - t107 * t56;
t26 = qJD(6) * t45 + t107 * t97 + t111 * t35;
t33 = qJDD(6) - t34;
t46 = t100 * t107 + t111 * t56;
t36 = -mrSges(7,1) * t45 + mrSges(7,2) * t46;
t54 = qJD(6) - t55;
t37 = -mrSges(7,2) * t54 + mrSges(7,3) * t45;
t14 = m(7) * (-t107 * t17 + t111 * t18) - t26 * mrSges(7,3) + t33 * mrSges(7,1) - t46 * t36 + t54 * t37;
t25 = -qJD(6) * t46 - t107 * t35 + t111 * t97;
t38 = mrSges(7,1) * t54 - mrSges(7,3) * t46;
t15 = m(7) * (t107 * t18 + t111 * t17) + t25 * mrSges(7,3) - t33 * mrSges(7,2) + t45 * t36 - t54 * t38;
t41 = -mrSges(6,1) * t55 + mrSges(6,2) * t56;
t49 = mrSges(6,1) * t100 - mrSges(6,3) * t56;
t10 = m(6) * t146 - t97 * mrSges(6,2) + t34 * mrSges(6,3) - t100 * t49 - t107 * t14 + t111 * t15 + t55 * t41;
t129 = -t108 * t24 + t112 * t22;
t120 = m(7) * (-pkin(5) * t97 - pkin(9) * t96 + t42 * t56 - t129) - t25 * mrSges(7,1) + t26 * mrSges(7,2) - t45 * t37 + t46 * t38;
t48 = -mrSges(6,2) * t100 + mrSges(6,3) * t55;
t11 = m(6) * t129 + t97 * mrSges(6,1) - t35 * mrSges(6,3) + t100 * t48 - t56 * t41 - t120;
t74 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t83;
t127 = m(5) * t124 + qJDD(3) * mrSges(5,3) + qJD(3) * t74 + t112 * t10 - t108 * t11;
t63 = mrSges(5,1) * t82 - mrSges(5,3) * t83;
t144 = -mrSges(4,1) * t82 - mrSges(4,2) * t83 - t63;
t147 = -mrSges(4,3) - mrSges(5,2);
t73 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t83;
t6 = m(4) * t145 - qJDD(3) * mrSges(4,2) - qJD(3) * t73 + t144 * t82 + t147 * t69 + t127;
t123 = -m(5) * t30 - t108 * t10 - t112 * t11;
t72 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t82;
t75 = -mrSges(5,2) * t82 + qJD(3) * mrSges(5,3);
t7 = m(4) * t132 + t144 * t83 + t147 * t70 + (mrSges(4,1) + mrSges(5,1)) * qJDD(3) + (t72 + t75) * qJD(3) + t123;
t4 = m(3) * t134 + t109 * t6 + t148 * t7 + (-m(3) * t84 - t128) * t105;
t5 = m(3) * t133 + t128 * t106 - t109 * t7 + t148 * t6;
t149 = t105 * t5 + t106 * t4;
t126 = m(6) * t118 - t34 * mrSges(6,1) + t35 * mrSges(6,2) + t107 * t15 + t111 * t14 - t55 * t48 + t56 * t49;
t119 = t70 * mrSges(5,3) + t83 * t74 + t126 - m(5) * ((pkin(3) * qJD(3) - (2 * qJD(4))) * t83 + t121) - t82 * t75 - t69 * mrSges(5,1);
t117 = m(4) * t122 + t69 * mrSges(4,1) + t70 * mrSges(4,2) + t82 * t72 + t83 * t73 - t119;
t116 = m(3) * t80 + t117;
t8 = (mrSges(2,1) - t130) * qJDD(1) + m(2) * t143 - t116 + (-mrSges(2,2) + t152) * t115;
t1 = m(2) * t131 - t115 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t105 * t4 + t106 * t5;
t2 = [-m(1) * g(1) + t1 * t113 - t110 * t8, t1, t5, t6, -t69 * mrSges(5,2) - t82 * t63 + t127, t10, t15; -m(1) * g(2) + t1 * t110 + t113 * t8, t8, t4, t7, -t119, t11, t14; (-m(1) - m(2)) * g(3) + t149, -m(2) * g(3) + t149, t130 * qJDD(1) - t115 * t152 + t116, t117, -qJDD(3) * mrSges(5,1) + t70 * mrSges(5,2) - qJD(3) * t75 + t83 * t63 - t123, t126, t120;];
f_new  = t2;
