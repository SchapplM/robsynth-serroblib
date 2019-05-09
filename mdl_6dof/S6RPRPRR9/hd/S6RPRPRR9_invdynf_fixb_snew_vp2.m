% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 19:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:29:25
% EndTime: 2019-05-05 19:29:52
% DurationCPUTime: 17.80s
% Computational Cost: add. (279740->218), mult. (932386->322), div. (0->0), fcn. (801710->16), ass. (0->128)
t103 = sin(pkin(13));
t107 = cos(pkin(13));
t113 = sin(qJ(3));
t109 = cos(pkin(7));
t117 = cos(qJ(3));
t147 = t109 * t117;
t105 = sin(pkin(7));
t151 = t105 * t117;
t104 = sin(pkin(12));
t106 = sin(pkin(6));
t110 = cos(pkin(6));
t108 = cos(pkin(12));
t159 = pkin(9) * t104;
t130 = -pkin(2) * t108 - t105 * t159;
t146 = qJD(1) * t106;
t156 = pkin(9) * qJDD(1);
t127 = qJD(1) * t130 * t146 + t109 * t156;
t141 = qJD(2) * t146;
t150 = t106 * t108;
t119 = qJD(1) ^ 2;
t114 = sin(qJ(1));
t118 = cos(qJ(1));
t140 = t114 * g(1) - t118 * g(2);
t155 = qJ(2) * t106;
t95 = qJDD(1) * pkin(1) + t119 * t155 + t140;
t157 = t110 * t95;
t131 = -g(3) * t150 - 0.2e1 * t104 * t141 + t108 * t157;
t149 = t106 * t109;
t92 = (t105 * t110 + t108 * t149) * qJD(1) * pkin(9);
t135 = -t118 * g(1) - t114 * g(2);
t96 = -t119 * pkin(1) + qJDD(1) * t155 + t135;
t59 = (pkin(2) * qJDD(1) + qJD(1) * t92) * t110 + (-t127 * t106 - t96) * t104 + t131;
t142 = t104 * t157 + (0.2e1 * t141 + t96) * t108;
t97 = (pkin(2) * t110 - t149 * t159) * qJD(1);
t60 = (-qJD(1) * t97 + t105 * t156) * t110 + (-g(3) * t104 + t127 * t108) * t106 + t142;
t137 = -t110 * g(3) + qJDD(2);
t67 = (-t95 + t130 * qJDD(1) + (t104 * t97 - t108 * t92) * qJD(1)) * t106 + t137;
t136 = -t113 * t60 + t59 * t147 + t67 * t151;
t148 = t109 * t113;
t152 = t105 * t113;
t122 = t110 * t152 + (t104 * t117 + t108 * t148) * t106;
t121 = t110 * t151 + (-t104 * t113 + t108 * t147) * t106;
t84 = t121 * qJD(1);
t79 = t84 * qJD(3) + t122 * qJDD(1);
t85 = t122 * qJD(1);
t126 = -t105 * t150 + t109 * t110;
t90 = t126 * qJDD(1) + qJDD(3);
t93 = t126 * qJD(1) + qJD(3);
t31 = (t84 * t93 - t79) * qJ(4) + (t84 * t85 + t90) * pkin(3) + t136;
t143 = t117 * t60 + t59 * t148 + t67 * t152;
t78 = -t85 * qJD(3) + t121 * qJDD(1);
t81 = t93 * pkin(3) - t85 * qJ(4);
t83 = t84 ^ 2;
t34 = -pkin(3) * t83 + qJ(4) * t78 - t81 * t93 + t143;
t76 = t103 * t84 + t107 * t85;
t160 = -0.2e1 * qJD(4) * t76 - t103 * t34 + t107 * t31;
t112 = sin(qJ(5));
t116 = cos(qJ(5));
t75 = -t103 * t85 + t107 * t84;
t144 = 0.2e1 * qJD(4) * t75 + t103 * t31 + t107 * t34;
t51 = -pkin(4) * t75 - pkin(10) * t76;
t89 = t93 ^ 2;
t25 = -pkin(4) * t89 + pkin(10) * t90 + t51 * t75 + t144;
t139 = -t105 * t59 + t109 * t67;
t124 = -t78 * pkin(3) - t83 * qJ(4) + t85 * t81 + qJDD(4) + t139;
t54 = -t103 * t79 + t107 * t78;
t55 = t103 * t78 + t107 * t79;
t27 = (-t75 * t93 - t55) * pkin(10) + (t76 * t93 - t54) * pkin(4) + t124;
t158 = t112 * t27 + t116 * t25;
t153 = t104 * t106;
t111 = sin(qJ(6));
t115 = cos(qJ(6));
t65 = -t112 * t76 + t116 * t93;
t66 = t112 * t93 + t116 * t76;
t44 = -pkin(5) * t65 - pkin(11) * t66;
t53 = qJDD(5) - t54;
t74 = qJD(5) - t75;
t73 = t74 ^ 2;
t21 = -pkin(5) * t73 + pkin(11) * t53 + t44 * t65 + t158;
t24 = -t90 * pkin(4) - t89 * pkin(10) + t76 * t51 - t160;
t40 = -qJD(5) * t66 - t112 * t55 + t116 * t90;
t41 = qJD(5) * t65 + t112 * t90 + t116 * t55;
t22 = (-t65 * t74 - t41) * pkin(11) + (t66 * t74 - t40) * pkin(5) + t24;
t45 = -t111 * t66 + t115 * t74;
t29 = qJD(6) * t45 + t111 * t53 + t115 * t41;
t46 = t111 * t74 + t115 * t66;
t36 = -mrSges(7,1) * t45 + mrSges(7,2) * t46;
t64 = qJD(6) - t65;
t37 = -mrSges(7,2) * t64 + mrSges(7,3) * t45;
t39 = qJDD(6) - t40;
t18 = m(7) * (-t111 * t21 + t115 * t22) - t29 * mrSges(7,3) + t39 * mrSges(7,1) - t46 * t36 + t64 * t37;
t28 = -qJD(6) * t46 - t111 * t41 + t115 * t53;
t38 = mrSges(7,1) * t64 - mrSges(7,3) * t46;
t19 = m(7) * (t111 * t22 + t115 * t21) + t28 * mrSges(7,3) - t39 * mrSges(7,2) + t45 * t36 - t64 * t38;
t43 = -mrSges(6,1) * t65 + mrSges(6,2) * t66;
t48 = mrSges(6,1) * t74 - mrSges(6,3) * t66;
t15 = m(6) * t158 - t53 * mrSges(6,2) + t40 * mrSges(6,3) - t111 * t18 + t115 * t19 + t65 * t43 - t74 * t48;
t132 = -t112 * t25 + t116 * t27;
t123 = m(7) * (-pkin(5) * t53 - pkin(11) * t73 + t66 * t44 - t132) - t28 * mrSges(7,1) + t29 * mrSges(7,2) - t45 * t37 + t46 * t38;
t47 = -mrSges(6,2) * t74 + mrSges(6,3) * t65;
t17 = m(6) * t132 + t53 * mrSges(6,1) - t41 * mrSges(6,3) - t66 * t43 + t74 * t47 - t123;
t68 = -t93 * mrSges(5,2) + t75 * mrSges(5,3);
t69 = t93 * mrSges(5,1) - t76 * mrSges(5,3);
t125 = m(5) * t124 - t54 * mrSges(5,1) + t55 * mrSges(5,2) + t112 * t15 + t116 * t17 - t75 * t68 + t76 * t69;
t80 = -t93 * mrSges(4,2) + t84 * mrSges(4,3);
t82 = t93 * mrSges(4,1) - t85 * mrSges(4,3);
t12 = m(4) * t139 - t78 * mrSges(4,1) + t79 * mrSges(4,2) - t84 * t80 + t85 * t82 + t125;
t129 = t110 * mrSges(3,1) - mrSges(3,3) * t153;
t50 = -mrSges(5,1) * t75 + mrSges(5,2) * t76;
t11 = m(5) * t144 - t90 * mrSges(5,2) + t54 * mrSges(5,3) - t112 * t17 + t116 * t15 + t75 * t50 - t93 * t69;
t120 = m(6) * t24 - t40 * mrSges(6,1) + t41 * mrSges(6,2) + t111 * t19 + t115 * t18 - t65 * t47 + t66 * t48;
t13 = m(5) * t160 + t90 * mrSges(5,1) - t55 * mrSges(5,3) - t76 * t50 + t93 * t68 - t120;
t77 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t10 = m(4) * t143 - t90 * mrSges(4,2) + t78 * mrSges(4,3) - t103 * t13 + t107 * t11 + t84 * t77 - t93 * t82;
t9 = m(4) * t136 + t90 * mrSges(4,1) - t79 * mrSges(4,3) + t103 * t11 + t107 * t13 - t85 * t77 + t93 * t80;
t134 = t10 * t113 + t117 * t9;
t133 = -mrSges(3,1) * t108 + mrSges(3,2) * t104;
t94 = t133 * t146;
t128 = -t110 * mrSges(3,2) + mrSges(3,3) * t150;
t99 = t128 * qJD(1);
t4 = m(3) * (-t104 * t96 + t131) - t105 * t12 + t134 * t109 + t129 * qJDD(1) + (t110 * t99 - t94 * t153) * qJD(1);
t98 = t129 * qJD(1);
t6 = m(3) * t137 + t109 * t12 + t134 * t105 + (-m(3) * t95 + t133 * qJDD(1) + (t104 * t98 - t108 * t99) * qJD(1)) * t106;
t8 = m(3) * (-g(3) * t153 + t142) + t117 * t10 - t113 * t9 + t128 * qJDD(1) + (-t110 * t98 + t94 * t150) * qJD(1);
t145 = t110 * t6 + t4 * t150 + t8 * t153;
t2 = m(2) * t135 - t119 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t104 * t4 + t108 * t8;
t1 = m(2) * t140 + qJDD(1) * mrSges(2,1) - t119 * mrSges(2,2) - t106 * t6 + (t104 * t8 + t108 * t4) * t110;
t3 = [-m(1) * g(1) - t1 * t114 + t118 * t2, t2, t8, t10, t11, t15, t19; -m(1) * g(2) + t1 * t118 + t114 * t2, t1, t4, t9, t13, t17, t18; (-m(1) - m(2)) * g(3) + t145, -m(2) * g(3) + t145, t6, t12, t125, t120, t123;];
f_new  = t3;
