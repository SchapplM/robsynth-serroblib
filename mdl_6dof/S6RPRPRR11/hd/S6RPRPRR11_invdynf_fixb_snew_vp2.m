% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR11
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
% Datum: 2019-05-05 20:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR11_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR11_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:11:08
% EndTime: 2019-05-05 20:11:34
% DurationCPUTime: 18.18s
% Computational Cost: add. (290508->218), mult. (922290->318), div. (0->0), fcn. (791628->16), ass. (0->127)
t109 = sin(pkin(7));
t112 = cos(pkin(12));
t114 = cos(pkin(6));
t110 = sin(pkin(6));
t113 = cos(pkin(7));
t155 = t110 * t113;
t167 = t109 * t114 + t112 * t155;
t117 = sin(qJ(3));
t121 = cos(qJ(3));
t108 = sin(pkin(12));
t159 = t108 * t110;
t166 = t117 * t159 - t121 * t167;
t164 = pkin(9) * t108;
t133 = -pkin(2) * t112 - t109 * t164;
t153 = qJD(1) * t110;
t161 = pkin(9) * qJDD(1);
t130 = qJD(1) * t133 * t153 + t113 * t161;
t145 = qJD(2) * t153;
t156 = t110 * t112;
t123 = qJD(1) ^ 2;
t118 = sin(qJ(1));
t122 = cos(qJ(1));
t144 = t118 * g(1) - g(2) * t122;
t160 = qJ(2) * t110;
t95 = qJDD(1) * pkin(1) + t123 * t160 + t144;
t162 = t114 * t95;
t134 = -g(3) * t156 - 0.2e1 * t108 * t145 + t112 * t162;
t92 = t167 * qJD(1) * pkin(9);
t139 = -g(1) * t122 - g(2) * t118;
t96 = -pkin(1) * t123 + qJDD(1) * t160 + t139;
t55 = (pkin(2) * qJDD(1) + qJD(1) * t92) * t114 + (-t130 * t110 - t96) * t108 + t134;
t149 = t108 * t162 + (0.2e1 * t145 + t96) * t112;
t97 = (pkin(2) * t114 - t155 * t164) * qJD(1);
t56 = (-qJD(1) * t97 + t109 * t161) * t114 + (-g(3) * t108 + t130 * t112) * t110 + t149;
t142 = -t114 * g(3) + qJDD(2);
t64 = (-t95 + t133 * qJDD(1) + (t108 * t97 - t112 * t92) * qJD(1)) * t110 + t142;
t165 = -t117 * t56 + (t109 * t64 + t113 * t55) * t121;
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t107 = sin(pkin(13));
t111 = cos(pkin(13));
t154 = t113 * t117;
t157 = t109 * t117;
t150 = t121 * t56 + t55 * t154 + t64 * t157;
t84 = t166 * qJD(1);
t126 = t114 * t157 + (t108 * t121 + t112 * t154) * t110;
t85 = t126 * qJD(1);
t71 = pkin(3) * t84 - qJ(4) * t85;
t129 = -t109 * t156 + t113 * t114;
t93 = t129 * qJD(1) + qJD(3);
t89 = t93 ^ 2;
t90 = t129 * qJDD(1) + qJDD(3);
t33 = -pkin(3) * t89 + qJ(4) * t90 - t71 * t84 + t150;
t143 = -t109 * t55 + t113 * t64;
t73 = qJD(3) * t85 + t166 * qJDD(1);
t74 = -t84 * qJD(3) + t126 * qJDD(1);
t36 = (t84 * t93 - t74) * qJ(4) + (t85 * t93 + t73) * pkin(3) + t143;
t79 = t107 * t93 + t111 * t85;
t140 = -0.2e1 * qJD(4) * t79 - t107 * t33 + t111 * t36;
t69 = t107 * t90 + t111 * t74;
t78 = -t107 * t85 + t111 * t93;
t24 = (t78 * t84 - t69) * pkin(10) + (t78 * t79 + t73) * pkin(4) + t140;
t151 = 0.2e1 * qJD(4) * t78 + t107 * t36 + t111 * t33;
t67 = pkin(4) * t84 - pkin(10) * t79;
t68 = -t107 * t74 + t111 * t90;
t77 = t78 ^ 2;
t26 = -pkin(4) * t77 + pkin(10) * t68 - t67 * t84 + t151;
t163 = t116 * t24 + t120 * t26;
t115 = sin(qJ(6));
t119 = cos(qJ(6));
t58 = -t116 * t79 + t120 * t78;
t59 = t116 * t78 + t120 * t79;
t46 = -pkin(5) * t58 - pkin(11) * t59;
t70 = qJDD(5) + t73;
t83 = qJD(5) + t84;
t82 = t83 ^ 2;
t21 = -pkin(5) * t82 + pkin(11) * t70 + t58 * t46 + t163;
t32 = -t90 * pkin(3) - t89 * qJ(4) + t85 * t71 + qJDD(4) - t165;
t125 = -t68 * pkin(4) - t77 * pkin(10) + t79 * t67 + t32;
t40 = -t59 * qJD(5) - t116 * t69 + t120 * t68;
t41 = t58 * qJD(5) + t116 * t68 + t120 * t69;
t22 = (-t58 * t83 - t41) * pkin(11) + (t59 * t83 - t40) * pkin(5) + t125;
t47 = -t115 * t59 + t119 * t83;
t30 = t47 * qJD(6) + t115 * t70 + t119 * t41;
t48 = t115 * t83 + t119 * t59;
t37 = -mrSges(7,1) * t47 + mrSges(7,2) * t48;
t39 = qJDD(6) - t40;
t57 = qJD(6) - t58;
t42 = -mrSges(7,2) * t57 + mrSges(7,3) * t47;
t18 = m(7) * (-t115 * t21 + t119 * t22) - t30 * mrSges(7,3) + t39 * mrSges(7,1) - t48 * t37 + t57 * t42;
t29 = -t48 * qJD(6) - t115 * t41 + t119 * t70;
t43 = mrSges(7,1) * t57 - mrSges(7,3) * t48;
t19 = m(7) * (t115 * t22 + t119 * t21) + t29 * mrSges(7,3) - t39 * mrSges(7,2) + t47 * t37 - t57 * t43;
t45 = -mrSges(6,1) * t58 + mrSges(6,2) * t59;
t50 = mrSges(6,1) * t83 - t59 * mrSges(6,3);
t13 = m(6) * t163 - t70 * mrSges(6,2) + t40 * mrSges(6,3) - t115 * t18 + t119 * t19 + t58 * t45 - t83 * t50;
t135 = -t116 * t26 + t120 * t24;
t127 = m(7) * (-pkin(5) * t70 - pkin(11) * t82 + t59 * t46 - t135) - t29 * mrSges(7,1) + t30 * mrSges(7,2) - t47 * t42 + t48 * t43;
t49 = -mrSges(6,2) * t83 + t58 * mrSges(6,3);
t15 = m(6) * t135 + t70 * mrSges(6,1) - t41 * mrSges(6,3) - t59 * t45 + t83 * t49 - t127;
t60 = -mrSges(5,1) * t78 + mrSges(5,2) * t79;
t65 = -mrSges(5,2) * t84 + mrSges(5,3) * t78;
t11 = m(5) * t140 + t73 * mrSges(5,1) - t69 * mrSges(5,3) + t116 * t13 + t120 * t15 - t79 * t60 + t84 * t65;
t66 = mrSges(5,1) * t84 - mrSges(5,3) * t79;
t12 = m(5) * t151 - t73 * mrSges(5,2) + t68 * mrSges(5,3) - t116 * t15 + t120 * t13 + t78 * t60 - t84 * t66;
t80 = -mrSges(4,2) * t93 - mrSges(4,3) * t84;
t81 = mrSges(4,1) * t93 - mrSges(4,3) * t85;
t10 = m(4) * t143 + t73 * mrSges(4,1) + t74 * mrSges(4,2) + t107 * t12 + t111 * t11 + t84 * t80 + t85 * t81;
t132 = mrSges(3,1) * t114 - mrSges(3,3) * t159;
t128 = -m(6) * t125 + t40 * mrSges(6,1) - t41 * mrSges(6,2) - t115 * t19 - t119 * t18 + t58 * t49 - t59 * t50;
t124 = m(5) * t32 - t68 * mrSges(5,1) + t69 * mrSges(5,2) - t78 * t65 + t79 * t66 - t128;
t72 = mrSges(4,1) * t84 + mrSges(4,2) * t85;
t14 = m(4) * t165 + t90 * mrSges(4,1) - t74 * mrSges(4,3) - t85 * t72 + t93 * t80 - t124;
t9 = m(4) * t150 - t90 * mrSges(4,2) - t73 * mrSges(4,3) - t107 * t11 + t111 * t12 - t84 * t72 - t93 * t81;
t138 = t117 * t9 + t121 * t14;
t137 = -mrSges(3,1) * t112 + mrSges(3,2) * t108;
t94 = t137 * t153;
t131 = -mrSges(3,2) * t114 + mrSges(3,3) * t156;
t99 = t131 * qJD(1);
t4 = m(3) * (-t108 * t96 + t134) - t109 * t10 + t138 * t113 + t132 * qJDD(1) + (t114 * t99 - t94 * t159) * qJD(1);
t98 = t132 * qJD(1);
t6 = m(3) * t142 + t113 * t10 + t138 * t109 + (-m(3) * t95 + t137 * qJDD(1) + (t108 * t98 - t112 * t99) * qJD(1)) * t110;
t8 = m(3) * (-g(3) * t159 + t149) + t121 * t9 - t117 * t14 + t131 * qJDD(1) + (-t114 * t98 + t94 * t156) * qJD(1);
t152 = t114 * t6 + t4 * t156 + t8 * t159;
t2 = m(2) * t139 - t123 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t108 * t4 + t112 * t8;
t1 = m(2) * t144 + qJDD(1) * mrSges(2,1) - t123 * mrSges(2,2) - t110 * t6 + (t108 * t8 + t112 * t4) * t114;
t3 = [-m(1) * g(1) - t1 * t118 + t122 * t2, t2, t8, t9, t12, t13, t19; -m(1) * g(2) + t1 * t122 + t118 * t2, t1, t4, t14, t11, t15, t18; (-m(1) - m(2)) * g(3) + t152, -m(2) * g(3) + t152, t6, t10, t124, -t128, t127;];
f_new  = t3;
