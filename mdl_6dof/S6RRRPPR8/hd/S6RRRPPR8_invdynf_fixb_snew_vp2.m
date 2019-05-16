% Calculate vector of cutting forces with Newton-Euler
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-05-07 06:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPPR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:09:23
% EndTime: 2019-05-07 06:09:28
% DurationCPUTime: 2.23s
% Computational Cost: add. (24968->216), mult. (54020->266), div. (0->0), fcn. (40150->10), ass. (0->107)
t114 = sin(pkin(6));
t121 = cos(qJ(2));
t146 = qJD(1) * t121;
t141 = t114 * t146;
t104 = -qJD(3) + t141;
t115 = cos(pkin(6));
t111 = t115 * qJD(1) + qJD(2);
t117 = sin(qJ(3));
t118 = sin(qJ(2));
t147 = qJD(1) * t114;
t142 = t118 * t147;
t162 = cos(qJ(3));
t85 = t117 * t111 + t162 * t142;
t71 = t104 * pkin(4) - t85 * qJ(5);
t170 = (2 * qJD(4)) + t71;
t109 = t111 ^ 2;
t110 = t115 * qJDD(1) + qJDD(2);
t123 = qJD(1) ^ 2;
t119 = sin(qJ(1));
t122 = cos(qJ(1));
t140 = t119 * g(1) - t122 * g(2);
t93 = t123 * t114 * pkin(8) + qJDD(1) * pkin(1) + t140;
t150 = t115 * t93;
t138 = -t122 * g(1) - t119 * g(2);
t144 = qJDD(1) * t114;
t94 = -t123 * pkin(1) + pkin(8) * t144 + t138;
t152 = t118 * t150 + t121 * t94;
t96 = (-pkin(2) * t121 - pkin(9) * t118) * t147;
t37 = -t109 * pkin(2) + t110 * pkin(9) + (-g(3) * t118 + t96 * t146) * t114 + t152;
t160 = t115 * g(3);
t97 = (qJD(2) * t146 + qJDD(1) * t118) * t114;
t98 = -qJD(2) * t142 + t121 * t144;
t38 = -t98 * pkin(2) - t97 * pkin(9) - t160 + (-t93 + (pkin(2) * t118 - pkin(9) * t121) * t111 * qJD(1)) * t114;
t156 = t117 * t38 + t162 * t37;
t167 = t104 ^ 2;
t84 = -t162 * t111 + t117 * t142;
t61 = t84 * pkin(3) - t85 * qJ(4);
t90 = -qJDD(3) + t98;
t137 = pkin(3) * t167 + t90 * qJ(4) + t84 * t61 - t156;
t73 = t104 * mrSges(5,1) + t85 * mrSges(5,2);
t166 = -2 * qJD(4);
t99 = t104 * t166;
t169 = m(5) * (-t137 + t99) - t104 * t73 - t90 * mrSges(5,3);
t151 = t104 * t84;
t148 = t114 * t121;
t143 = -g(3) * t148 - t118 * t94 + t121 * t150;
t36 = -t110 * pkin(2) - t109 * pkin(9) + t96 * t142 - t143;
t58 = t85 * qJD(3) - t162 * t110 + t117 * t97;
t59 = -t84 * qJD(3) + t117 * t110 + t162 * t97;
t132 = t58 * pkin(3) + t36 + (-t151 - t59) * qJ(4);
t116 = sin(qJ(6));
t120 = cos(qJ(6));
t83 = t84 ^ 2;
t127 = -t83 * qJ(5) + t170 * t85 + qJDD(5) - t132;
t163 = -pkin(4) - pkin(10);
t16 = (pkin(5) * t84 + (pkin(3) + pkin(10)) * t85) * t104 + t163 * t58 + t127 + t59 * pkin(5);
t139 = -t117 * t37 + t162 * t38;
t130 = -qJ(4) * t167 + t85 * t61 + qJDD(4) - t139;
t164 = -2 * qJD(5);
t124 = (-t59 + t151) * qJ(5) + t130 + (t84 * pkin(4) + t164) * t85;
t64 = t85 * pkin(5) - t84 * pkin(10);
t17 = -t167 * pkin(5) - t85 * t64 + (pkin(3) - t163) * t90 + t124;
t67 = t120 * t104 - t116 * t84;
t31 = t67 * qJD(6) + t116 * t90 + t120 * t58;
t68 = t116 * t104 + t120 * t84;
t39 = -t67 * mrSges(7,1) + t68 * mrSges(7,2);
t82 = qJD(6) + t85;
t44 = -t82 * mrSges(7,2) + t67 * mrSges(7,3);
t52 = qJDD(6) + t59;
t14 = m(7) * (-t116 * t17 + t120 * t16) - t31 * mrSges(7,3) + t52 * mrSges(7,1) - t68 * t39 + t82 * t44;
t30 = -t68 * qJD(6) - t116 * t58 + t120 * t90;
t45 = t82 * mrSges(7,1) - t68 * mrSges(7,3);
t15 = m(7) * (t116 * t16 + t120 * t17) + t30 * mrSges(7,3) - t52 * mrSges(7,2) + t67 * t39 - t82 * t45;
t161 = pkin(3) * t104;
t69 = t104 * mrSges(6,1) - t84 * mrSges(6,3);
t74 = -t104 * mrSges(6,2) - t85 * mrSges(6,3);
t135 = t116 * t15 + t120 * t14 + m(6) * (-t58 * pkin(4) + t85 * t161 + t127) + t58 * mrSges(6,2) + t59 * mrSges(6,1) + t84 * t69 + t85 * t74;
t75 = -t84 * mrSges(5,2) - t104 * mrSges(5,3);
t131 = m(5) * ((t166 - t161) * t85 + t132) + t58 * mrSges(5,1) + t84 * t75 - t135;
t70 = t104 * mrSges(4,2) - t84 * mrSges(4,3);
t72 = -t104 * mrSges(4,1) - t85 * mrSges(4,3);
t168 = m(4) * t36 + t58 * mrSges(4,1) + (t72 - t73) * t85 + (mrSges(4,2) - mrSges(5,3)) * t59 + t84 * t70 + t131;
t159 = mrSges(5,1) - mrSges(6,2);
t157 = -mrSges(4,3) - mrSges(5,2);
t62 = t84 * mrSges(5,1) - t85 * mrSges(5,3);
t155 = -t84 * mrSges(4,1) - t85 * mrSges(4,2) - t62;
t154 = -t69 + t75;
t149 = t114 * t118;
t129 = t83 * pkin(4) - t58 * qJ(5) + t137;
t134 = -t30 * mrSges(7,1) - t67 * t44 + m(7) * (-t90 * pkin(5) - t167 * pkin(10) - t104 * t71 + t99 + ((2 * qJD(5)) + t64) * t84 - t129) + t31 * mrSges(7,2) + t68 * t45;
t60 = t85 * mrSges(6,1) + t84 * mrSges(6,2);
t128 = -m(6) * (t170 * t104 + t84 * t164 + t129) + t134 + t84 * t60 + t58 * mrSges(6,3);
t10 = t128 + (mrSges(4,2) - mrSges(6,1)) * t90 + t157 * t58 + t155 * t84 + (t72 - t74) * t104 + m(4) * t156 + t169;
t136 = t116 * t14 - t120 * t15 - m(6) * ((pkin(3) + pkin(4)) * t90 + t124) + t85 * t60 + t59 * mrSges(6,3);
t133 = m(5) * (t90 * pkin(3) + t130) - t136;
t9 = m(4) * t139 + t155 * t85 + t157 * t59 + (-mrSges(4,1) - t159) * t90 + (-t70 - t154) * t104 - t133;
t91 = t111 * mrSges(3,1) - mrSges(3,3) * t142;
t95 = (-mrSges(3,1) * t121 + mrSges(3,2) * t118) * t147;
t4 = m(3) * (-g(3) * t149 + t152) + t98 * mrSges(3,3) - t110 * mrSges(3,2) + t95 * t141 - t111 * t91 + t162 * t10 - t117 * t9;
t92 = -t111 * mrSges(3,2) + mrSges(3,3) * t141;
t6 = m(3) * (-t114 * t93 - t160) + t97 * mrSges(3,2) - t98 * mrSges(3,1) + t117 * t10 + t162 * t9 + (t118 * t91 - t121 * t92) * t147;
t8 = m(3) * t143 + t110 * mrSges(3,1) - t97 * mrSges(3,3) + t111 * t92 - t95 * t142 - t168;
t145 = t115 * t6 + t8 * t148 + t4 * t149;
t125 = t90 * mrSges(6,1) + t104 * t74 - t128;
t2 = m(2) * t138 - t123 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t118 * t8 + t121 * t4;
t1 = m(2) * t140 + qJDD(1) * mrSges(2,1) - t123 * mrSges(2,2) - t114 * t6 + (t118 * t4 + t121 * t8) * t115;
t3 = [-m(1) * g(1) - t119 * t1 + t122 * t2, t2, t4, t10, -t58 * mrSges(5,2) - t84 * t62 - t125 + t169, -t90 * mrSges(6,2) - t104 * t69 - t136, t15; -m(1) * g(2) + t122 * t1 + t119 * t2, t1, t8, t9, -t59 * mrSges(5,3) - t85 * t73 + t131, t125, t14; (-m(1) - m(2)) * g(3) + t145, -m(2) * g(3) + t145, t6, t168, t59 * mrSges(5,2) + t154 * t104 + t159 * t90 + t85 * t62 + t133, t135, t134;];
f_new  = t3;
