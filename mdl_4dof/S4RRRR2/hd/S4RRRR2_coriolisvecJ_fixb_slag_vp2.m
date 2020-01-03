% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:08
% EndTime: 2019-12-31 17:23:10
% DurationCPUTime: 1.02s
% Computational Cost: add. (1650->193), mult. (3022->296), div. (0->0), fcn. (1711->6), ass. (0->122)
t100 = qJD(3) + qJD(4);
t104 = sin(qJ(4));
t105 = sin(qJ(3));
t107 = cos(qJ(4));
t108 = cos(qJ(3));
t82 = -t104 * t105 + t107 * t108;
t45 = t100 * t82;
t172 = t45 / 0.2e1;
t101 = qJD(1) + qJD(2);
t70 = t82 * t101;
t171 = -t70 / 0.2e1;
t109 = cos(qJ(2));
t154 = pkin(1) * qJD(1);
t135 = t109 * t154;
t163 = -pkin(7) - pkin(6);
t93 = t163 * t105;
t99 = t108 * pkin(7);
t94 = pkin(6) * t108 + t99;
t55 = t104 * t93 + t107 * t94;
t83 = t104 * t108 + t105 * t107;
t131 = qJD(3) * t163;
t84 = t105 * t131;
t85 = t108 * t131;
t170 = -qJD(4) * t55 - t104 * t84 + t107 * t85 + t83 * t135;
t54 = -t104 * t94 + t107 * t93;
t169 = qJD(4) * t54 + t104 * t85 + t107 * t84 - t82 * t135;
t106 = sin(qJ(2));
t89 = pkin(6) * t101 + t106 * t154;
t168 = (t105 ^ 2 + t108 ^ 2) * t89;
t71 = t83 * t101;
t39 = -mrSges(5,1) * t70 + mrSges(5,2) * t71;
t98 = -t108 * pkin(3) - pkin(2);
t73 = t101 * t98 - t135;
t167 = -m(5) * t73 - t39;
t164 = t71 / 0.2e1;
t96 = pkin(1) * t106 + pkin(6);
t161 = -pkin(7) - t96;
t160 = mrSges(5,3) * t70;
t159 = Ifges(5,4) * t71;
t158 = pkin(1) * t109;
t157 = t71 * mrSges(5,3);
t156 = Ifges(4,4) * t105;
t153 = pkin(1) * qJD(2);
t152 = qJD(3) * pkin(3);
t129 = pkin(7) * t101 + t89;
t62 = t129 * t108;
t151 = t104 * t62;
t132 = qJD(1) * t153;
t124 = t109 * t132;
t117 = t105 * t124;
t140 = qJD(3) * t108;
t50 = -t140 * t89 - t117;
t150 = t105 * t50;
t149 = t107 * t62;
t141 = qJD(3) * t105;
t92 = t108 * t124;
t49 = -t141 * t89 + t92;
t148 = t108 * t49;
t146 = Ifges(4,5) * qJD(3);
t145 = Ifges(4,6) * qJD(3);
t144 = t101 * t105;
t143 = t101 * t108;
t142 = qJD(2) * t106;
t139 = qJD(4) * t104;
t138 = qJD(4) * t107;
t137 = -qJD(1) - t101;
t136 = -qJD(2) + t101;
t134 = t109 * t153;
t133 = pkin(3) * t141;
t130 = t101 * t141;
t128 = qJD(3) * t161;
t125 = t106 * t132;
t61 = t129 * t105;
t51 = -t61 + t152;
t26 = t107 * t51 - t151;
t27 = t104 * t51 + t149;
t121 = qJD(3) * t129;
t40 = -t105 * t121 + t92;
t41 = -t108 * t121 - t117;
t4 = -qJD(4) * t27 - t104 * t40 + t107 * t41;
t122 = -t26 * t45 - t4 * t83;
t120 = -mrSges(4,1) * t108 + mrSges(4,2) * t105;
t80 = t161 * t105;
t81 = t108 * t96 + t99;
t42 = -t104 * t81 + t107 * t80;
t43 = t104 * t80 + t107 * t81;
t87 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t144;
t88 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t143;
t119 = t105 * t88 + t108 * t87;
t118 = t105 * t87 - t108 * t88;
t116 = (Ifges(4,2) * t108 + t156) * t101;
t115 = (mrSges(4,1) * t105 + mrSges(4,2) * t108) * qJD(3);
t3 = qJD(4) * t26 + t104 * t41 + t107 * t40;
t31 = Ifges(5,2) * t70 + Ifges(5,6) * t100 + t159;
t65 = Ifges(5,4) * t70;
t32 = Ifges(5,1) * t71 + Ifges(5,5) * t100 + t65;
t36 = t45 * t101;
t46 = t100 * t83;
t37 = t46 * t101;
t112 = t4 * mrSges(5,1) - t3 * mrSges(5,2) + t26 * t160 + t31 * t164 - t73 * (mrSges(5,1) * t71 + mrSges(5,2) * t70) - t71 * (Ifges(5,1) * t70 - t159) / 0.2e1 - t100 * (Ifges(5,5) * t70 - Ifges(5,6) * t71) / 0.2e1 - Ifges(5,6) * t37 + Ifges(5,5) * t36 + (-Ifges(5,2) * t71 + t32 + t65) * t171;
t111 = m(4) * (t148 - t150) - t119 * qJD(3);
t68 = t116 + t145;
t95 = Ifges(4,4) * t143;
t69 = Ifges(4,1) * t144 + t146 + t95;
t76 = pkin(3) * t130 + t125;
t90 = -t101 * pkin(2) - t135;
t110 = (Ifges(4,1) * t108 - t156) * t130 + t120 * t125 + qJD(3) ^ 2 * (Ifges(4,5) * t108 - Ifges(4,6) * t105) / 0.2e1 + t100 * (Ifges(5,5) * t45 - Ifges(5,6) * t46) / 0.2e1 + t76 * (-mrSges(5,1) * t82 + mrSges(5,2) * t83) + t73 * (mrSges(5,1) * t46 + mrSges(5,2) * t45) + t32 * t172 - t46 * t31 / 0.2e1 + mrSges(4,3) * t148 + t90 * t115 - (t116 + t68) * t141 / 0.2e1 + (t46 * t171 - t82 * t37) * Ifges(5,2) + (t45 * t164 + t36 * t83) * Ifges(5,1) + (-t27 * t46 + t3 * t82) * mrSges(5,3) + (-t46 * t164 + t70 * t172 + t82 * t36 - t37 * t83) * Ifges(5,4) + (t69 + (0.3e1 * Ifges(4,4) * t108 + (Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t105) * t101) * t140 / 0.2e1;
t97 = -pkin(2) - t158;
t91 = t98 - t158;
t86 = pkin(1) * t142 + t133;
t77 = t120 * t101;
t72 = t101 * t115;
t57 = -t105 * t134 + t108 * t128;
t56 = t105 * t128 + t108 * t134;
t48 = mrSges(5,1) * t100 - t157;
t47 = -mrSges(5,2) * t100 + t160;
t29 = -t107 * t61 - t151;
t28 = t104 * t61 - t149;
t9 = mrSges(5,1) * t37 + mrSges(5,2) * t36;
t8 = -qJD(4) * t43 - t104 * t56 + t107 * t57;
t7 = qJD(4) * t42 + t104 * t57 + t107 * t56;
t1 = [m(5) * (t26 * t8 + t27 * t7 + t3 * t43 + t4 * t42 + t73 * t86 + t76 * t91) + t110 + t111 * t96 + (-t36 * t42 - t37 * t43 + t122) * mrSges(5,3) + t91 * t9 + t97 * t72 + t86 * t39 + t7 * t47 + t8 * t48 - mrSges(4,3) * t150 + ((m(4) * (qJD(1) * t97 + t90) + t77 + t137 * mrSges(3,1)) * t106 + (m(4) * t168 + t137 * mrSges(3,2) - t118) * t109) * t153; t170 * t48 + t169 * t47 + (-t36 * t54 - t37 * t55 + t122) * mrSges(5,3) + t110 + t111 * pkin(6) + t98 * t9 - pkin(2) * t72 + (-t50 * mrSges(4,3) + t152 * t39) * t105 + ((mrSges(3,2) * t136 + t118) * t109 + (mrSges(3,1) * t136 + t167 - t77) * t106 + (-pkin(2) * t142 - t106 * t90 - t109 * t168) * m(4)) * t154 + (t73 * t133 + t169 * t27 + t170 * t26 + t3 * t55 + t4 * t54 + t76 * t98) * m(5); t112 + ((t146 / 0.2e1 - t90 * mrSges(4,2) - t95 / 0.2e1 - t69 / 0.2e1) * t108 + (-t145 / 0.2e1 - t90 * mrSges(4,1) + t68 / 0.2e1 + (t156 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t108) * t101 + t167 * pkin(3)) * t105) * t101 + (m(5) * (t104 * t3 + t107 * t4 + t138 * t27 - t139 * t26) + t47 * t138 - t48 * t139 + (-t104 * t37 - t107 * t36) * mrSges(5,3)) * pkin(3) - m(5) * (t26 * t28 + t27 * t29) + t119 * t89 - t29 * t47 - t28 * t48 - t49 * mrSges(4,2) + t50 * mrSges(4,1) + t27 * t157; t112 - t26 * t47 + (t48 + t157) * t27;];
tauc = t1(:);
