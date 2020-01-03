% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:31
% EndTime: 2020-01-03 12:03:38
% DurationCPUTime: 2.63s
% Computational Cost: add. (4594->284), mult. (7935->397), div. (0->0), fcn. (5542->8), ass. (0->140)
t159 = cos(qJ(2));
t186 = qJD(1) * pkin(1);
t177 = t159 * t186;
t163 = qJD(3) - t177;
t152 = sin(pkin(9));
t153 = cos(pkin(9));
t155 = sin(qJ(4));
t158 = cos(qJ(4));
t136 = t152 * t158 + t153 * t155;
t140 = (-pkin(7) - qJ(3)) * t152;
t146 = t153 * pkin(7);
t141 = qJ(3) * t153 + t146;
t95 = t155 * t140 + t158 * t141;
t208 = -qJD(4) * t95 - t163 * t136;
t180 = t153 * t158;
t135 = -t152 * t155 + t180;
t178 = qJD(4) * t158;
t207 = -t135 * t177 + t140 * t178 + qJD(3) * t180 + (-qJD(3) * t152 - qJD(4) * t141) * t155;
t124 = t135 * qJD(4);
t194 = pkin(8) * t124;
t222 = -t194 + t208;
t125 = t136 * qJD(4);
t121 = t125 * pkin(8);
t221 = t121 - t207;
t150 = qJD(4) + qJD(5);
t151 = qJD(1) + qJD(2);
t112 = t135 * t151;
t113 = t136 * t151;
t154 = sin(qJ(5));
t157 = cos(qJ(5));
t167 = t157 * t112 - t113 * t154;
t79 = t112 * t154 + t113 * t157;
t196 = Ifges(6,4) * t79;
t156 = sin(qJ(2));
t173 = t156 * t186;
t138 = qJ(3) * t151 + t173;
t171 = pkin(7) * t151 + t138;
t101 = t171 * t152;
t102 = t171 * t153;
t65 = -t101 * t155 + t102 * t158;
t47 = pkin(8) * t112 + t65;
t185 = t154 * t47;
t64 = -t158 * t101 - t102 * t155;
t46 = -pkin(8) * t113 + t64;
t45 = qJD(4) * pkin(4) + t46;
t20 = t157 * t45 - t185;
t184 = t157 * t47;
t21 = t154 * t45 + t184;
t108 = t151 * t125;
t172 = qJD(2) * t186;
t164 = t159 * t172;
t134 = t151 * qJD(3) + t164;
t39 = t134 * t180 - t101 * t178 + (-qJD(4) * t102 - t134 * t152) * t155;
t25 = -pkin(8) * t108 + t39;
t107 = t151 * t124;
t40 = -qJD(4) * t65 - t134 * t136;
t26 = -pkin(8) * t107 + t40;
t3 = qJD(5) * t20 + t154 * t26 + t157 * t25;
t73 = Ifges(6,4) * t167;
t31 = t79 * Ifges(6,1) + t150 * Ifges(6,5) + t73;
t36 = qJD(5) * t167 + t107 * t157 - t108 * t154;
t37 = -qJD(5) * t79 - t107 * t154 - t108 * t157;
t4 = -qJD(5) * t21 - t154 * t25 + t157 * t26;
t145 = -t153 * pkin(3) - pkin(2);
t111 = t145 * t151 + t163;
t83 = -t112 * pkin(4) + t111;
t220 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + Ifges(6,5) * t36 + Ifges(6,6) * t37 - (Ifges(6,5) * t167 - Ifges(6,6) * t79) * t150 / 0.2e1 + (t167 * t20 + t21 * t79) * mrSges(6,3) - (-Ifges(6,2) * t79 + t31 + t73) * t167 / 0.2e1 - t83 * (mrSges(6,1) * t79 + mrSges(6,2) * t167) - (Ifges(6,1) * t167 - t196) * t79 / 0.2e1;
t219 = t124 / 0.2e1;
t218 = -t125 / 0.2e1;
t217 = t167 / 0.2e1;
t179 = t152 ^ 2 + t153 ^ 2;
t169 = t179 * t134;
t30 = Ifges(6,2) * t167 + t150 * Ifges(6,6) + t196;
t215 = t30 / 0.2e1;
t193 = pkin(8) * t136;
t94 = t158 * t140 - t141 * t155;
t84 = t94 - t193;
t129 = t135 * pkin(8);
t85 = t129 + t95;
t42 = t154 * t84 + t157 * t85;
t214 = -qJD(5) * t42 + t221 * t154 + t222 * t157;
t41 = -t154 * t85 + t157 * t84;
t213 = qJD(5) * t41 + t222 * t154 - t221 * t157;
t161 = -mrSges(4,1) * t153 + mrSges(4,2) * t152;
t206 = mrSges(5,1) * t112 - mrSges(5,2) * t113 - t161 * t151;
t144 = pkin(1) * t156 + qJ(3);
t126 = (-pkin(7) - t144) * t152;
t127 = t144 * t153 + t146;
t89 = t155 * t126 + t158 * t127;
t205 = -t64 * t124 - t40 * t136;
t90 = t135 * t157 - t136 * t154;
t48 = qJD(5) * t90 + t124 * t157 - t125 * t154;
t91 = t135 * t154 + t136 * t157;
t204 = -t20 * t48 - t4 * t91;
t201 = t79 / 0.2e1;
t198 = t113 / 0.2e1;
t195 = pkin(1) * t159;
t192 = t125 * pkin(4);
t189 = mrSges(4,3) * t151;
t188 = Ifges(5,4) * t113;
t187 = pkin(1) * qJD(2);
t176 = t156 * t187;
t175 = t159 * t187;
t12 = -t37 * mrSges(6,1) + t36 * mrSges(6,2);
t170 = t179 * mrSges(4,3);
t72 = t108 * mrSges(5,1) + t107 * mrSges(5,2);
t168 = t179 * t138;
t88 = t158 * t126 - t127 * t155;
t166 = qJD(3) * t179;
t165 = t156 * t172;
t70 = t88 - t193;
t71 = t129 + t89;
t28 = -t154 * t71 + t157 * t70;
t29 = t154 * t70 + t157 * t71;
t110 = -t135 * pkin(4) + t145;
t142 = qJD(3) + t175;
t51 = t126 * t178 + t142 * t180 + (-qJD(4) * t127 - t142 * t152) * t155;
t52 = -qJD(4) * t89 - t136 * t142;
t49 = -qJD(5) * t91 - t124 * t154 - t125 * t157;
t74 = t112 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t188;
t109 = Ifges(5,4) * t112;
t75 = t113 * Ifges(5,1) + Ifges(5,5) * qJD(4) + t109;
t92 = pkin(4) * t108 + t165;
t160 = t150 * (Ifges(6,5) * t48 + Ifges(6,6) * t49) / 0.2e1 + t75 * t219 + t74 * t218 + t92 * (-mrSges(6,1) * t90 + mrSges(6,2) * t91) + t83 * (-mrSges(6,1) * t49 + mrSges(6,2) * t48) + t48 * t31 / 0.2e1 + t49 * t215 + qJD(4) * (Ifges(5,5) * t124 - Ifges(5,6) * t125) / 0.2e1 + t111 * (mrSges(5,1) * t125 + mrSges(5,2) * t124) + mrSges(4,3) * t169 + (-mrSges(5,1) * t135 + mrSges(5,2) * t136 + t161) * t165 + (t49 * t217 + t90 * t37) * Ifges(6,2) + (-t135 * t108 + t112 * t218) * Ifges(5,2) + (t48 * t201 + t36 * t91) * Ifges(6,1) + (t107 * t136 + t124 * t198) * Ifges(5,1) + (t21 * t49 + t3 * t90) * mrSges(6,3) + (-t65 * t125 + t39 * t135) * mrSges(5,3) + (t49 * t201 + t48 * t217 + t90 * t36 + t37 * t91) * Ifges(6,4) + (t135 * t107 - t108 * t136 + t112 * t219 - t125 * t198) * Ifges(5,4);
t139 = t145 - t195;
t137 = -t151 * pkin(2) + t163;
t103 = t176 + t192;
t99 = t110 - t195;
t98 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t113;
t97 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t112;
t67 = mrSges(6,1) * t150 - mrSges(6,3) * t79;
t66 = -mrSges(6,2) * t150 + mrSges(6,3) * t167;
t44 = t52 - t194;
t43 = -t121 + t51;
t38 = -mrSges(6,1) * t167 + mrSges(6,2) * t79;
t23 = t157 * t46 - t185;
t22 = -t154 * t46 - t184;
t6 = -qJD(5) * t29 - t154 * t43 + t157 * t44;
t5 = qJD(5) * t28 + t154 * t44 + t157 * t43;
t1 = [t139 * t72 + t99 * t12 + t103 * t38 + t51 * t97 + t52 * t98 + t6 * t67 + t5 * t66 + t160 + m(6) * (t103 * t83 + t20 * t6 + t21 * t5 + t28 * t4 + t29 * t3 + t92 * t99) - mrSges(3,1) * t165 + m(5) * (t39 * t89 + t40 * t88 + t51 * t65 + t52 * t64) + m(4) * (t142 * t168 + t144 * t169) + t179 * t142 * t189 + (-t151 * t175 - t164) * mrSges(3,2) + (-t151 * mrSges(3,1) + m(5) * (qJD(1) * t139 + t111) + m(4) * (t137 + (-pkin(2) - t195) * qJD(1)) - t206) * t176 + (-t28 * t36 + t29 * t37 + t204) * mrSges(6,3) + (-t88 * t107 - t89 * t108 + t205) * mrSges(5,3); t145 * t72 + t110 * t12 + t38 * t192 + t160 + t208 * t98 + t207 * t97 + (-t41 * t36 + t42 * t37 + t204) * mrSges(6,3) + (-t94 * t107 - t95 * t108 + t205) * mrSges(5,3) + t214 * t67 + t213 * t66 + ((-mrSges(3,2) * qJD(2) + (mrSges(3,2) - t170) * t151) * t159 + (-t38 + (-qJD(2) + t151) * mrSges(3,1) + t206) * t156) * t186 + t166 * t189 + (t110 * t92 + t3 * t42 + t4 * t41 + (-t173 + t192) * t83 + t213 * t21 + t214 * t20) * m(6) + (-t111 * t173 + t145 * t165 + t207 * t65 + t208 * t64 + t39 * t95 + t40 * t94) * m(5) + (-pkin(2) * t165 + qJ(3) * t169 + t138 * t166 - (t137 * t156 + t159 * t168) * t186) * m(4); -t112 * t97 + t113 * t98 - t167 * t66 + t79 * t67 + (m(4) + m(5)) * t165 - m(5) * (t112 * t65 - t113 * t64) + t72 + t12 + (-m(4) * t168 - t170 * t151) * t151 + (-t167 * t21 + t20 * t79 + t92) * m(6); (t112 * t64 + t113 * t65) * mrSges(5,3) - t111 * (mrSges(5,1) * t113 + mrSges(5,2) * t112) - qJD(4) * (Ifges(5,5) * t112 - Ifges(5,6) * t113) / 0.2e1 + Ifges(5,5) * t107 - Ifges(5,6) * t108 - t64 * t97 + t65 * t98 - t23 * t66 - t22 * t67 - t39 * mrSges(5,2) + t40 * mrSges(5,1) + t74 * t198 - (-Ifges(5,2) * t113 + t109 + t75) * t112 / 0.2e1 + t79 * t215 - m(6) * (t20 * t22 + t21 * t23) - t113 * (Ifges(5,1) * t112 - t188) / 0.2e1 + (-t113 * t38 + (-t154 * t67 + t157 * t66) * qJD(5) + (t154 * t37 - t157 * t36) * mrSges(6,3) + (-t113 * t83 + t154 * t3 + t157 * t4 + (-t154 * t20 + t157 * t21) * qJD(5)) * m(6)) * pkin(4) + t220; -t20 * t66 + t30 * t201 + t21 * t67 + t220;];
tauc = t1(:);
