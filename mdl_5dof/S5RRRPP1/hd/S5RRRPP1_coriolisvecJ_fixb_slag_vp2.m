% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:24
% EndTime: 2019-12-31 20:49:30
% DurationCPUTime: 2.52s
% Computational Cost: add. (2403->287), mult. (4385->380), div. (0->0), fcn. (2535->6), ass. (0->156)
t236 = Ifges(5,1) + Ifges(6,1);
t234 = Ifges(5,5) + Ifges(6,4);
t238 = mrSges(6,1) + mrSges(5,1);
t237 = -Ifges(6,5) + Ifges(5,4);
t142 = cos(qJ(3));
t139 = sin(pkin(8));
t140 = sin(qJ(3));
t190 = t139 * t140;
t195 = cos(pkin(8));
t150 = t195 * t142 - t190;
t108 = t150 * qJD(3);
t136 = qJD(1) + qJD(2);
t88 = t136 * t108;
t226 = t88 / 0.2e1;
t230 = -Ifges(5,6) + Ifges(6,6);
t95 = t150 * t136;
t215 = Ifges(6,5) * t95;
t90 = Ifges(5,4) * t95;
t172 = t195 * t140;
t112 = t139 * t142 + t172;
t96 = t112 * t136;
t233 = t234 * qJD(3) + t236 * t96 - t215 + t90;
t141 = sin(qJ(2));
t200 = pkin(1) * qJD(1);
t182 = t141 * t200;
t117 = pkin(7) * t136 + t182;
t143 = cos(qJ(2));
t199 = pkin(1) * qJD(2);
t179 = qJD(1) * t199;
t166 = t143 * t179;
t121 = t142 * t166;
t187 = qJD(3) * t140;
t70 = -t117 * t187 + t121;
t186 = qJD(3) * t142;
t71 = -t117 * t186 - t140 * t166;
t157 = -t140 * t71 + t142 * t70;
t209 = t95 * mrSges(5,3);
t210 = t95 * mrSges(6,2);
t75 = qJD(3) * mrSges(6,3) + t210;
t204 = -qJD(3) * mrSges(5,2) + t209 + t75;
t208 = t96 * mrSges(5,3);
t217 = mrSges(6,2) * t96;
t203 = t238 * qJD(3) - t208 - t217;
t134 = t142 * qJD(4);
t205 = -qJ(4) - pkin(7);
t173 = qJD(3) * t205;
t103 = t140 * t173 + t134;
t149 = -t140 * qJD(4) + t142 * t173;
t181 = t143 * t200;
t229 = -t103 * t139 + t112 * t181 + t149 * t195;
t228 = t103 * t195 + t139 * t149 - t150 * t181;
t192 = t136 * t140;
t115 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t192;
t191 = t136 * t142;
t116 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t191;
t227 = t115 * t142 + t116 * t140;
t225 = t95 / 0.2e1;
t224 = -t95 / 0.2e1;
t222 = t96 / 0.2e1;
t129 = pkin(1) * t141 + pkin(7);
t135 = t142 * qJ(4);
t110 = t129 * t142 + t135;
t189 = -qJ(4) - t129;
t169 = t189 * t140;
t58 = t110 * t139 - t169 * t195;
t168 = qJ(4) * t136 + t117;
t154 = qJD(3) * t168;
t145 = (-qJD(4) * t136 - t166) * t140 - t142 * t154;
t52 = t134 * t136 - t140 * t154 + t121;
t6 = t139 * t52 - t145 * t195;
t221 = t58 * t6;
t122 = pkin(7) * t142 + t135;
t67 = t122 * t139 - t172 * t205;
t220 = t6 * t67;
t216 = Ifges(5,4) * t96;
t214 = pkin(1) * t143;
t213 = pkin(3) * t139;
t212 = t6 * t112;
t211 = t88 * mrSges(6,2);
t207 = -qJD(3) / 0.2e1;
t206 = qJD(3) / 0.2e1;
t7 = t139 * t145 + t195 * t52;
t84 = t168 * t142;
t65 = t195 * t84;
t83 = t168 * t140;
t69 = qJD(3) * pkin(3) - t83;
t30 = t139 * t69 + t65;
t202 = Ifges(4,4) * t140;
t198 = t139 * t84;
t125 = t141 * t179;
t177 = t136 * t187;
t104 = pkin(3) * t177 + t125;
t188 = qJD(3) * t136;
t185 = -qJD(1) - t136;
t184 = -qJD(2) + t136;
t183 = pkin(3) * t192;
t180 = t143 * t199;
t132 = pkin(3) * t187;
t131 = -t142 * pkin(3) - pkin(2);
t178 = t195 * pkin(3);
t107 = t112 * qJD(3);
t87 = t136 * t107;
t41 = t87 * mrSges(5,1) + t88 * mrSges(5,2);
t40 = t87 * mrSges(6,1) - t88 * mrSges(6,3);
t170 = t117 * (t140 ^ 2 + t142 ^ 2);
t167 = qJD(3) * t189;
t59 = t110 * t195 + t139 * t169;
t164 = t58 * t88 - t59 * t87;
t68 = t122 * t195 + t190 * t205;
t163 = t67 * t88 - t68 * t87;
t24 = qJD(3) * qJ(5) + t30;
t5 = qJD(3) * qJD(5) + t7;
t162 = -t24 * t107 + t150 * t5;
t29 = t195 * t69 - t198;
t161 = -t29 * t108 + t212;
t160 = -mrSges(4,1) * t142 + mrSges(4,2) * t140;
t159 = mrSges(4,1) * t140 + mrSges(4,2) * t142;
t158 = Ifges(4,5) * t142 - Ifges(4,6) * t140;
t155 = t140 * t115 - t142 * t116;
t153 = t140 * (Ifges(4,1) * t142 - t202);
t152 = (t142 * Ifges(4,2) + t202) * t136;
t151 = t159 * qJD(3);
t60 = -pkin(4) * t150 - t112 * qJ(5) + t131;
t32 = pkin(4) * t107 - qJ(5) * t108 - qJD(5) * t112 + t132;
t92 = t131 * t136 + qJD(4) - t181;
t146 = (-qJD(4) - t180) * t140 + t142 * t167;
t118 = -t136 * pkin(2) - t181;
t21 = -t95 * pkin(4) - t96 * qJ(5) + t92;
t23 = -qJD(3) * pkin(4) + qJD(5) - t29;
t89 = Ifges(6,5) * t96;
t43 = Ifges(6,6) * qJD(3) - Ifges(6,3) * t95 + t89;
t44 = Ifges(5,2) * t95 + Ifges(5,6) * qJD(3) + t216;
t8 = pkin(4) * t87 - qJ(5) * t88 - qJD(5) * t96 + t104;
t97 = Ifges(4,6) * qJD(3) + t152;
t124 = Ifges(4,4) * t191;
t98 = Ifges(4,1) * t192 + Ifges(4,5) * qJD(3) + t124;
t144 = t153 * t188 + qJD(3) ^ 2 * t158 / 0.2e1 + t118 * t151 + mrSges(6,2) * t212 + t160 * t125 + (t236 * t88 - t237 * t87) * t112 / 0.2e1 - (t152 + t97) * t187 / 0.2e1 + t157 * mrSges(4,3) + (t43 / 0.2e1 - t44 / 0.2e1 + t92 * mrSges(5,1) + t21 * mrSges(6,1) - t30 * mrSges(5,3) + Ifges(6,3) * t224 - Ifges(5,2) * t225 - t237 * t222 + t230 * t206) * t107 + (t104 * mrSges(5,2) - t8 * mrSges(6,3) + (-Ifges(5,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t87 + t236 * t226) * t112 + (-t104 * mrSges(5,1) - t8 * mrSges(6,1) + t7 * mrSges(5,3) + (-Ifges(6,3) - Ifges(5,2)) * t87 + 0.2e1 * t237 * t226) * t150 + (t98 + (0.3e1 * Ifges(4,4) * t142 + (Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t140) * t136) * t186 / 0.2e1 + (t233 / 0.2e1 + t92 * mrSges(5,2) + t23 * mrSges(6,2) - t21 * mrSges(6,3) + Ifges(5,4) * t225 + Ifges(6,5) * t224 + t234 * t206 + t236 * t222) * t108;
t133 = t141 * t199;
t130 = -pkin(2) - t214;
t128 = -t178 - pkin(4);
t126 = qJ(5) + t213;
t119 = t131 - t214;
t114 = t133 + t132;
t105 = t160 * t136;
t99 = t136 * t151;
t63 = t140 * t167 + t142 * t180 + t134;
t55 = t60 - t214;
t51 = -mrSges(5,1) * t95 + mrSges(5,2) * t96;
t50 = -mrSges(6,1) * t95 - mrSges(6,3) * t96;
t39 = pkin(4) * t96 - qJ(5) * t95 + t183;
t38 = -t195 * t83 - t198;
t37 = -t139 * t83 + t65;
t22 = t133 + t32;
t12 = t139 * t146 + t195 * t63;
t11 = t139 * t63 - t146 * t195;
t1 = [m(5) * (t104 * t119 - t11 * t29 + t114 * t92 + t12 * t30 + t59 * t7 + t221) + m(6) * (t11 * t23 + t12 * t24 + t21 * t22 + t5 * t59 + t55 * t8 + t221) + (t161 + t164) * mrSges(5,3) + (t162 + t164) * mrSges(6,2) + ((m(4) * (qJD(1) * t130 + t118) + t105 + t185 * mrSges(3,1)) * t141 + (m(4) * t170 + mrSges(3,2) * t185 - t155) * t143) * t199 + t130 * t99 + t114 * t51 + t119 * t41 + t22 * t50 + t55 * t40 - t203 * t11 + t204 * t12 + (m(4) * t157 - qJD(3) * t227) * t129 + t144; (t140 * pkin(3) * t51 - pkin(7) * t227) * qJD(3) + (t161 + t163) * mrSges(5,3) + t131 * t41 - pkin(2) * t99 + t32 * t50 + t60 * t40 + (t162 + t163) * mrSges(6,2) + t144 + ((mrSges(3,2) * t184 + t155) * t143 + (mrSges(3,1) * t184 - t105 - t50 - t51) * t141) * t200 + t228 * t204 + t229 * t203 + (t5 * t68 + t60 * t8 + t220 + t228 * t24 - t229 * t23 + (-t182 + t32) * t21) * m(6) + (t104 * t131 + t68 * t7 + t220 + (t132 - t182) * t92 + t228 * t30 + t229 * t29) * m(5) + (-pkin(2) * t125 + t157 * pkin(7) - (t118 * t141 + t143 * t170) * t200) * m(4); (-mrSges(5,3) * t178 + t234) * t88 + (t230 * t96 + t234 * t95) * t207 - t238 * t6 - t21 * (mrSges(6,1) * t96 - mrSges(6,3) * t95) - t92 * (t96 * mrSges(5,1) + t95 * mrSges(5,2)) + (Ifges(6,3) * t96 + t215) * t225 + (-Ifges(5,2) * t96 + t233 + t90) * t224 + (-t92 * t183 + t29 * t37 - t30 * t38 + (t139 * t7 - t195 * t6) * pkin(3)) * m(5) + t29 * t209 + t97 * t192 / 0.2e1 - t158 * t188 / 0.2e1 - t51 * t183 - Ifges(4,6) * t177 - t70 * mrSges(4,2) + t71 * mrSges(4,1) + qJD(5) * t75 - t39 * t50 + t203 * t37 - t204 * t38 + (-mrSges(6,2) * t126 - mrSges(5,3) * t213 + t230) * t87 + t227 * t117 - (-Ifges(4,2) * t192 + t124 + t98) * t191 / 0.2e1 + t5 * mrSges(6,3) - t7 * mrSges(5,2) - (t236 * t95 - t216 + t43 + t89) * t96 / 0.2e1 + t30 * t208 - t23 * t210 + t128 * t211 + t24 * t217 + t44 * t222 + (Ifges(4,5) * t186 - t118 * t159 - t153 * t136 / 0.2e1) * t136 + (t126 * t5 + t128 * t6 - t21 * t39 - t23 * t37 + (-t38 + qJD(5)) * t24) * m(6); t203 * t96 - t204 * t95 + t40 + t41 + (-t23 * t96 - t24 * t95 + t8) * m(6) + (t29 * t96 - t30 * t95 + t104) * m(5); t211 - qJD(3) * t75 + t96 * t50 + 0.2e1 * (t6 / 0.2e1 + t24 * t207 + t21 * t222) * m(6);];
tauc = t1(:);
