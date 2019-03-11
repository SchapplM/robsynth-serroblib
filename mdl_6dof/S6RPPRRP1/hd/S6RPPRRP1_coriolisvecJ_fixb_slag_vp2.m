% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:57:37
% EndTime: 2019-03-09 01:57:46
% DurationCPUTime: 4.30s
% Computational Cost: add. (4701->376), mult. (11648->486), div. (0->0), fcn. (8136->8), ass. (0->173)
t256 = Ifges(6,4) + Ifges(7,4);
t257 = Ifges(6,1) + Ifges(7,1);
t247 = Ifges(6,5) + Ifges(7,5);
t255 = Ifges(6,2) + Ifges(7,2);
t246 = Ifges(7,6) + Ifges(6,6);
t258 = Ifges(5,2) / 0.2e1;
t129 = sin(pkin(10));
t131 = cos(pkin(10));
t134 = sin(qJ(4));
t136 = cos(qJ(4));
t118 = t129 * t136 + t131 * t134;
t111 = t118 * qJD(1);
t133 = sin(qJ(5));
t135 = cos(qJ(5));
t93 = qJD(4) * t135 - t111 * t133;
t254 = t256 * t93;
t253 = t256 * t135;
t252 = t256 * t133;
t94 = qJD(4) * t133 + t111 * t135;
t251 = t256 * t94;
t144 = -cos(pkin(9)) * pkin(1) - pkin(3) * t131 - pkin(2);
t107 = qJD(1) * t144 + qJD(3);
t117 = t129 * t134 - t136 * t131;
t110 = t117 * qJD(1);
t122 = sin(pkin(9)) * pkin(1) + qJ(3);
t119 = t122 * qJD(1);
t126 = t131 * qJD(2);
t198 = pkin(7) * qJD(1);
t91 = t126 + (-t119 - t198) * t129;
t101 = t129 * qJD(2) + t131 * t119;
t92 = t131 * t198 + t101;
t50 = t134 * t91 + t136 * t92;
t47 = qJD(4) * pkin(8) + t50;
t58 = pkin(4) * t110 - pkin(8) * t111 + t107;
t20 = t133 * t58 + t135 * t47;
t11 = qJ(6) * t93 + t20;
t19 = -t133 * t47 + t135 * t58;
t148 = t133 * t20 + t135 * t19;
t161 = mrSges(7,1) * t133 + mrSges(7,2) * t135;
t163 = mrSges(6,1) * t133 + mrSges(6,2) * t135;
t216 = t135 / 0.2e1;
t219 = -t133 / 0.2e1;
t240 = qJD(5) + t110;
t221 = -t240 / 0.2e1;
t223 = t94 / 0.2e1;
t226 = -t93 / 0.2e1;
t233 = t135 * t257 - t252;
t234 = -t133 * t255 + t253;
t235 = -t133 * t246 + t135 * t247;
t237 = t247 * t240 + t257 * t94 + t254;
t245 = t240 * t246 + t255 * t93 + t251;
t49 = -t134 * t92 + t136 * t91;
t46 = -qJD(4) * pkin(4) - t49;
t28 = -pkin(5) * t93 + qJD(6) + t46;
t10 = -qJ(6) * t94 + t19;
t9 = pkin(5) * t240 + t10;
t231 = t148 * mrSges(6,3) + (t11 * t133 + t135 * t9) * mrSges(7,3) - t161 * t28 - t163 * t46 + t234 * t226 - t233 * t223 + t235 * t221 - t245 * t219 - t237 * t216;
t249 = t111 * Ifges(5,1) / 0.2e1;
t250 = t107 * mrSges(5,2) - t49 * mrSges(5,3) - Ifges(5,4) * t110 + Ifges(5,5) * qJD(4) - t231 + t249;
t248 = t110 * t258;
t171 = qJD(1) * (t129 ^ 2 + t131 ^ 2);
t244 = mrSges(4,3) * t171;
t243 = t133 * t247 + t135 * t246;
t242 = t135 * t255 + t252;
t241 = t133 * t257 + t253;
t224 = -t94 / 0.2e1;
t113 = t118 * qJD(4);
t103 = qJD(1) * t113;
t112 = t117 * qJD(4);
t102 = qJD(1) * t112;
t64 = qJD(5) * t93 - t102 * t135;
t65 = -qJD(5) * t94 + t102 * t133;
t239 = t103 * t246 + t255 * t65 + t256 * t64;
t238 = t103 * t247 + t256 * t65 + t257 * t64;
t209 = pkin(7) + t122;
t114 = t209 * t129;
t115 = t209 * t131;
t236 = -t136 * t114 - t115 * t134;
t181 = qJD(5) * t135;
t182 = qJD(5) * t133;
t142 = t117 * qJD(3);
t32 = -qJD(1) * t142 + qJD(4) * t49;
t75 = pkin(4) * t103 + pkin(8) * t102;
t4 = t133 * t75 + t135 * t32 + t58 * t181 - t182 * t47;
t5 = -qJD(5) * t20 - t133 * t32 + t135 * t75;
t168 = -t133 * t5 + t135 * t4;
t176 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t177 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t178 = -Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1;
t229 = t176 * t240 + t177 * t93 - t178 * t94 - t11 * mrSges(7,2) - t20 * mrSges(6,2) - t50 * mrSges(5,3) - Ifges(5,6) * qJD(4) - t111 * Ifges(5,4) + t248 + t107 * mrSges(5,1) + t19 * mrSges(6,1) + t9 * mrSges(7,1) - t246 * t226 - t247 * t224 - (Ifges(7,3) + Ifges(6,3)) * t221;
t228 = t64 / 0.2e1;
t227 = t65 / 0.2e1;
t222 = t103 / 0.2e1;
t143 = t118 * qJD(3);
t33 = qJD(1) * t143 + qJD(4) * t50;
t213 = t33 * t236;
t208 = -qJ(6) - pkin(8);
t41 = mrSges(7,1) * t103 - mrSges(7,3) * t64;
t42 = mrSges(6,1) * t103 - mrSges(6,3) * t64;
t207 = -t41 - t42;
t43 = -mrSges(7,2) * t103 + mrSges(7,3) * t65;
t44 = -mrSges(6,2) * t103 + mrSges(6,3) * t65;
t206 = t43 + t44;
t86 = pkin(4) * t111 + pkin(8) * t110;
t26 = t133 * t86 + t135 * t49;
t205 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t93 + mrSges(6,2) * t94 + t111 * mrSges(5,3);
t66 = -mrSges(7,2) * t240 + mrSges(7,3) * t93;
t67 = -mrSges(6,2) * t240 + mrSges(6,3) * t93;
t204 = t66 + t67;
t68 = mrSges(7,1) * t240 - mrSges(7,3) * t94;
t69 = mrSges(6,1) * t240 - mrSges(6,3) * t94;
t203 = t68 + t69;
t85 = -t114 * t134 + t115 * t136;
t76 = t135 * t85;
t77 = pkin(4) * t117 - pkin(8) * t118 + t144;
t30 = t133 * t77 + t76;
t193 = t117 * t33;
t189 = qJ(6) * t135;
t188 = t112 * t133;
t187 = t112 * t135;
t185 = t118 * t133;
t184 = t133 * t110;
t52 = qJD(4) * t236 - t142;
t87 = pkin(4) * t113 + pkin(8) * t112;
t180 = t133 * t87 + t135 * t52 + t77 * t181;
t56 = -mrSges(7,1) * t93 + mrSges(7,2) * t94;
t179 = t56 + t205;
t175 = t118 * t181;
t22 = -t65 * mrSges(7,1) + t64 * mrSges(7,2);
t174 = t103 * mrSges(5,1) - t102 * mrSges(5,2);
t25 = -t133 * t49 + t135 * t86;
t173 = -t133 * t52 + t135 * t87;
t29 = -t133 * t85 + t135 * t77;
t172 = qJD(5) * t208;
t1 = pkin(5) * t103 - qJ(6) * t64 - qJD(6) * t94 + t5;
t2 = qJ(6) * t65 + qJD(6) * t93 + t4;
t170 = -t1 * t135 - t133 * t2;
t169 = -t1 * t133 + t135 * t2;
t167 = -t4 * t133 - t5 * t135;
t166 = t11 * t135 - t133 * t9;
t164 = mrSges(6,1) * t135 - mrSges(6,2) * t133;
t162 = mrSges(7,1) * t135 - mrSges(7,2) * t133;
t147 = t133 * t19 - t135 * t20;
t146 = -(-t119 * t129 + t126) * t129 + t101 * t131;
t145 = qJ(6) * t112 - qJD(6) * t118;
t141 = t5 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2);
t53 = qJD(4) * t85 + t143;
t124 = -pkin(5) * t135 - pkin(4);
t121 = t208 * t135;
t120 = t208 * t133;
t109 = -qJD(6) * t133 + t135 * t172;
t108 = qJD(6) * t135 + t133 * t172;
t99 = Ifges(6,3) * t103;
t98 = Ifges(7,3) * t103;
t95 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t110;
t63 = Ifges(6,5) * t64;
t62 = Ifges(7,5) * t64;
t61 = Ifges(6,6) * t65;
t60 = Ifges(7,6) * t65;
t54 = pkin(5) * t185 - t236;
t34 = -pkin(5) * t184 + t50;
t27 = (t175 - t188) * pkin(5) + t53;
t24 = -qJ(6) * t185 + t30;
t23 = -mrSges(6,1) * t65 + mrSges(6,2) * t64;
t21 = pkin(5) * t117 - t118 * t189 + t29;
t18 = qJ(6) * t184 + t26;
t13 = -pkin(5) * t65 + t33;
t12 = pkin(5) * t111 + t110 * t189 + t25;
t8 = -qJD(5) * t30 + t173;
t7 = -t182 * t85 + t180;
t6 = -qJ(6) * t175 + (-qJD(5) * t85 + t145) * t133 + t180;
t3 = pkin(5) * t113 + t145 * t135 + (-t76 + (qJ(6) * t118 - t77) * t133) * qJD(5) + t173;
t14 = [t21 * t41 + t29 * t42 + t24 * t43 + t30 * t44 + t54 * t22 + t27 * t56 + t6 * t66 + t7 * t67 + t3 * t68 + t8 * t69 - t236 * t23 + t52 * t95 + t144 * t174 + t205 * t53 + (t102 * t236 - t103 * t85) * mrSges(5,3) + m(6) * (t19 * t8 + t20 * t7 + t29 * t5 + t30 * t4 + t46 * t53 - t213) + m(5) * (t32 * t85 - t49 * t53 + t50 * t52 - t213) + m(7) * (t1 * t21 + t11 * t6 + t13 * t54 + t2 * t24 + t27 * t28 + t3 * t9) + (m(4) * (t122 * t171 + t146) + 0.2e1 * t244) * qJD(3) + (t248 + t229) * t113 - (t249 + t250) * t112 + (-t32 * mrSges(5,3) + t62 / 0.2e1 + t60 / 0.2e1 + t98 / 0.2e1 + t63 / 0.2e1 + t61 / 0.2e1 + t99 / 0.2e1 + Ifges(5,4) * t102 + t177 * t65 - t178 * t64 + (Ifges(5,2) + t176) * t103 + t141) * t117 + (t13 * t161 - Ifges(5,1) * t102 - Ifges(5,4) * t103 + (mrSges(5,3) + t163) * t33 + t170 * mrSges(7,3) + t167 * mrSges(6,3) + (mrSges(6,3) * t147 - mrSges(7,3) * t166 + t162 * t28 + t164 * t46 + t242 * t226 + t241 * t224 + t243 * t221 - t245 * t135 / 0.2e1) * qJD(5) + t233 * t228 + t234 * t227 + t235 * t222 + t238 * t216 + (qJD(5) * t237 + t239) * t219) * t118; (-t102 * mrSges(5,3) + t22 + t23) * t117 + t179 * t113 - (-t133 * t203 + t135 * t204 + t95) * t112 + m(5) * (-t112 * t50 - t113 * t49 + t193) + m(6) * (t113 * t46 - t187 * t20 + t188 * t19 + t193) + m(7) * (-t11 * t187 + t113 * t28 + t117 * t13 + t188 * t9) + (-t103 * mrSges(5,3) + t206 * t135 + t207 * t133 + (-t133 * t204 - t135 * t203) * qJD(5) + m(5) * t32 + m(6) * (-t181 * t19 - t182 * t20 + t168) + m(7) * (-t11 * t182 - t181 * t9 + t169)) * t118; t110 * t95 - t179 * t111 + (t204 * t240 - t207) * t135 + (-t203 * t240 + t206) * t133 - m(5) * (-t110 * t50 - t111 * t49) + t174 + (-m(4) * t146 - t244) * qJD(1) + (-t111 * t28 + t166 * t240 - t170) * m(7) + (-t46 * t111 - t147 * t240 - t167) * m(6); t238 * t133 / 0.2e1 + t239 * t216 + ((-m(6) * t148 - t133 * t67 - t135 * t69) * qJD(5) - t133 * t42 + t135 * t44 + m(6) * t168) * pkin(8) + t168 * mrSges(6,3) + t169 * mrSges(7,3) - t13 * t162 + (-mrSges(5,1) - t164) * t33 + (-t18 + t108) * t66 - t229 * t111 + (-t231 + (m(7) * t28 + t56) * t133 * pkin(5)) * qJD(5) + t241 * t228 + t242 * t227 + t243 * t222 - m(7) * (t11 * t18 + t12 * t9 + t28 * t34) + (-t12 + t109) * t68 + m(7) * (t1 * t120 + t108 * t11 + t109 * t9 - t121 * t2 + t124 * t13) + (-pkin(4) * t33 - t19 * t25 - t20 * t26 - t46 * t50) * m(6) - pkin(4) * t23 - t32 * mrSges(5,2) - t34 * t56 - t26 * t67 - t25 * t69 - t49 * t95 - Ifges(5,5) * t102 - Ifges(5,6) * t103 + t120 * t41 - t121 * t43 + t124 * t22 - ((t258 - Ifges(5,1) / 0.2e1) * t111 - t250) * t110 - t205 * t50; t141 + (t11 * t94 + t9 * t93) * mrSges(7,3) + (t19 * t93 + t20 * t94) * mrSges(6,3) + (-t56 * t94 + t41) * pkin(5) + t98 + t99 + t63 + t62 + t61 + t60 - t10 * t66 - t19 * t67 + t11 * t68 + t20 * t69 - t28 * (mrSges(7,1) * t94 + mrSges(7,2) * t93) - t46 * (mrSges(6,1) * t94 + mrSges(6,2) * t93) + (-(t10 - t9) * t11 + (-t28 * t94 + t1) * pkin(5)) * m(7) + (t257 * t93 - t251) * t224 + t245 * t223 + (-t246 * t94 + t247 * t93) * t221 + (-t255 * t94 + t237 + t254) * t226; -t93 * t66 + t94 * t68 + 0.2e1 * (t13 / 0.2e1 + t11 * t226 + t9 * t223) * m(7) + t22;];
tauc  = t14(:);
