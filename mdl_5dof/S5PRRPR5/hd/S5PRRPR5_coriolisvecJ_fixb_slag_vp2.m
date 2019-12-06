% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:19
% EndTime: 2019-12-05 16:25:39
% DurationCPUTime: 5.37s
% Computational Cost: add. (2741->368), mult. (7416->528), div. (0->0), fcn. (5299->10), ass. (0->192)
t125 = sin(qJ(2));
t121 = sin(pkin(5));
t181 = qJD(1) * t121;
t163 = t125 * t181;
t124 = sin(qJ(3));
t176 = qJD(3) * t124;
t246 = pkin(3) * t176 - t163;
t120 = sin(pkin(10));
t127 = cos(qJ(3));
t199 = -qJ(4) - pkin(7);
t156 = qJD(3) * t199;
t132 = -qJD(4) * t124 + t127 * t156;
t186 = t120 * t124;
t189 = cos(pkin(10));
t136 = t127 * t189 - t186;
t128 = cos(qJ(2));
t162 = t128 * t181;
t174 = qJD(4) * t127;
t88 = t124 * t156 + t174;
t245 = -t120 * t132 + t136 * t162 - t189 * t88;
t155 = t189 * t124;
t100 = t120 * t127 + t155;
t92 = t100 * qJD(3);
t93 = t136 * qJD(3);
t244 = pkin(4) * t92 - pkin(8) * t93 + t246;
t90 = t136 * qJD(2);
t243 = -t90 / 0.2e1;
t123 = sin(qJ(5));
t126 = cos(qJ(5));
t118 = -pkin(3) * t127 - pkin(2);
t59 = -pkin(4) * t136 - pkin(8) * t100 + t118;
t108 = t199 * t127;
t71 = -t108 * t189 + t186 * t199;
t24 = t123 * t59 + t126 * t71;
t242 = -qJD(5) * t24 + t123 * t245 + t126 * t244;
t23 = -t123 * t71 + t126 * t59;
t241 = qJD(5) * t23 + t123 * t244 - t126 * t245;
t240 = t100 * t162 - t120 * t88 + t132 * t189;
t177 = qJD(2) * t127;
t91 = -qJD(2) * t155 - t120 * t177;
t201 = t91 * mrSges(5,3);
t73 = qJD(3) * t126 + t123 * t91;
t74 = qJD(3) * t123 - t126 * t91;
t238 = qJD(3) * mrSges(5,1) + mrSges(6,1) * t73 - mrSges(6,2) * t74 + t201;
t237 = -Ifges(5,6) / 0.2e1;
t85 = qJD(2) * t118 + qJD(4) - t162;
t236 = t85 * mrSges(5,2);
t235 = Ifges(5,5) * qJD(3);
t171 = qJD(2) * qJD(3);
t158 = t124 * t171;
t234 = qJD(5) - t90;
t175 = qJD(3) * t127;
t104 = qJD(2) * pkin(7) + t163;
t122 = cos(pkin(5));
t180 = qJD(1) * t124;
t161 = t122 * t180;
t76 = t104 * t127 + t161;
t233 = -t76 * qJD(3) + (-qJ(4) * t175 + (-qJD(4) - t162) * t124) * qJD(2);
t183 = t122 * t127;
t113 = qJD(1) * t183;
t184 = t121 * t128;
t159 = qJD(2) * t184;
t153 = t127 * t159;
t55 = qJD(1) * t153 + qJD(3) * t113 - t104 * t176;
t39 = (-qJ(4) * t176 + t174) * qJD(2) + t55;
t14 = t120 * t233 + t189 * t39;
t83 = qJD(2) * t92;
t84 = qJD(2) * t93;
t178 = qJD(2) * t125;
t160 = t121 * t178;
t89 = pkin(3) * t158 + qJD(1) * t160;
t31 = pkin(4) * t83 - pkin(8) * t84 + t89;
t154 = qJ(4) * qJD(2) + t104;
t67 = t127 * t154 + t161;
t61 = t189 * t67;
t66 = -t124 * t154 + t113;
t63 = qJD(3) * pkin(3) + t66;
t26 = t120 * t63 + t61;
t22 = qJD(3) * pkin(8) + t26;
t34 = -pkin(4) * t90 + pkin(8) * t91 + t85;
t7 = -t123 * t22 + t126 * t34;
t1 = qJD(5) * t7 + t123 * t31 + t126 * t14;
t8 = t123 * t34 + t126 * t22;
t2 = -qJD(5) * t8 - t123 * t14 + t126 * t31;
t232 = t1 * t126 - t123 * t2;
t40 = qJD(5) * t73 + t126 * t84;
t41 = -qJD(5) * t74 - t123 * t84;
t231 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t40 + Ifges(6,6) * t41;
t230 = qJD(3) * t237;
t213 = Ifges(5,4) * t91;
t229 = Ifges(5,2) * t243 + t213 / 0.2e1;
t147 = mrSges(6,1) * t123 + mrSges(6,2) * t126;
t192 = t120 * t67;
t25 = t189 * t63 - t192;
t21 = -qJD(3) * pkin(4) - t25;
t138 = t21 * t147;
t142 = Ifges(6,5) * t126 - Ifges(6,6) * t123;
t194 = Ifges(6,4) * t126;
t144 = -Ifges(6,2) * t123 + t194;
t195 = Ifges(6,4) * t123;
t146 = Ifges(6,1) * t126 - t195;
t212 = Ifges(6,4) * t74;
t17 = Ifges(6,2) * t73 + Ifges(6,6) * t234 + t212;
t72 = Ifges(6,4) * t73;
t18 = Ifges(6,1) * t74 + Ifges(6,5) * t234 + t72;
t214 = t126 / 0.2e1;
t215 = -t123 / 0.2e1;
t220 = t74 / 0.2e1;
t228 = t234 * t142 / 0.2e1 + t146 * t220 + t73 * t144 / 0.2e1 + t138 + t18 * t214 + t17 * t215;
t227 = t85 * mrSges(5,1) + t7 * mrSges(6,1) - t8 * mrSges(6,2) + t229 + t230;
t129 = qJD(2) ^ 2;
t225 = t40 / 0.2e1;
t224 = t41 / 0.2e1;
t222 = -t73 / 0.2e1;
t221 = -t74 / 0.2e1;
t219 = t83 / 0.2e1;
t218 = -t234 / 0.2e1;
t86 = Ifges(5,4) * t90;
t211 = pkin(3) * t120;
t13 = t120 * t39 - t189 * t233;
t185 = t121 * t125;
t94 = -t124 * t185 + t183;
t95 = t122 * t124 + t127 * t185;
t51 = t120 * t95 - t189 * t94;
t208 = t13 * t51;
t70 = -t108 * t120 - t155 * t199;
t207 = t13 * t70;
t206 = t73 * Ifges(6,6);
t205 = t74 * Ifges(6,5);
t204 = t234 * Ifges(6,3);
t203 = t90 * mrSges(5,3);
t200 = t91 * Ifges(5,1);
t196 = Ifges(4,4) * t124;
t191 = t123 * t90;
t190 = t126 * t90;
t188 = Ifges(4,5) * qJD(3);
t187 = Ifges(4,6) * qJD(3);
t179 = qJD(2) * t124;
t173 = qJD(5) * t123;
t172 = qJD(5) * t126;
t170 = pkin(3) * t179;
t168 = mrSges(4,3) * t179;
t167 = mrSges(4,3) * t177;
t164 = t189 * pkin(3);
t46 = mrSges(5,1) * t83 + mrSges(5,2) * t84;
t151 = -t1 * t123 - t126 * t2;
t150 = -t123 * t8 - t126 * t7;
t149 = t123 * t7 - t126 * t8;
t148 = mrSges(6,1) * t126 - mrSges(6,2) * t123;
t145 = Ifges(6,1) * t123 + t194;
t143 = Ifges(6,2) * t126 + t195;
t141 = Ifges(6,5) * t123 + Ifges(6,6) * t126;
t56 = -t104 * t175 + (-qJD(3) * t122 - t159) * t180;
t140 = -t124 * t56 + t127 * t55;
t52 = t120 * t94 + t189 * t95;
t35 = -t123 * t52 - t126 * t184;
t139 = t123 * t184 - t126 * t52;
t137 = (mrSges(4,1) * t124 + mrSges(4,2) * t127) * qJD(2);
t119 = Ifges(4,4) * t177;
t117 = -t164 - pkin(4);
t116 = pkin(8) + t211;
t107 = -qJD(3) * mrSges(4,2) + t167;
t106 = qJD(3) * mrSges(4,1) - t168;
t105 = -qJD(2) * pkin(2) - t162;
t98 = qJD(3) * t137;
t97 = Ifges(4,1) * t179 + t119 + t188;
t96 = t187 + (t127 * Ifges(4,2) + t196) * qJD(2);
t81 = Ifges(6,3) * t83;
t77 = -qJD(3) * mrSges(5,2) + t203;
t75 = -t104 * t124 + t113;
t65 = qJD(3) * t94 + t153;
t64 = -qJD(3) * t95 - t124 * t159;
t57 = -mrSges(5,1) * t90 - mrSges(5,2) * t91;
t54 = -t200 + t86 + t235;
t47 = -pkin(4) * t91 - pkin(8) * t90 + t170;
t43 = mrSges(6,1) * t234 - mrSges(6,3) * t74;
t42 = -mrSges(6,2) * t234 + mrSges(6,3) * t73;
t30 = t189 * t66 - t192;
t29 = t120 * t64 + t189 * t65;
t28 = t120 * t66 + t61;
t27 = t120 * t65 - t189 * t64;
t20 = -mrSges(6,2) * t83 + mrSges(6,3) * t41;
t19 = mrSges(6,1) * t83 - mrSges(6,3) * t40;
t16 = t204 + t205 + t206;
t15 = -mrSges(6,1) * t41 + mrSges(6,2) * t40;
t12 = t123 * t47 + t126 * t30;
t11 = -t123 * t30 + t126 * t47;
t10 = Ifges(6,1) * t40 + Ifges(6,4) * t41 + Ifges(6,5) * t83;
t9 = Ifges(6,4) * t40 + Ifges(6,2) * t41 + Ifges(6,6) * t83;
t6 = qJD(5) * t139 - t123 * t29 + t126 * t160;
t5 = qJD(5) * t35 + t123 * t160 + t126 * t29;
t3 = [t64 * t106 + t65 * t107 + t51 * t15 + t35 * t19 - t139 * t20 + t29 * t77 + t5 * t42 + t6 * t43 - t238 * t27 + (t51 * t84 - t52 * t83) * mrSges(5,3) + (-t124 * t95 - t127 * t94) * mrSges(4,3) * t171 + ((-mrSges(3,2) * t129 - t46 - t98) * t128 + (-mrSges(3,1) * t129 + (qJD(2) * (-mrSges(4,1) * t127 + mrSges(4,2) * t124) + t57) * qJD(2)) * t125) * t121 + m(4) * (t55 * t95 + t56 * t94 + t64 * t75 + t65 * t76 + (t105 - t162) * t160) + m(6) * (-t1 * t139 + t2 * t35 + t21 * t27 + t5 * t8 + t6 * t7 + t208) + m(5) * (t208 + t14 * t52 - t25 * t27 + t26 * t29 + (-t128 * t89 + t178 * t85) * t121); (-t200 / 0.2e1 + t86 / 0.2e1 + t236 + t54 / 0.2e1 + t150 * mrSges(6,3) + t228) * t93 + t118 * t46 - pkin(2) * t98 + t70 * t15 + t23 * t19 + t24 * t20 + (-0.3e1 / 0.2e1 * t124 ^ 2 + 0.3e1 / 0.2e1 * t127 ^ 2) * Ifges(4,4) * t171 + (t16 / 0.2e1 + t204 / 0.2e1 + t206 / 0.2e1 + t205 / 0.2e1 + t227 + t229) * t92 - (-t14 * mrSges(5,3) + t81 / 0.2e1 - Ifges(5,4) * t84 + t89 * mrSges(5,1) + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t83 + t231) * t136 + ((t106 * t124 - t107 * t127) * t128 - t125 * t57) * t181 + (Ifges(5,5) * t93 / 0.2e1 + t92 * t237 + (t105 * mrSges(4,2) + t97 / 0.2e1 - pkin(7) * t106 - t75 * mrSges(4,3) + t188 / 0.2e1) * t127 + (t105 * mrSges(4,1) - t96 / 0.2e1 - pkin(7) * t107 - t76 * mrSges(4,3) + pkin(3) * t57 - t187 / 0.2e1 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t177) * t124) * qJD(3) + (-t25 * t93 - t26 * t92 + t70 * t84 - t71 * t83) * mrSges(5,3) + t242 * t43 - t245 * t77 + t241 * t42 + (t10 * t214 + t9 * t215 + t142 * t219 + t146 * t225 + t144 * t224 + t89 * mrSges(5,2) - Ifges(5,4) * t83 + Ifges(5,1) * t84 + (mrSges(5,3) + t147) * t13 + t151 * mrSges(6,3) + (t21 * t148 + t141 * t218 + t143 * t222 + t145 * t221 - t126 * t17 / 0.2e1 + t18 * t215 + t149 * mrSges(6,3)) * qJD(5)) * t100 + t140 * mrSges(4,3) + t238 * t240 + (t1 * t24 + t2 * t23 - t21 * t240 + t241 * t8 + t242 * t7 + t207) * m(6) + (((-t124 * t76 - t127 * t75) * qJD(3) + t140) * pkin(7) - (pkin(2) * t178 + t105 * t125 + (-t124 * t75 + t127 * t76) * t128) * t181) * m(4) + (t118 * t89 + t14 * t71 + t240 * t25 - t245 * t26 + t246 * t85 + t207) * m(5); (m(6) * t117 - mrSges(5,1) - t148) * t13 + (-t236 - t235 / 0.2e1 + t142 * t218 + t146 * t221 + t144 * t222 - t138) * t90 + ((-t173 + t191) * t8 + (-t172 + t190) * t7 + t232) * mrSges(6,3) + (m(6) * t232 - t123 * t19 + t126 * t20 - t172 * t43 - t173 * t42) * t116 + (-m(6) * t21 + t238) * t28 + (m(6) * t116 * t150 + t228) * qJD(5) - (-Ifges(4,2) * t179 + t119 + t97) * t177 / 0.2e1 + (Ifges(5,1) * t90 + t16 + t213) * t91 / 0.2e1 + t123 * t10 / 0.2e1 + t117 * t15 - Ifges(5,6) * t83 + Ifges(5,5) * t84 - t30 * t77 + t56 * mrSges(4,1) - t55 * mrSges(4,2) - t12 * t42 - t11 * t43 - t14 * mrSges(5,2) + (t106 + t168) * t76 + (-Ifges(6,5) * t221 - Ifges(6,6) * t222 - Ifges(6,3) * t218 + t227 + t230) * t91 + (-t107 + t167) * t75 + t143 * t224 + t145 * t225 + t9 * t214 + t141 * t219 + t25 * t203 + (-t164 * t84 - t211 * t83) * mrSges(5,3) + ((t120 * t14 - t13 * t189) * pkin(3) - t170 * t85 + t25 * t28 - t26 * t30) * m(5) - m(6) * (t11 * t7 + t12 * t8) - t26 * t201 + t17 * t191 / 0.2e1 - t18 * t190 / 0.2e1 + t96 * t179 / 0.2e1 - t105 * t137 - t124 * t129 * (Ifges(4,1) * t127 - t196) / 0.2e1 + (Ifges(5,2) * t91 + t54 + t86) * t243 + t171 * Ifges(4,5) * t127 / 0.2e1 - Ifges(4,6) * t158 / 0.2e1 - t57 * t170; -t90 * t77 - t238 * t91 + (t234 * t42 + t19) * t126 + (-t234 * t43 + t20) * t123 + t46 + (-t149 * t234 + t21 * t91 - t151) * m(6) + (-t25 * t91 - t26 * t90 + t89) * m(5); t81 - t21 * (mrSges(6,1) * t74 + mrSges(6,2) * t73) + (Ifges(6,1) * t73 - t212) * t221 + t17 * t220 + (Ifges(6,5) * t73 - Ifges(6,6) * t74) * t218 - t7 * t42 + t8 * t43 + (t7 * t73 + t74 * t8) * mrSges(6,3) + (-Ifges(6,2) * t74 + t18 + t72) * t222 + t231;];
tauc = t3(:);
