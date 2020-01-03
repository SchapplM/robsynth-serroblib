% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR14_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR14_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:31
% EndTime: 2019-12-31 18:34:39
% DurationCPUTime: 3.56s
% Computational Cost: add. (2496->334), mult. (5652->479), div. (0->0), fcn. (3520->6), ass. (0->166)
t210 = qJ(2) * (m(4) + m(3)) + mrSges(3,3);
t101 = cos(qJ(5));
t157 = cos(pkin(8));
t181 = cos(qJ(3));
t123 = t157 * t181;
t117 = qJD(1) * t123;
t100 = sin(qJ(3));
t153 = qJD(1) * t100;
t98 = sin(pkin(8));
t67 = t153 * t98 - t117;
t99 = sin(qJ(5));
t49 = qJD(3) * t101 + t67 * t99;
t208 = t49 * Ifges(6,6);
t109 = -t100 * t157 - t181 * t98;
t68 = t109 * qJD(1);
t202 = qJD(5) - t68;
t207 = t202 * Ifges(6,3);
t206 = t68 * Ifges(5,2);
t171 = t67 * mrSges(5,3);
t50 = qJD(3) * t99 - t101 * t67;
t164 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t49 + mrSges(6,2) * t50 - t171;
t61 = qJD(3) * t68;
t172 = t61 * mrSges(5,3);
t27 = qJD(5) * t49 + t101 * t61;
t28 = -qJD(5) * t50 - t61 * t99;
t7 = -mrSges(6,1) * t28 + mrSges(6,2) * t27;
t205 = t7 + t172;
t204 = Ifges(5,5) * qJD(3);
t203 = Ifges(5,6) * qJD(3);
t102 = -pkin(1) - pkin(6);
t155 = qJ(4) - t102;
t151 = qJD(5) * t101;
t156 = qJD(5) * t99;
t84 = qJD(1) * t102 + qJD(2);
t63 = (-qJ(4) * qJD(1) + t84) * t100;
t55 = t157 * t63;
t137 = qJD(1) * t181;
t77 = t181 * t84;
t64 = -qJ(4) * t137 + t77;
t57 = qJD(3) * pkin(3) + t64;
t30 = t98 * t57 + t55;
t26 = qJD(3) * pkin(7) + t30;
t79 = pkin(3) * t153 + qJD(1) * qJ(2) + qJD(4);
t33 = -pkin(4) * t68 + pkin(7) * t67 + t79;
t8 = t101 * t33 - t26 * t99;
t9 = t101 * t26 + t33 * t99;
t201 = -t8 * t151 - t9 * t156;
t149 = qJD(1) * qJD(3);
t134 = t100 * t149;
t60 = -qJD(3) * t117 + t134 * t98;
t17 = -mrSges(6,1) * t60 - mrSges(6,3) * t27;
t18 = mrSges(6,2) * t60 + mrSges(6,3) * t28;
t200 = t101 * t18 - t99 * t17;
t199 = t9 * t101 - t8 * t99;
t135 = t181 * qJD(4);
t152 = qJD(3) * t100;
t141 = t84 * t152;
t104 = -t141 + (qJ(4) * t152 - t135) * qJD(1);
t136 = qJD(3) * t181;
t129 = t84 * t136;
t150 = t100 * qJD(4);
t47 = t129 + (-qJ(4) * t136 - t150) * qJD(1);
t19 = -t104 * t157 + t47 * t98;
t74 = t100 * t98 - t123;
t176 = t19 * t74;
t168 = t98 * t63;
t29 = t157 * t57 - t168;
t69 = -qJD(3) * t123 + t152 * t98;
t70 = t109 * qJD(3);
t198 = t29 * t70 - t30 * t69 + t176;
t25 = -qJD(3) * pkin(4) - t29;
t197 = -m(6) * t25 - t164;
t20 = t104 * t98 + t157 * t47;
t128 = qJD(1) * t136;
t96 = qJD(1) * qJD(2);
t78 = pkin(3) * t128 + t96;
t22 = -pkin(4) * t60 - pkin(7) * t61 + t78;
t1 = qJD(5) * t8 + t101 * t20 + t22 * t99;
t2 = -qJD(5) * t9 + t101 * t22 - t20 * t99;
t196 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t195 = qJD(1) ^ 2;
t194 = t27 / 0.2e1;
t193 = t28 / 0.2e1;
t192 = -t49 / 0.2e1;
t191 = -t50 / 0.2e1;
t190 = t50 / 0.2e1;
t189 = -t60 / 0.2e1;
t188 = -t202 / 0.2e1;
t187 = -t67 / 0.2e1;
t184 = t99 / 0.2e1;
t183 = pkin(3) * t98;
t180 = Ifges(5,4) * t67;
t179 = Ifges(6,4) * t50;
t178 = Ifges(6,4) * t99;
t112 = t155 * t181;
t80 = t155 * t100;
t45 = t112 * t157 - t80 * t98;
t177 = t19 * t45;
t173 = t60 * mrSges(5,3);
t170 = t68 * mrSges(5,3);
t167 = t99 * mrSges(6,3);
t163 = Ifges(4,4) * t100;
t162 = Ifges(6,4) * t101;
t82 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t137;
t161 = t100 * t82;
t160 = t101 * mrSges(6,3);
t158 = t101 * t74;
t92 = t100 * pkin(3) + qJ(2);
t85 = pkin(3) * t136 + qJD(2);
t148 = Ifges(6,5) * t27 + Ifges(6,6) * t28 - Ifges(6,3) * t60;
t15 = Ifges(6,2) * t49 + Ifges(6,6) * t202 + t179;
t143 = t15 * t184;
t142 = Ifges(4,4) * t181;
t48 = Ifges(6,4) * t49;
t16 = Ifges(6,1) * t50 + Ifges(6,5) * t202 + t48;
t140 = -t101 * t16 / 0.2e1;
t139 = t157 * pkin(3);
t138 = -t60 * mrSges(5,1) + t61 * mrSges(5,2);
t132 = t151 / 0.2e1;
t130 = pkin(3) * t137;
t126 = t1 * t101 - t2 * t99;
t124 = t8 * t101 + t9 * t99;
t122 = mrSges(6,1) * t99 + mrSges(6,2) * t101;
t121 = Ifges(6,1) * t101 - t178;
t120 = -Ifges(6,2) * t99 + t162;
t119 = Ifges(6,5) * t101 - Ifges(6,6) * t99;
t31 = -mrSges(6,2) * t202 + mrSges(6,3) * t49;
t32 = mrSges(6,1) * t202 - mrSges(6,3) * t50;
t118 = -t101 * t32 - t99 * t31;
t43 = -pkin(4) * t109 + pkin(7) * t74 + t92;
t46 = -t112 * t98 - t157 * t80;
t13 = t101 * t46 + t43 * t99;
t12 = t101 * t43 - t46 * t99;
t116 = t101 * t70 + t156 * t74;
t115 = t151 * t74 - t70 * t99;
t114 = -Ifges(4,5) * t100 - Ifges(4,6) * t181;
t111 = qJ(2) * (mrSges(4,1) * t181 - mrSges(4,2) * t100);
t110 = t100 * (-Ifges(4,2) * t181 - t163);
t108 = (Ifges(4,1) * t181 - t163) * qJD(1);
t107 = (-Ifges(4,2) * t100 + t142) * qJD(1);
t76 = (t100 * mrSges(4,1) + mrSges(4,2) * t181) * qJD(1);
t106 = t181 * (-Ifges(4,1) * t100 - t142);
t105 = t152 * t155 - t135;
t91 = -t139 - pkin(4);
t81 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t153;
t72 = Ifges(4,5) * qJD(3) + t108;
t71 = Ifges(4,6) * qJD(3) + t107;
t65 = Ifges(5,4) * t68;
t62 = -qJD(3) * t112 - t150;
t53 = -qJD(3) * mrSges(5,2) + t170;
t42 = -mrSges(5,1) * t68 - mrSges(5,2) * t67;
t41 = -t67 * Ifges(5,1) + t204 + t65;
t40 = -t180 + t203 + t206;
t39 = -t67 * pkin(4) - t68 * pkin(7) + t130;
t38 = -pkin(4) * t69 - pkin(7) * t70 + t85;
t37 = t157 * t64 - t168;
t36 = t64 * t98 + t55;
t35 = t105 * t98 + t157 * t62;
t14 = t50 * Ifges(6,5) + t207 + t208;
t11 = t101 * t37 + t39 * t99;
t10 = t101 * t39 - t37 * t99;
t6 = t27 * Ifges(6,1) + t28 * Ifges(6,4) - t60 * Ifges(6,5);
t5 = t27 * Ifges(6,4) + t28 * Ifges(6,2) - t60 * Ifges(6,6);
t4 = -qJD(5) * t13 + t101 * t38 - t35 * t99;
t3 = qJD(5) * t12 + t101 * t35 + t38 * t99;
t21 = [(t115 * t9 - t116 * t8 + t158 * t2) * mrSges(6,3) + m(5) * (t20 * t46 + t30 * t35 + t78 * t92 + t79 * t85 + t177) + m(6) * (t1 * t13 + t12 * t2 + t3 * t9 + t4 * t8 + t177) - (t107 + t71) * t136 / 0.2e1 - (t108 + t72) * t152 / 0.2e1 + (0.2e1 * t111 - t110 + t106) * t149 + t85 * t42 + 0.2e1 * qJD(2) * t76 + t35 * t53 + t3 * t31 + t4 * t32 + t12 * t17 + t13 * t18 + t202 * (Ifges(6,5) * t116 + Ifges(6,6) * t115) / 0.2e1 + (t61 * t109 + t68 * t70 / 0.2e1 + t69 * t187) * Ifges(5,4) + (t109 * t20 - t198) * mrSges(5,3) - (t78 * mrSges(5,1) + t148 / 0.2e1 - Ifges(5,2) * t60 + Ifges(6,3) * t189 + Ifges(6,6) * t193 + Ifges(6,5) * t194 + t196) * t109 + (-Ifges(5,4) * t60 + (qJD(5) * t16 + t5) * t184 - t78 * mrSges(5,2) - Ifges(5,1) * t61 + t1 * t167 - t119 * t189 - t120 * t193 - t121 * t194 + t132 * t15) * t74 - (t143 + t140 - t79 * mrSges(5,2) - t204 / 0.2e1 - t41 / 0.2e1 - Ifges(5,1) * t187) * t70 + t25 * (-mrSges(6,1) * t115 + mrSges(6,2) * t116) + qJD(3) ^ 2 * t114 / 0.2e1 + 0.2e1 * t210 * t96 - t6 * t158 / 0.2e1 + t92 * t138 - t122 * t176 + (-m(5) * t29 - t197) * (-t105 * t157 + t62 * t98) + t205 * t45 + (-t79 * mrSges(5,1) + t203 / 0.2e1 + t40 / 0.2e1 - t14 / 0.2e1 - t208 / 0.2e1 - t207 / 0.2e1 - t8 * mrSges(6,1) + t9 * mrSges(6,2) + t206 / 0.2e1 - Ifges(6,5) * t190) * t69 + t46 * t173 + (Ifges(6,1) * t116 + Ifges(6,4) * t115) * t190 + t49 * (Ifges(6,4) * t116 + Ifges(6,2) * t115) / 0.2e1 + (t136 * t81 - t152 * t82) * t102; t205 * t74 - t164 * t70 + (t181 * t81 - t161) * qJD(3) + m(6) * (-t25 * t70 + t176) + m(5) * t198 - (t173 + t118 * qJD(5) + m(6) * (t126 + t201) + m(5) * t20 + t200) * t109 - t210 * t195 + (-m(5) * t79 - m(6) * t124 + t118 - t42 - t76) * qJD(1) + (-m(6) * t199 - t101 * t31 + t99 * t32 - t53) * t69; t202 * t25 * t122 + t68 * t143 + t68 * t140 + t16 * t132 - (Ifges(5,2) * t67 + t41 + t65) * t68 / 0.2e1 + (Ifges(5,1) * t68 + t14 + t180) * t67 / 0.2e1 + (-t106 / 0.2e1 + t110 / 0.2e1 - t111) * t195 + t101 * t5 / 0.2e1 + t91 * t7 - t79 * (-mrSges(5,1) * t67 + mrSges(5,2) * t68) + (m(6) * t91 - mrSges(6,1) * t101 + mrSges(6,2) * t99 - mrSges(5,1)) * t19 - Ifges(4,5) * t134 - t42 * t130 - mrSges(4,2) * t129 - qJD(3) * (Ifges(5,5) * t68 + Ifges(5,6) * t67) / 0.2e1 + Ifges(5,6) * t60 + Ifges(5,5) * t61 - t37 * t53 - t11 * t31 - t10 * t32 - t20 * mrSges(5,2) + (t119 * t202 + t120 * t49 + t121 * t50) * qJD(5) / 0.2e1 - m(6) * (t10 * t8 + t11 * t9) - Ifges(4,6) * t128 + (-t31 * t156 - t32 * t151 + m(6) * (-qJD(5) * t124 + t126) + t200) * (pkin(7) + t183) + t201 * mrSges(6,3) - t8 * (-mrSges(6,1) * t67 - t160 * t68) - t15 * t156 / 0.2e1 - t114 * t149 / 0.2e1 + t72 * t153 / 0.2e1 - mrSges(4,1) * t141 - t81 * t77 + t71 * t137 / 0.2e1 - t30 * t171 - t139 * t172 - t9 * (mrSges(6,2) * t67 - t167 * t68) - t2 * t167 + t197 * t36 + t1 * t160 + t84 * t161 + t29 * t170 + t173 * t183 + t6 * t184 + t40 * t187 + (-Ifges(6,3) * t67 + t119 * t68) * t188 + (Ifges(6,5) * t99 + Ifges(6,6) * t101) * t189 + (-Ifges(6,5) * t67 + t121 * t68) * t191 + (-Ifges(6,6) * t67 + t120 * t68) * t192 + (Ifges(6,2) * t101 + t178) * t193 + (Ifges(6,1) * t99 + t162) * t194 + ((-t157 * t19 + t20 * t98) * pkin(3) - t130 * t79 + t29 * t36 - t30 * t37) * m(5); -t68 * t53 + t164 * t67 + (-t202 * t32 + t18) * t99 + (t202 * t31 + t17) * t101 + t138 + (t1 * t99 + t2 * t101 + t199 * t202 + t25 * t67) * m(6) + (-t29 * t67 - t30 * t68 + t78) * m(5); -t25 * (mrSges(6,1) * t50 + mrSges(6,2) * t49) + (Ifges(6,1) * t49 - t179) * t191 + t15 * t190 + (Ifges(6,5) * t49 - Ifges(6,6) * t50) * t188 - t8 * t31 + t9 * t32 + (t49 * t8 + t50 * t9) * mrSges(6,3) + t148 + (-Ifges(6,2) * t50 + t16 + t48) * t192 + t196;];
tauc = t21(:);
