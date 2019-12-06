% Calculate vector of inverse dynamics joint torques for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:11
% EndTime: 2019-12-05 15:14:19
% DurationCPUTime: 3.52s
% Computational Cost: add. (1898->293), mult. (4288->416), div. (0->0), fcn. (3092->14), ass. (0->149)
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t139 = -t117 * mrSges(5,1) + t114 * mrSges(5,2);
t212 = -mrSges(4,1) + t139;
t108 = qJ(4) + qJ(5);
t103 = sin(t108);
t104 = cos(t108);
t95 = pkin(4) * t117 + pkin(3);
t211 = m(5) * pkin(3) + m(6) * t95 + mrSges(6,1) * t104 - mrSges(6,2) * t103 - t212;
t119 = -pkin(7) - pkin(6);
t210 = -m(5) * pkin(6) + m(6) * t119 + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t115 = sin(qJ(3));
t118 = cos(qJ(3));
t157 = qJD(1) * qJD(3);
t209 = qJDD(1) * t115 + t118 * t157;
t208 = qJDD(1) * t118 - t115 * t157;
t109 = sin(pkin(9));
t111 = cos(pkin(9));
t207 = -t109 * t115 + t118 * t111;
t156 = qJD(3) * qJD(4);
t80 = qJDD(3) * t117 - t114 * t156;
t206 = t80 / 0.2e1;
t113 = sin(qJ(5));
t116 = cos(qJ(5));
t84 = t119 * t114;
t85 = t119 * t117;
t52 = t113 * t84 - t116 * t85;
t65 = t207 * qJD(1);
t76 = t113 * t117 + t114 * t116;
t147 = qJD(4) * t119;
t78 = t114 * t147;
t79 = t117 * t147;
t205 = -qJD(5) * t52 - t113 * t78 + t116 * t79 + t76 * t65;
t132 = t113 * t114 - t116 * t117;
t51 = t113 * t85 + t116 * t84;
t204 = qJD(5) * t51 + t113 * t79 + t116 * t78 + t132 * t65;
t74 = t109 * t118 + t111 * t115;
t36 = t132 * t74;
t197 = m(6) * pkin(4);
t203 = -mrSges(5,1) - t197;
t161 = qJD(3) * t114;
t149 = mrSges(5,3) * t161;
t82 = qJD(4) * mrSges(5,1) - t149;
t160 = qJD(3) * t117;
t148 = mrSges(5,3) * t160;
t83 = -qJD(4) * mrSges(5,2) + t148;
t133 = -t114 * t82 + t117 * t83;
t202 = -qJD(3) * mrSges(4,2) + t133;
t62 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t80;
t81 = qJDD(3) * t114 + t117 * t156;
t63 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t81;
t201 = -t114 * t63 + t117 * t62;
t102 = t117 * qJD(2);
t159 = qJD(4) * t114;
t42 = t208 * t109 + t209 * t111;
t37 = qJDD(3) * pkin(6) + t42;
t66 = t74 * qJD(1);
t61 = qJD(3) * pkin(6) + t66;
t12 = qJD(4) * t102 + t114 * qJDD(2) + t117 * t37 - t159 * t61;
t162 = qJD(2) * t114;
t49 = t117 * t61 + t162;
t13 = -qJD(4) * t49 + t117 * qJDD(2) - t114 * t37;
t200 = -t114 * t13 + t117 * t12;
t107 = qJD(4) + qJD(5);
t71 = t132 * qJD(3);
t72 = t76 * qJD(3);
t41 = mrSges(6,1) * t71 + mrSges(6,2) * t72;
t199 = t212 * qJD(3) + t41;
t195 = t72 / 0.2e1;
t194 = m(3) + m(4);
t106 = pkin(9) + qJ(3);
t98 = sin(t106);
t193 = g(3) * t98;
t142 = pkin(7) * qJD(3) + t61;
t45 = t117 * t142 + t162;
t182 = t113 * t45;
t44 = -t114 * t142 + t102;
t39 = qJD(4) * pkin(4) + t44;
t10 = t116 * t39 - t182;
t190 = t10 * mrSges(6,3);
t189 = t72 * mrSges(6,3);
t188 = t72 * Ifges(6,4);
t112 = cos(pkin(8));
t169 = t104 * t112;
t110 = sin(pkin(8));
t170 = t104 * t110;
t171 = t103 * t112;
t172 = t103 * t110;
t99 = cos(t106);
t187 = (-t172 * t99 - t169) * mrSges(6,1) + (-t170 * t99 + t171) * mrSges(6,2);
t186 = (-t171 * t99 + t170) * mrSges(6,1) + (-t169 * t99 - t172) * mrSges(6,2);
t185 = Ifges(5,4) * t114;
t184 = Ifges(5,4) * t117;
t178 = t116 * t45;
t167 = t110 * t114;
t166 = t110 * t117;
t165 = t112 * t114;
t164 = t112 * t117;
t158 = qJD(4) * t117;
t153 = pkin(4) * t161;
t152 = pkin(4) * t159;
t151 = m(5) + m(6) + t194;
t138 = mrSges(5,1) * t114 + mrSges(5,2) * t117;
t137 = -mrSges(6,1) * t103 - mrSges(6,2) * t104;
t136 = t117 * Ifges(5,2) + t185;
t135 = Ifges(5,5) * t117 - Ifges(5,6) * t114;
t11 = t113 * t39 + t178;
t48 = -t114 * t61 + t102;
t134 = -t114 * t48 + t117 * t49;
t60 = -qJD(3) * pkin(3) - t65;
t129 = t60 * t138;
t128 = t114 * (Ifges(5,1) * t117 - t185);
t43 = -t209 * t109 + t208 * t111;
t127 = t76 * qJD(5);
t126 = t132 * qJD(5);
t38 = -qJDD(3) * pkin(3) - t43;
t123 = (-t114 * t49 - t117 * t48) * qJD(4) + t200;
t105 = qJDD(4) + qJDD(5);
t8 = qJDD(4) * pkin(4) - pkin(7) * t81 + t13;
t9 = pkin(7) * t80 + t12;
t2 = qJD(5) * t10 + t113 * t8 + t116 * t9;
t27 = -qJD(3) * t126 + t113 * t80 + t116 * t81;
t28 = -qJD(3) * t127 - t113 * t81 + t116 * t80;
t3 = -qJD(5) * t11 - t113 * t9 + t116 * t8;
t33 = -t71 * Ifges(6,2) + t107 * Ifges(6,6) + t188;
t64 = Ifges(6,4) * t71;
t34 = t72 * Ifges(6,1) + t107 * Ifges(6,5) - t64;
t53 = -qJD(3) * t95 - t65;
t122 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - t71 * t190 + t33 * t195 - t53 * (mrSges(6,1) * t72 - mrSges(6,2) * t71) - t72 * (-Ifges(6,1) * t71 - t188) / 0.2e1 + Ifges(6,6) * t28 + Ifges(6,5) * t27 - t107 * (-Ifges(6,5) * t71 - Ifges(6,6) * t72) / 0.2e1 + Ifges(6,3) * t105 + (-Ifges(6,2) * t72 + t34 - t64) * t71 / 0.2e1;
t47 = -qJD(4) * t76 - t127;
t96 = Ifges(5,4) * t160;
t70 = Ifges(5,1) * t161 + Ifges(5,5) * qJD(4) + t96;
t69 = Ifges(5,6) * qJD(4) + qJD(3) * t136;
t68 = t74 * qJD(3);
t67 = t207 * qJD(3);
t59 = mrSges(6,1) * t107 - t189;
t58 = -mrSges(6,2) * t107 - mrSges(6,3) * t71;
t50 = -mrSges(5,1) * t80 + mrSges(5,2) * t81;
t46 = -qJD(4) * t132 - t126;
t35 = t76 * t74;
t24 = -pkin(4) * t80 + t38;
t21 = -mrSges(6,2) * t105 + mrSges(6,3) * t28;
t20 = mrSges(6,1) * t105 - mrSges(6,3) * t27;
t15 = t116 * t44 - t182;
t14 = -t113 * t44 - t178;
t6 = -mrSges(6,1) * t28 + mrSges(6,2) * t27;
t5 = t107 * t36 - t76 * t67;
t4 = -t132 * t67 + t47 * t74;
t1 = [-t35 * t20 - t36 * t21 + t4 * t58 + t5 * t59 - (-qJDD(3) * mrSges(4,1) + t50 + t6) * t207 + t199 * t68 + t202 * t67 + (-qJDD(3) * mrSges(4,2) + (-t114 * t83 - t117 * t82) * qJD(4) + t201) * t74 + (-m(2) - t151) * g(3) + m(5) * (t123 * t74 + t134 * t67 - t207 * t38 + t60 * t68) + m(6) * (t10 * t5 + t11 * t4 - t2 * t36 - t207 * t24 - t3 * t35 + t53 * t68) + m(4) * (t207 * t43 + t42 * t74 - t65 * t68 + t66 * t67) + (m(2) + m(3) * (t109 ^ 2 + t111 ^ 2)) * qJDD(1); t114 * t62 + t117 * t63 - t132 * t20 + t76 * t21 + t46 * t58 + t47 * t59 + t194 * qJDD(2) + t133 * qJD(4) + m(5) * (qJD(4) * t134 + t114 * t12 + t117 * t13) + m(6) * (t10 * t47 + t11 * t46 - t132 * t3 + t2 * t76) + (-g(1) * t110 + g(2) * t112) * t151; (g(1) * t112 + g(2) * t110) * (t210 * t99 + t211 * t98) + (t210 * t98 - t211 * t99) * g(3) + t38 * t139 + t117 * (Ifges(5,4) * t81 + Ifges(5,2) * t80) / 0.2e1 + qJDD(4) * (Ifges(5,5) * t114 + Ifges(5,6) * t117) + t107 * (Ifges(6,5) * t46 + Ifges(6,6) * t47) / 0.2e1 - t95 * t6 + t136 * t206 + t11 * t47 * mrSges(6,3) - t71 * (Ifges(6,4) * t46 + Ifges(6,2) * t47) / 0.2e1 + t46 * t34 / 0.2e1 + t47 * t33 / 0.2e1 + t51 * t20 + t52 * t21 + t53 * (-mrSges(6,1) * t47 + mrSges(6,2) * t46) - t42 * mrSges(4,2) + t43 * mrSges(4,1) + (t129 + t135 * qJD(4) / 0.2e1) * qJD(4) + (t117 * (-Ifges(5,2) * t114 + t184) + t128) * t156 / 0.2e1 + t41 * t152 + t81 * t184 / 0.2e1 - (-mrSges(6,1) * t24 + mrSges(6,3) * t2 + Ifges(6,4) * t27 + Ifges(6,2) * t28 + Ifges(6,6) * t105) * t132 + Ifges(4,3) * qJDD(3) + (Ifges(6,1) * t46 + Ifges(6,4) * t47) * t195 + t70 * t158 / 0.2e1 - t69 * t159 / 0.2e1 - t46 * t190 + (-m(5) * t60 - m(6) * t53 - t199) * t66 + (-t158 * t48 - t159 * t49 + t200) * mrSges(5,3) + (m(5) * t123 - t158 * t82 - t159 * t83 + t201) * pkin(6) + (-m(5) * t134 - t202) * t65 + t204 * t58 + t205 * t59 + (t10 * t205 + t11 * t204 + t152 * t53 + t2 * t52 - t24 * t95 + t3 * t51) * m(6) + (Ifges(5,1) * t81 + Ifges(5,4) * t206) * t114 + (mrSges(6,2) * t24 - mrSges(6,3) * t3 + Ifges(6,1) * t27 + Ifges(6,4) * t28 + Ifges(6,5) * t105) * t76 + (-m(5) * t38 - t50) * pkin(3); Ifges(5,6) * t80 + Ifges(5,5) * t81 + t122 - t15 * t58 - t14 * t59 - t12 * mrSges(5,2) + t13 * mrSges(5,1) + t11 * t189 + (t113 * t2 + t116 * t3 + (-t10 * t113 + t11 * t116) * qJD(5)) * t197 - t41 * t153 - m(6) * (t10 * t14 + t11 * t15 + t153 * t53) - t135 * t156 / 0.2e1 + t69 * t161 / 0.2e1 + Ifges(5,3) * qJDD(4) + (t82 + t149) * t49 + (-t83 + t148) * t48 + (t114 * t197 - t137 + t138) * t193 - (-Ifges(5,2) * t161 + t70 + t96) * t160 / 0.2e1 + (-t129 - t128 * qJD(3) / 0.2e1) * qJD(3) + (-(-t166 * t99 + t165) * mrSges(5,2) - t187 + t203 * (-t167 * t99 - t164)) * g(2) + (-(-t164 * t99 - t167) * mrSges(5,2) - t186 + t203 * (-t165 * t99 + t166)) * g(1) + ((-t113 * t59 + t116 * t58) * qJD(5) + t113 * t21 + t116 * t20) * pkin(4); (t59 + t189) * t11 - t137 * t193 + t122 - g(1) * t186 - t10 * t58 - g(2) * t187;];
tau = t1;
