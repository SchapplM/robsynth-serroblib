% Calculate time derivative of joint inertia matrix for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:12:14
% EndTime: 2019-03-09 05:12:20
% DurationCPUTime: 2.69s
% Computational Cost: add. (5056->260), mult. (10888->376), div. (0->0), fcn. (11273->8), ass. (0->124)
t209 = -2 * mrSges(4,3);
t208 = -2 * Ifges(4,4);
t207 = Ifges(4,1) - Ifges(4,2);
t206 = mrSges(6,2) - mrSges(5,1);
t99 = cos(qJ(6));
t161 = qJD(6) * t99;
t183 = sin(qJ(3));
t185 = cos(qJ(3));
t95 = sin(pkin(10));
t96 = cos(pkin(10));
t118 = t183 * t95 - t185 * t96;
t113 = qJD(3) * t118;
t79 = t183 * t96 + t185 * t95;
t114 = qJD(3) * t79;
t98 = sin(qJ(4));
t115 = t98 * t118;
t184 = cos(qJ(4));
t141 = qJD(4) * t184;
t47 = -qJD(4) * t115 - t113 * t98 + t114 * t184 + t141 * t79;
t110 = t184 * t118;
t64 = t79 * t98 + t110;
t97 = sin(qJ(6));
t120 = t64 * t161 + t47 * t97;
t162 = qJD(6) * t97;
t119 = t64 * t162 - t47 * t99;
t180 = Ifges(7,4) * t99;
t84 = -Ifges(7,2) * t97 + t180;
t181 = Ifges(7,4) * t97;
t85 = Ifges(7,1) * t99 - t181;
t205 = t84 * t99 + t85 * t97;
t166 = t97 ^ 2 + t99 ^ 2;
t204 = (-mrSges(5,2) * t184 + (-t166 * mrSges(7,3) + t206) * t98) * pkin(3) * qJD(4);
t167 = pkin(7) + qJ(2);
t199 = t167 * t185;
t203 = -t183 * qJD(2) - qJD(3) * t199;
t200 = t167 * t183;
t202 = -t185 * qJD(2) + qJD(3) * t200;
t83 = mrSges(7,1) * t97 + mrSges(7,2) * t99;
t201 = mrSges(6,3) + t83;
t65 = t184 * t79 - t115;
t147 = -t96 * pkin(2) - pkin(1);
t68 = pkin(3) * t118 + t147;
t105 = -t65 * qJ(5) + t68;
t188 = pkin(4) + pkin(9);
t23 = t188 * t64 + t105;
t66 = -t199 * t95 - t200 * t96;
t112 = -t79 * pkin(8) + t66;
t61 = t184 * t112;
t67 = t199 * t96 - t200 * t95;
t62 = -pkin(8) * t118 + t67;
t31 = t62 * t98 - t61;
t24 = pkin(5) * t65 + t31;
t12 = -t23 * t97 + t24 * t99;
t13 = t23 * t99 + t24 * t97;
t111 = pkin(3) * t114;
t164 = qJD(4) * t98;
t46 = t114 * t98 + t164 * t79 - (-qJD(3) - qJD(4)) * t110;
t19 = t47 * pkin(4) + t46 * qJ(5) - t65 * qJD(5) + t111;
t10 = t47 * pkin(9) + t19;
t59 = -t202 * t96 + t203 * t95;
t101 = -pkin(8) * t114 + t59;
t60 = t202 * t95 + t203 * t96;
t102 = pkin(8) * t113 + t60;
t108 = t98 * t112;
t18 = qJD(4) * t108 + t101 * t98 - t184 * t102 + t141 * t62;
t9 = -t46 * pkin(5) + t18;
t1 = qJD(6) * t12 + t10 * t99 + t9 * t97;
t186 = t1 * t97;
t163 = qJD(6) * t13;
t2 = -t10 * t97 + t9 * t99 - t163;
t20 = mrSges(7,2) * t46 - mrSges(7,3) * t119;
t21 = -mrSges(7,1) * t46 - mrSges(7,3) * t120;
t104 = m(7) * (t186 + t2 * t99 + (-t12 * t97 + t13 * t99) * qJD(6)) + t99 * t21 + t97 * t20;
t171 = t97 * mrSges(7,3);
t48 = mrSges(7,1) * t65 - t171 * t64;
t169 = t99 * mrSges(7,3);
t49 = -mrSges(7,2) * t65 + t169 * t64;
t198 = t49 * t161 - t48 * t162 + t104;
t194 = 2 * m(6);
t193 = 0.2e1 * m(7);
t80 = -mrSges(7,1) * t161 + mrSges(7,2) * t162;
t191 = -0.2e1 * t80;
t187 = pkin(3) * t98;
t179 = t46 * mrSges(6,1);
t176 = t64 * t97;
t175 = t64 * t99;
t137 = pkin(3) * t141;
t86 = t137 + qJD(5);
t88 = qJ(5) + t187;
t172 = t86 * t88;
t154 = pkin(3) * t164;
t153 = t184 * pkin(3);
t148 = t175 / 0.2e1;
t143 = t47 * mrSges(5,1) - t46 * mrSges(5,2);
t142 = -t47 * mrSges(6,2) + t46 * mrSges(6,3);
t138 = t65 * t154;
t89 = -t153 - pkin(4);
t132 = mrSges(7,1) * t99 - mrSges(7,2) * t97;
t131 = Ifges(7,1) * t97 + t180;
t130 = Ifges(7,2) * t99 + t181;
t129 = Ifges(7,5) * t97 + Ifges(7,6) * t99;
t128 = t12 * t99 + t13 * t97;
t17 = -qJD(4) * t61 - t184 * t101 - t98 * t102 + t164 * t62;
t32 = t184 * t62 + t108;
t127 = -t32 * t17 + t31 * t18;
t126 = t166 * t154;
t125 = qJ(5) * t86 + qJD(5) * t88;
t117 = Ifges(7,5) * t120 - Ifges(7,6) * t119 - Ifges(7,3) * t46;
t81 = t130 * qJD(6);
t82 = t131 * qJD(6);
t109 = -qJD(6) * t205 + t97 * t81 - t99 * t82;
t107 = mrSges(4,1) * t114 - mrSges(4,2) * t113;
t25 = -t64 * pkin(5) + t32;
t28 = Ifges(7,6) * t65 + t130 * t64;
t29 = Ifges(7,5) * t65 + t131 * t64;
t5 = Ifges(7,4) * t120 - Ifges(7,2) * t119 - t46 * Ifges(7,6);
t6 = Ifges(7,1) * t120 - Ifges(7,4) * t119 - t46 * Ifges(7,5);
t8 = -pkin(5) * t47 - t17;
t103 = t12 * mrSges(7,3) * t162 - t28 * t161 / 0.2e1 - t97 * t5 / 0.2e1 + t99 * t6 / 0.2e1 - t82 * t176 / 0.2e1 - t81 * t148 + t8 * t83 - t25 * t80 + (-Ifges(7,5) * t99 / 0.2e1 + Ifges(7,6) * t97 / 0.2e1 + Ifges(6,4) - Ifges(5,5)) * t46 + t206 * t18 - (t64 * t84 + t29) * t162 / 0.2e1 + (t85 * t148 - t65 * t129 / 0.2e1) * qJD(6) + (Ifges(6,5) - Ifges(5,6) + t205 / 0.2e1) * t47;
t87 = -pkin(9) + t89;
t38 = t132 * t64;
t30 = t64 * pkin(4) + t105;
t14 = mrSges(7,1) * t119 + mrSges(7,2) * t120;
t3 = [0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * qJD(2) * (t95 ^ 2 + t96 ^ 2) - t119 * t28 + t120 * t29 + 0.2e1 * t147 * t107 + 0.2e1 * ((-Ifges(5,4) - Ifges(6,6)) * t65 + (Ifges(6,3) + Ifges(5,2)) * t64) * t47 + 0.2e1 * m(5) * (t111 * t68 + t127) + 0.2e1 * t30 * t142 + 0.2e1 * t68 * t143 + (-t113 * t66 + t118 * t59 + t60 * t79) * t209 - (t118 * t208 + t207 * t79) * t113 + (-t118 * t207 + t208 * t79 + t209 * t67) * t114 + ((-(2 * Ifges(5,1)) - (2 * Ifges(6,2)) - Ifges(7,3)) * t65 + (0.2e1 * Ifges(5,4) + 0.2e1 * Ifges(6,6) - t129) * t64) * t46 + 0.2e1 * t19 * (-mrSges(6,2) * t64 - mrSges(6,3) * t65) + 0.2e1 * t1 * t49 + 0.2e1 * t2 * t48 - 0.2e1 * t8 * t38 + 0.2e1 * t25 * t14 + 0.2e1 * t13 * t20 + 0.2e1 * t12 * t21 + 0.2e1 * (mrSges(5,1) * t64 + mrSges(5,2) * t65) * t111 + 0.2e1 * m(4) * (t59 * t67 + t60 * t66) + t65 * t117 + (t1 * t13 + t12 * t2 + t25 * t8) * t193 + (t19 * t30 + t127) * t194 + 0.2e1 * (mrSges(5,3) + mrSges(6,1)) * (t17 * t64 + t18 * t65 - t31 * t46 - t32 * t47) + t5 * t175 + t6 * t176; m(7) * (-qJD(6) * t128 + t1 * t99 - t2 * t97) - t49 * t162 + t99 * t20 - t48 * t161 - t97 * t21 + m(6) * t19 + m(5) * t111 + t107 + t142 + t143; 0; t103 + m(7) * (t25 * t86 + t8 * t88) + m(6) * (-t17 * t88 + t18 * t89 + t32 * t86) + m(5) * (-t184 * t18 - t17 * t98 + (t184 * t32 + t31 * t98) * qJD(4)) * pkin(3) - t13 * mrSges(7,3) * t161 - t86 * t38 + t88 * t14 - t59 * mrSges(4,2) + t60 * mrSges(4,1) - t17 * mrSges(6,3) + t17 * mrSges(5,2) - t89 * t179 - t2 * t169 - t1 * t171 - Ifges(4,5) * t113 - Ifges(4,6) * t114 + (m(6) * t31 + m(7) * t128 + t99 * t48 + t97 * t49) * t154 + (-t47 * t88 - t64 * t86 + t138) * mrSges(6,1) + t198 * t87 + (-t137 * t64 + t153 * t46 - t187 * t47 + t138) * mrSges(5,3); 0; t88 * t191 + 0.2e1 * t201 * t86 + 0.2e1 * t204 + (t126 * t87 + t172) * t193 + (t154 * t89 + t172) * t194 + t109; t103 - t198 * t188 - qJD(5) * t38 + qJ(5) * t14 + (-t186 + (-t2 - t163) * t99) * mrSges(7,3) + (mrSges(5,2) - mrSges(6,3)) * t17 + m(7) * (qJ(5) * t8 + qJD(5) * t25) + m(6) * (-pkin(4) * t18 - qJ(5) * t17 + qJD(5) * t32) + (pkin(4) * t46 - qJ(5) * t47 - qJD(5) * t64) * mrSges(6,1); 0; (-qJ(5) - t88) * t80 + t204 + m(7) * (-t126 * t188 + t125) + m(6) * (-pkin(4) * t154 + t125) + t109 + t201 * (qJD(5) + t86); qJ(5) * t191 + 0.2e1 * ((m(6) + m(7)) * qJ(5) + t201) * qJD(5) + t109; -t179 + (-t97 * t48 + t99 * t49) * qJD(6) + m(6) * t18 + t104; 0; (m(7) * t166 + m(6)) * t154; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t117; t80; t132 * t154 + (-t83 * t87 - t129) * qJD(6); (t188 * t83 - t129) * qJD(6); -t83 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
