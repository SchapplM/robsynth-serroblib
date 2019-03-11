% Calculate joint inertia matrix for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:05:54
% EndTime: 2019-03-10 01:05:58
% DurationCPUTime: 1.68s
% Computational Cost: add. (2728->314), mult. (5170->421), div. (0->0), fcn. (5607->8), ass. (0->120)
t222 = Ifges(6,1) + Ifges(7,1);
t221 = Ifges(6,4) + Ifges(7,4);
t220 = Ifges(6,5) + Ifges(7,5);
t219 = Ifges(6,2) + Ifges(7,2);
t218 = Ifges(6,6) + Ifges(7,6);
t144 = sin(qJ(4));
t148 = cos(qJ(4));
t143 = sin(qJ(5));
t147 = cos(qJ(5));
t113 = -t143 * t144 + t147 * t148;
t115 = t143 * t148 + t144 * t147;
t165 = t113 * t218 + t115 * t220;
t217 = Ifges(5,5) * t144 + Ifges(5,6) * t148 + t165;
t187 = Ifges(6,3) + Ifges(7,3);
t145 = sin(qJ(3));
t146 = sin(qJ(2));
t149 = cos(qJ(3));
t150 = cos(qJ(2));
t114 = t145 * t146 - t149 * t150;
t116 = t145 * t150 + t146 * t149;
t60 = t115 * t116;
t61 = t113 * t116;
t216 = t114 * t218 - t219 * t60 + t221 * t61;
t215 = t114 * t220 - t221 * t60 + t222 * t61;
t214 = t113 * t219 + t115 * t221;
t213 = t113 * t221 + t115 * t222;
t134 = -pkin(2) * t150 - pkin(1);
t69 = pkin(3) * t114 - pkin(9) * t116 + t134;
t198 = -pkin(8) - pkin(7);
t127 = t198 * t150;
t168 = t198 * t146;
t87 = -t149 * t127 + t145 * t168;
t32 = -t144 * t87 + t148 * t69;
t33 = t144 * t69 + t148 * t87;
t163 = -t144 * t32 + t148 * t33;
t178 = t116 * t148;
t212 = Ifges(5,5) * t178 + Ifges(5,3) * t114;
t211 = (t147 * mrSges(6,1) + (-mrSges(6,2) - mrSges(7,2)) * t143) * pkin(4);
t180 = qJ(6) * t115;
t130 = pkin(2) * t145 + pkin(9);
t109 = (-pkin(10) - t130) * t144;
t138 = t148 * pkin(10);
t110 = t130 * t148 + t138;
t66 = t109 * t147 - t110 * t143;
t41 = t66 - t180;
t67 = t109 * t143 + t110 * t147;
t98 = t113 * qJ(6);
t42 = t98 + t67;
t210 = mrSges(6,1) * t66 + mrSges(7,1) * t41 - t67 * mrSges(6,2) - t42 * mrSges(7,2);
t125 = (-pkin(10) - pkin(9)) * t144;
t126 = pkin(9) * t148 + t138;
t83 = t125 * t147 - t126 * t143;
t53 = t83 - t180;
t86 = t125 * t143 + t126 * t147;
t54 = t98 + t86;
t209 = mrSges(6,1) * t83 + mrSges(7,1) * t53 - t86 * mrSges(6,2) - t54 * mrSges(7,2);
t84 = -t127 * t145 - t149 * t168;
t208 = t84 ^ 2;
t207 = 0.2e1 * mrSges(7,1);
t72 = -mrSges(7,1) * t113 + mrSges(7,2) * t115;
t206 = 0.2e1 * t72;
t73 = -mrSges(6,1) * t113 + mrSges(6,2) * t115;
t205 = 0.2e1 * t73;
t204 = 0.2e1 * t84;
t203 = 0.2e1 * t134;
t202 = m(6) * pkin(4);
t201 = m(7) * pkin(5);
t194 = pkin(2) * t149;
t193 = pkin(4) * t143;
t192 = pkin(4) * t147;
t16 = pkin(4) * t114 - pkin(10) * t178 + t32;
t179 = t116 * t144;
t22 = -pkin(10) * t179 + t33;
t8 = t143 * t16 + t147 * t22;
t186 = Ifges(5,4) * t144;
t185 = Ifges(5,4) * t148;
t184 = t115 * mrSges(7,3);
t174 = t144 ^ 2 + t148 ^ 2;
t173 = t146 ^ 2 + t150 ^ 2;
t172 = 2 * mrSges(6,3);
t171 = 0.2e1 * mrSges(7,3);
t7 = -t143 * t22 + t147 * t16;
t3 = pkin(5) * t114 - qJ(6) * t61 + t7;
t36 = mrSges(7,1) * t114 - mrSges(7,3) * t61;
t169 = m(7) * t3 + t36;
t133 = -pkin(4) * t148 - pkin(3);
t24 = mrSges(7,1) * t60 + mrSges(7,2) * t61;
t166 = t174 * t130;
t51 = pkin(4) * t179 + t84;
t164 = mrSges(5,1) * t144 + mrSges(5,2) * t148;
t70 = -mrSges(5,2) * t114 - mrSges(5,3) * t179;
t71 = mrSges(5,1) * t114 - mrSges(5,3) * t178;
t162 = -t144 * t71 + t148 * t70;
t161 = 0.2e1 * mrSges(5,3) * t174;
t90 = -pkin(5) * t113 + t133;
t160 = t114 * t187 - t218 * t60 + t220 * t61;
t159 = (mrSges(4,1) * t149 - mrSges(4,2) * t145) * pkin(2);
t123 = Ifges(5,2) * t148 + t186;
t124 = Ifges(5,1) * t144 + t185;
t158 = t113 * t214 + t115 * t213 + t148 * t123 + t144 * t124 + Ifges(4,3);
t4 = -qJ(6) * t60 + t8;
t157 = mrSges(6,1) * t7 + mrSges(7,1) * t3 - t8 * mrSges(6,2) - t4 * mrSges(7,2) + t160;
t131 = pkin(5) + t192;
t156 = (-mrSges(6,3) * t192 - mrSges(7,3) * t131) * t115 + (mrSges(7,3) + mrSges(6,3)) * t113 * t193 + t217;
t122 = -mrSges(5,1) * t148 + mrSges(5,2) * t144;
t23 = pkin(5) * t60 + t51;
t43 = Ifges(5,6) * t114 + (-Ifges(5,2) * t144 + t185) * t116;
t44 = Ifges(5,5) * t114 + (Ifges(5,1) * t148 - t186) * t116;
t155 = t124 * t178 / 0.2e1 - t123 * t179 / 0.2e1 + t148 * t43 / 0.2e1 + t144 * t44 / 0.2e1 + Ifges(4,5) * t116 - t87 * mrSges(4,2) + t23 * t72 + t51 * t73 + (t122 - mrSges(4,1)) * t84 - t214 * t60 / 0.2e1 + t213 * t61 / 0.2e1 + t163 * mrSges(5,3) + (-t7 * mrSges(6,3) - t3 * mrSges(7,3) + t215 / 0.2e1) * t115 + (t216 / 0.2e1 + mrSges(6,3) * t8 + mrSges(7,3) * t4) * t113 + (-Ifges(4,6) + t217 / 0.2e1) * t114;
t153 = pkin(4) ^ 2;
t137 = t143 ^ 2 * t153;
t132 = -pkin(3) - t194;
t121 = t133 - t194;
t88 = t90 - t194;
t68 = t164 * t116;
t37 = mrSges(6,1) * t114 - mrSges(6,3) * t61;
t35 = -mrSges(6,2) * t114 - mrSges(6,3) * t60;
t34 = -mrSges(7,2) * t114 - mrSges(7,3) * t60;
t25 = mrSges(6,1) * t60 + mrSges(6,2) * t61;
t1 = [t146 * (Ifges(3,1) * t146 + Ifges(3,4) * t150) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t150 + mrSges(3,2) * t146) + t150 * (Ifges(3,4) * t146 + Ifges(3,2) * t150) + 0.2e1 * t33 * t70 + 0.2e1 * t32 * t71 + t68 * t204 + 0.2e1 * t51 * t25 + 0.2e1 * t4 * t34 + 0.2e1 * t8 * t35 + 0.2e1 * t3 * t36 + 0.2e1 * t7 * t37 + 0.2e1 * t23 * t24 + Ifges(2,3) + t215 * t61 - t216 * t60 + 0.2e1 * t173 * pkin(7) * mrSges(3,3) + (mrSges(4,2) * t203 + mrSges(4,3) * t204 + Ifges(4,1) * t116 - t144 * t43 + t148 * t44) * t116 + (mrSges(4,1) * t203 - 0.2e1 * t87 * mrSges(4,3) + Ifges(4,2) * t114 + (-Ifges(5,6) * t144 - (2 * Ifges(4,4))) * t116 + t160 + t212) * t114 + m(3) * (pkin(7) ^ 2 * t173 + pkin(1) ^ 2) + m(4) * (t134 ^ 2 + t87 ^ 2 + t208) + m(5) * (t32 ^ 2 + t33 ^ 2 + t208) + m(6) * (t51 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(7) * (t23 ^ 2 + t3 ^ 2 + t4 ^ 2); (-mrSges(3,1) * t146 - mrSges(3,2) * t150) * pkin(7) + Ifges(3,5) * t146 + Ifges(3,6) * t150 + t132 * t68 + t121 * t25 + t88 * t24 + t66 * t37 + t67 * t35 + t41 * t36 + t42 * t34 + t155 + m(5) * (t130 * t163 + t132 * t84) + t162 * t130 + m(6) * (t121 * t51 + t66 * t7 + t67 * t8) + m(7) * (t23 * t88 + t3 * t41 + t4 * t42) + (m(4) * (t145 * t87 - t149 * t84) + (-t114 * t145 - t116 * t149) * mrSges(4,3)) * pkin(2); m(4) * (t145 ^ 2 + t149 ^ 2) * pkin(2) ^ 2 + 0.2e1 * t132 * t122 + t121 * t205 + t88 * t206 + 0.2e1 * t159 + (t113 * t42 - t115 * t41) * t171 + (t113 * t67 - t115 * t66) * t172 + m(5) * (t130 ^ 2 * t174 + t132 ^ 2) + t158 + m(7) * (t41 ^ 2 + t42 ^ 2 + t88 ^ 2) + m(6) * (t121 ^ 2 + t66 ^ 2 + t67 ^ 2) + t130 * t161 + Ifges(3,3); t162 * pkin(9) + m(6) * (t133 * t51 + t7 * t83 + t8 * t86) + m(7) * (t23 * t90 + t3 * t53 + t4 * t54) + t133 * t25 + t90 * t24 + t83 * t37 + t86 * t35 - pkin(3) * t68 + t53 * t36 + t54 * t34 + t155 + m(5) * (-pkin(3) * t84 + pkin(9) * t163); m(5) * (-pkin(3) * t132 + pkin(9) * t166) + t158 + t159 + (t132 - pkin(3)) * t122 + ((-t41 - t53) * t115 + (t42 + t54) * t113) * mrSges(7,3) + ((-t66 - t83) * t115 + (t67 + t86) * t113) * mrSges(6,3) + (t88 + t90) * t72 + (pkin(9) * t174 + t166) * mrSges(5,3) + (t133 + t121) * t73 + m(6) * (t121 * t133 + t66 * t83 + t67 * t86) + m(7) * (t41 * t53 + t42 * t54 + t88 * t90); -0.2e1 * pkin(3) * t122 + t133 * t205 + t90 * t206 + (t113 * t54 - t53 * t115) * t171 + (t86 * t113 - t83 * t115) * t172 + pkin(9) * t161 + m(7) * (t53 ^ 2 + t54 ^ 2 + t90 ^ 2) + m(6) * (t133 ^ 2 + t83 ^ 2 + t86 ^ 2) + m(5) * (pkin(9) ^ 2 * t174 + pkin(3) ^ 2) + t158; -Ifges(5,6) * t179 + t32 * mrSges(5,1) - t33 * mrSges(5,2) + t157 + t169 * t131 + ((m(6) * t7 + t37) * t147 + (m(6) * t8 + m(7) * t4 + t34 + t35) * t143) * pkin(4) + t212; m(7) * (t131 * t41 + t193 * t42) - t164 * t130 + t156 + (t143 * t67 + t147 * t66) * t202 + t210; m(7) * (t131 * t53 + t193 * t54) - t164 * pkin(9) + t156 + (t143 * t86 + t147 * t83) * t202 + t209; t131 * t207 + Ifges(5,3) + m(6) * (t147 ^ 2 * t153 + t137) + m(7) * (t131 ^ 2 + t137) + 0.2e1 * t211 + t187; pkin(5) * t169 + t157; (m(7) * t41 - t184) * pkin(5) + t165 + t210; (m(7) * t53 - t184) * pkin(5) + t165 + t209; t131 * t201 + (pkin(5) + t131) * mrSges(7,1) + t211 + t187; (t207 + t201) * pkin(5) + t187; m(7) * t23 + t24; m(7) * t88 + t72; m(7) * t90 + t72; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
