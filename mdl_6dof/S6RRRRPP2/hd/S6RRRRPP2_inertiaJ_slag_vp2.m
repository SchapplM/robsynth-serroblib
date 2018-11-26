% Calculate joint inertia matrix for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:05
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:04:51
% EndTime: 2018-11-23 18:04:53
% DurationCPUTime: 1.63s
% Computational Cost: add. (1519->293), mult. (2858->374), div. (0->0), fcn. (2727->6), ass. (0->106)
t195 = Ifges(7,4) + Ifges(6,5);
t194 = Ifges(6,4) + Ifges(5,5);
t193 = Ifges(7,2) + Ifges(6,3);
t118 = cos(qJ(4));
t192 = t195 * t118;
t115 = sin(qJ(4));
t191 = t195 * t115;
t190 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t189 = Ifges(6,2) + Ifges(5,3);
t116 = sin(qJ(3));
t117 = sin(qJ(2));
t119 = cos(qJ(3));
t120 = cos(qJ(2));
t73 = t116 * t117 - t119 * t120;
t74 = t116 * t120 + t117 * t119;
t188 = (t115 * t193 + t192) * t74 + (-Ifges(7,6) + Ifges(6,6)) * t73;
t187 = -t118 * t193 + t191;
t186 = Ifges(5,6) * t118 + t115 * t194;
t185 = t115 ^ 2 + t118 ^ 2;
t160 = Ifges(5,4) * t115;
t168 = t73 * Ifges(7,5);
t184 = -t168 + t194 * t73 + (t118 * t190 - t160 + t191) * t74;
t159 = Ifges(5,4) * t118;
t183 = t115 * t190 + t159 - t192;
t182 = pkin(4) * t118 + qJ(5) * t115;
t181 = (mrSges(5,3) + mrSges(6,2)) * t185;
t173 = -pkin(8) - pkin(7);
t142 = t173 * t117;
t90 = t173 * t120;
t47 = -t116 * t90 - t119 * t142;
t180 = t47 ^ 2;
t179 = -2 * mrSges(7,3);
t178 = 0.2e1 * t47;
t80 = -mrSges(6,1) * t118 - mrSges(6,3) * t115;
t177 = 0.2e1 * t80;
t81 = mrSges(7,1) * t118 + mrSges(7,2) * t115;
t176 = 0.2e1 * t81;
t97 = -pkin(2) * t120 - pkin(1);
t175 = 0.2e1 * t97;
t121 = -pkin(4) - pkin(5);
t31 = pkin(3) * t73 - pkin(9) * t74 + t97;
t49 = t116 * t142 - t119 * t90;
t12 = t115 * t31 + t118 * t49;
t6 = qJ(5) * t73 + t12;
t169 = t118 * t6;
t167 = -Ifges(5,6) - Ifges(7,6);
t166 = pkin(9) - qJ(6);
t153 = t115 * t74;
t33 = -mrSges(5,2) * t73 - mrSges(5,3) * t153;
t37 = -mrSges(6,2) * t153 + mrSges(6,3) * t73;
t165 = t33 + t37;
t151 = t118 * t74;
t35 = mrSges(5,1) * t73 - mrSges(5,3) * t151;
t36 = -t73 * mrSges(6,1) + mrSges(6,2) * t151;
t164 = -t35 + t36;
t95 = pkin(2) * t116 + pkin(9);
t163 = t185 * pkin(9) * t95;
t161 = t185 * t95 ^ 2;
t38 = t115 * t49;
t11 = t118 * t31 - t38;
t154 = t11 * t115;
t101 = t115 * mrSges(6,2);
t152 = t118 * t12;
t150 = -qJ(6) + t95;
t149 = qJ(5) * t118;
t147 = t185 * pkin(9) ^ 2;
t145 = t117 ^ 2 + t120 ^ 2;
t144 = t115 * t179;
t143 = pkin(3) + t182;
t96 = -pkin(2) * t119 - pkin(3);
t29 = -mrSges(7,1) * t153 + mrSges(7,2) * t151;
t34 = -mrSges(7,1) * t73 - mrSges(7,3) * t151;
t86 = Ifges(5,2) * t118 + t160;
t138 = t115 * t183 + t118 * t86 + Ifges(4,3);
t137 = Ifges(6,6) * t153 + t151 * t194 + t189 * t73;
t7 = -pkin(4) * t73 - t11;
t136 = t115 * t7 + t169;
t135 = t115 * mrSges(5,1) + t118 * mrSges(5,2);
t134 = t115 * mrSges(6,1) - t118 * mrSges(6,3);
t133 = pkin(4) * t115 - t149;
t132 = t152 - t154;
t59 = t96 - t182;
t131 = (mrSges(4,1) * t119 - mrSges(4,2) * t116) * pkin(2);
t129 = 0.2e1 * t181;
t128 = -m(6) * t133 - t134 - t135;
t103 = Ifges(7,6) * t118;
t127 = t103 + mrSges(6,2) * t149 + (-mrSges(7,3) * qJ(5) - Ifges(6,6)) * t118 + (-mrSges(6,2) * pkin(4) - mrSges(7,3) * t121 - Ifges(7,5)) * t115 + t186;
t13 = t133 * t74 + t47;
t22 = t73 * Ifges(5,6) + (-Ifges(5,2) * t115 + t159) * t74;
t8 = (t115 * t121 + t149) * t74 - t47;
t82 = -mrSges(5,1) * t118 + mrSges(5,2) * t115;
t126 = t118 * t22 / 0.2e1 - t73 * (Ifges(7,5) * t115 - t103) / 0.2e1 + t13 * t80 + t8 * t81 - Ifges(4,6) * t73 + Ifges(4,5) * t74 - t49 * mrSges(4,2) + mrSges(5,3) * t152 + t7 * t101 + mrSges(6,2) * t169 + (t82 - mrSges(4,1)) * t47 + (-Ifges(6,6) * t118 + t186) * t73 / 0.2e1 - t188 * t118 / 0.2e1 + t184 * t115 / 0.2e1 + t183 * t151 / 0.2e1 + (-t86 / 0.2e1 + t187 / 0.2e1) * t153;
t122 = qJ(5) ^ 2;
t109 = t118 * pkin(5);
t83 = t166 * t118;
t79 = t166 * t115;
t63 = t150 * t118;
t62 = t150 * t115;
t60 = t109 + t143;
t53 = t109 - t59;
t32 = mrSges(7,2) * t73 + mrSges(7,3) * t153;
t30 = t135 * t74;
t28 = t134 * t74;
t2 = qJ(6) * t153 + t6;
t1 = t38 + t121 * t73 + (-qJ(6) * t74 - t31) * t118;
t3 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t120 + mrSges(3,2) * t117) + t120 * (Ifges(3,4) * t117 + Ifges(3,2) * t120) + t117 * (Ifges(3,1) * t117 + Ifges(3,4) * t120) + t30 * t178 + 0.2e1 * t13 * t28 + 0.2e1 * t8 * t29 + 0.2e1 * t2 * t32 + 0.2e1 * t12 * t33 + 0.2e1 * t1 * t34 + 0.2e1 * t11 * t35 + 0.2e1 * t7 * t36 + 0.2e1 * t6 * t37 + Ifges(2,3) + 0.2e1 * t145 * pkin(7) * mrSges(3,3) + m(3) * (pkin(7) ^ 2 * t145 + pkin(1) ^ 2) + m(4) * (t49 ^ 2 + t97 ^ 2 + t180) + m(5) * (t11 ^ 2 + t12 ^ 2 + t180) + m(7) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) + m(6) * (t13 ^ 2 + t6 ^ 2 + t7 ^ 2) + (mrSges(4,1) * t175 - 0.2e1 * t49 * mrSges(4,3) + (Ifges(7,3) + Ifges(4,2)) * t73 + t137) * t73 + (mrSges(4,2) * t175 + mrSges(4,3) * t178 + Ifges(4,1) * t74 - 0.2e1 * Ifges(4,4) * t73 + (-t168 + t184) * t118 + (t167 * t73 + t188 - t22) * t115) * t74; m(6) * (t13 * t59 + t136 * t95) + m(5) * (t132 * t95 + t47 * t96) + (-t2 * mrSges(7,3) + t165 * t95) * t118 + (-t11 * mrSges(5,3) - t1 * mrSges(7,3) + t164 * t95) * t115 + Ifges(3,6) * t120 + Ifges(3,5) * t117 + t96 * t30 + m(7) * (t1 * t62 + t2 * t63 + t53 * t8) + t126 + t59 * t28 + t62 * t34 + t63 * t32 + t53 * t29 + (-mrSges(3,1) * t117 - mrSges(3,2) * t120) * pkin(7) + (m(4) * (t116 * t49 - t119 * t47) + (-t116 * t73 - t119 * t74) * mrSges(4,3)) * pkin(2); t62 * t144 + t53 * t176 + t59 * t177 + 0.2e1 * t96 * t82 + Ifges(3,3) + 0.2e1 * t131 + (t179 * t63 - t187) * t118 + t129 * t95 + m(7) * (t53 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(6) * (t59 ^ 2 + t161) + m(5) * (t96 ^ 2 + t161) + m(4) * (t116 ^ 2 + t119 ^ 2) * pkin(2) ^ 2 + t138; (-t1 * t115 - t118 * t2) * mrSges(7,3) + m(6) * (pkin(9) * t136 - t13 * t143) + m(5) * (-pkin(3) * t47 + pkin(9) * t132) - mrSges(5,3) * t154 + t126 + t79 * t34 + t83 * t32 - t143 * t28 + t60 * t29 - pkin(3) * t30 + (t115 * t164 + t118 * t165) * pkin(9) + m(7) * (t1 * t79 + t2 * t83 + t60 * t8); (t96 - pkin(3)) * t82 + (t60 + t53) * t81 + (-t143 + t59) * t80 - t187 * t118 + t131 + m(5) * (-pkin(3) * t96 + t163) + m(6) * (-t143 * t59 + t163) + m(7) * (t53 * t60 + t62 * t79 + t63 * t83) + ((-t63 - t83) * t118 + (-t62 - t79) * t115) * mrSges(7,3) + t138 + (pkin(9) + t95) * t181; t79 * t144 - 0.2e1 * pkin(3) * t82 + t60 * t176 - t143 * t177 + (t179 * t83 - t187) * t118 + m(6) * (t143 ^ 2 + t147) + m(7) * (t60 ^ 2 + t79 ^ 2 + t83 ^ 2) + m(5) * (pkin(3) ^ 2 + t147) + t129 * pkin(9) + t138; t11 * mrSges(5,1) - t7 * mrSges(6,1) - t1 * mrSges(7,1) - t12 * mrSges(5,2) + t2 * mrSges(7,2) + t6 * mrSges(6,3) + Ifges(7,3) * t73 - pkin(4) * t36 + t121 * t34 + (t32 + t37) * qJ(5) + m(7) * (qJ(5) * t2 + t1 * t121) + m(6) * (-pkin(4) * t7 + qJ(5) * t6) + (-Ifges(7,5) * t118 + t115 * t167) * t74 + t137; m(7) * (qJ(5) * t63 + t121 * t62) - t62 * mrSges(7,1) + t63 * mrSges(7,2) + t128 * t95 + t127; m(7) * (qJ(5) * t83 + t121 * t79) - t79 * mrSges(7,1) + t83 * mrSges(7,2) + t128 * pkin(9) + t127; 0.2e1 * pkin(4) * mrSges(6,1) - 0.2e1 * t121 * mrSges(7,1) + Ifges(7,3) + 0.2e1 * (mrSges(7,2) + mrSges(6,3)) * qJ(5) + m(7) * (t121 ^ 2 + t122) + m(6) * (pkin(4) ^ 2 + t122) + t189; m(6) * t7 + m(7) * t1 + t34 + t36; m(7) * t62 + t101 + (m(6) * t95 - mrSges(7,3)) * t115; m(7) * t79 + t101 + (m(6) * pkin(9) - mrSges(7,3)) * t115; -m(6) * pkin(4) + m(7) * t121 - mrSges(6,1) - mrSges(7,1); m(6) + m(7); m(7) * t8 + t29; m(7) * t53 + t81; m(7) * t60 + t81; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
