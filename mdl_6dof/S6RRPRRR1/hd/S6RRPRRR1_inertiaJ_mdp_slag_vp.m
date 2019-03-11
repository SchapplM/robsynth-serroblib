% Calculate joint inertia matrix for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR1_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:14:31
% EndTime: 2019-03-09 13:14:33
% DurationCPUTime: 0.52s
% Computational Cost: add. (1463->145), mult. (2714->196), div. (0->0), fcn. (3282->10), ass. (0->83)
t150 = sin(qJ(6));
t154 = cos(qJ(6));
t183 = MDP(29) * t150 + MDP(30) * t154;
t175 = t154 * MDP(32);
t164 = -MDP(33) * t150 + t175;
t160 = 0.2e1 * t164;
t148 = sin(pkin(11));
t149 = cos(pkin(11));
t153 = sin(qJ(2));
t157 = cos(qJ(2));
t128 = -t148 * t153 + t149 * t157;
t129 = t148 * t157 + t149 * t153;
t152 = sin(qJ(4));
t156 = cos(qJ(4));
t119 = -t128 * t156 + t129 * t152;
t141 = -pkin(2) * t157 - pkin(1);
t124 = -pkin(3) * t128 + t141;
t107 = pkin(4) * t119 + t124;
t193 = 0.2e1 * t107;
t192 = 0.2e1 * t124;
t191 = 0.2e1 * t157;
t190 = pkin(2) * t148;
t189 = pkin(5) * t150;
t155 = cos(qJ(5));
t188 = t155 * pkin(4);
t187 = -qJ(3) - pkin(7);
t151 = sin(qJ(5));
t120 = t128 * t152 + t129 * t156;
t133 = t187 * t153;
t134 = t187 * t157;
t121 = t133 * t149 + t134 * t148;
t114 = -pkin(8) * t129 + t121;
t122 = t133 * t148 - t134 * t149;
t115 = pkin(8) * t128 + t122;
t172 = t114 * t156 - t115 * t152;
t95 = -pkin(9) * t120 + t172;
t166 = -t114 * t152 - t115 * t156;
t96 = -pkin(9) * t119 - t166;
t91 = t151 * t96 - t155 * t95;
t88 = t91 * t150;
t186 = t91 * t154;
t185 = t150 * t154;
t138 = pkin(2) * t149 + pkin(3);
t127 = t138 * t152 + t156 * t190;
t184 = t155 * t127;
t105 = t119 * t155 + t120 * t151;
t182 = t105 * MDP(31);
t106 = -t119 * t151 + t120 * t155;
t181 = t106 * MDP(26);
t126 = t138 * t156 - t152 * t190;
t125 = pkin(4) + t126;
t171 = -t125 * t155 + t127 * t151;
t180 = t171 * MDP(25);
t113 = -t125 * t151 - t184;
t179 = t113 * MDP(26);
t178 = t119 * MDP(18);
t177 = t126 * MDP(18);
t176 = t127 * MDP(19);
t174 = MDP(28) * t185;
t146 = t150 ^ 2;
t173 = MDP(27) * t146 + MDP(24) + 0.2e1 * t174;
t170 = MDP(17) + t173;
t169 = -pkin(5) * t106 - pkin(10) * t105;
t110 = -pkin(5) + t171;
t111 = pkin(10) - t113;
t168 = -t105 * t111 + t106 * t110;
t139 = pkin(4) * t151 + pkin(10);
t140 = -pkin(5) - t188;
t167 = -t105 * t139 + t106 * t140;
t165 = MDP(29) * t154 - MDP(30) * t150;
t163 = -MDP(32) * t150 - MDP(33) * t154;
t162 = (MDP(25) * t155 - MDP(26) * t151) * pkin(4);
t147 = t154 ^ 2;
t92 = t151 * t95 + t155 * t96;
t161 = -t91 * MDP(25) - t92 * MDP(26) + (MDP(22) + (-t146 + t147) * MDP(28) + MDP(27) * t185) * t106 + (-MDP(23) + t183) * t105;
t159 = MDP(15) * t120 - MDP(16) * t119 + MDP(18) * t172 + MDP(19) * t166 + t161;
t145 = pkin(5) * t154;
t136 = t140 * t150;
t108 = t110 * t150;
t93 = pkin(5) * t105 - pkin(10) * t106 + t107;
t87 = t150 * t93 + t154 * t92;
t86 = -t150 * t92 + t154 * t93;
t1 = [MDP(1) + pkin(1) * MDP(9) * t191 + (t121 ^ 2 + t122 ^ 2 + t141 ^ 2) * MDP(12) + t178 * t192 + t181 * t193 + (-0.2e1 * MDP(10) * pkin(1) + MDP(4) * t153 + MDP(5) * t191) * t153 + (MDP(13) * t120 - 0.2e1 * MDP(14) * t119 + MDP(19) * t192) * t120 + (MDP(27) * t147 + MDP(20) - 0.2e1 * t174) * t106 ^ 2 + (MDP(25) * t193 + t182 + 0.2e1 * (-MDP(21) + t165) * t106) * t105 + 0.2e1 * (-t121 * t129 + t122 * t128) * MDP(11) + 0.2e1 * (t105 * t86 + t106 * t88) * MDP(32) + 0.2e1 * (-t105 * t87 + t106 * t186) * MDP(33); t159 + (-MDP(10) * t157 - MDP(9) * t153) * pkin(7) + ((t128 * t148 - t129 * t149) * MDP(11) + (t121 * t149 + t122 * t148) * MDP(12)) * pkin(2) + (t150 * t168 - t186) * MDP(32) + (t154 * t168 + t88) * MDP(33) + t153 * MDP(6) + t157 * MDP(7); MDP(8) + (t148 ^ 2 + t149 ^ 2) * MDP(12) * pkin(2) ^ 2 - t110 * t160 + 0.2e1 * t177 - 0.2e1 * t176 - 0.2e1 * t180 + 0.2e1 * t179 + t170; t141 * MDP(12) + t178 + t120 * MDP(19) + t181 + (MDP(25) + t164) * t105; 0; MDP(12); t159 + (t150 * t167 - t186) * MDP(32) + (t154 * t167 + t88) * MDP(33); t177 - t176 + (-t171 + t188) * MDP(25) + (-t184 + (-pkin(4) - t125) * t151) * MDP(26) + (t136 + t108) * MDP(33) + (-t110 - t140) * t175 + t170; 0; -t140 * t160 + 0.2e1 * t162 + t170; (t150 * t169 - t186) * MDP(32) + (t154 * t169 + t88) * MDP(33) + t161; -t180 + t179 + (-t110 * t154 + t145) * MDP(32) + (t108 - t189) * MDP(33) + t173; 0; (-t140 * t154 + t145) * MDP(32) + (t136 - t189) * MDP(33) + t162 + t173; pkin(5) * t160 + t173; t86 * MDP(32) - t87 * MDP(33) + t106 * t165 + t182; t111 * t163 + t183; t164; t139 * t163 + t183; pkin(10) * t163 + t183; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
