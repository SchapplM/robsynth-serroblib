% Calculate joint inertia matrix for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRRP5_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:25:55
% EndTime: 2019-03-09 00:25:57
% DurationCPUTime: 0.83s
% Computational Cost: add. (1081->231), mult. (2619->345), div. (0->0), fcn. (2906->12), ass. (0->99)
t146 = sin(pkin(7));
t203 = 0.2e1 * t146;
t150 = sin(qJ(5));
t154 = cos(qJ(5));
t161 = t150 * MDP(24) + t154 * MDP(25);
t202 = t150 * MDP(21) + t154 * MDP(22) - t161 * pkin(11) - MDP(15);
t201 = 0.2e1 * MDP(24);
t200 = 0.2e1 * MDP(25);
t199 = 2 * MDP(26);
t152 = sin(qJ(3));
t198 = pkin(2) * t152;
t156 = cos(qJ(3));
t197 = pkin(2) * t156;
t196 = pkin(10) * t150;
t195 = pkin(10) * t154;
t155 = cos(qJ(4));
t194 = pkin(10) * t155;
t193 = -qJ(6) - pkin(11);
t192 = MDP(27) * pkin(5);
t191 = pkin(3) * MDP(17);
t190 = pkin(3) * MDP(18);
t189 = pkin(10) * MDP(18);
t151 = sin(qJ(4));
t188 = qJ(6) * t151;
t148 = cos(pkin(7));
t184 = t146 * t156;
t167 = pkin(9) * t184;
t123 = t167 + (pkin(10) + t198) * t148;
t124 = (-pkin(3) * t156 - pkin(10) * t152 - pkin(2)) * t146;
t110 = -t151 * t123 + t155 * t124;
t108 = pkin(4) * t184 - t110;
t187 = t108 * t150;
t186 = t108 * t154;
t185 = t146 * t152;
t183 = t148 * MDP(9);
t157 = cos(qJ(2));
t182 = t148 * t157;
t181 = t152 * MDP(7);
t180 = MDP(16) * t156;
t147 = sin(pkin(6));
t149 = cos(pkin(6));
t153 = sin(qJ(2));
t114 = t149 * t185 + (t152 * t182 + t153 * t156) * t147;
t125 = -t147 * t157 * t146 + t149 * t148;
t106 = t114 * t155 + t125 * t151;
t113 = -t149 * t184 + (t152 * t153 - t156 * t182) * t147;
t102 = t106 * t154 + t113 * t150;
t179 = t102 * MDP(25);
t127 = t151 * t148 + t155 * t185;
t115 = t150 * t127 + t154 * t184;
t178 = t115 * MDP(22);
t116 = t154 * t127 - t150 * t184;
t177 = t116 * MDP(19);
t176 = t116 * MDP(21);
t126 = -t155 * t148 + t151 * t185;
t175 = t126 * MDP(23);
t174 = t127 * MDP(13);
t173 = t127 * MDP(14);
t141 = -t154 * pkin(5) - pkin(4);
t172 = t141 * MDP(27);
t171 = t150 * MDP(22);
t170 = t154 * MDP(19);
t169 = t155 * MDP(23);
t168 = t154 * t194;
t166 = t150 * t154 * MDP(20);
t165 = pkin(10) * MDP(17) - MDP(14);
t164 = -MDP(26) * pkin(5) + MDP(21);
t138 = pkin(9) * t185;
t122 = t138 + (-pkin(3) - t197) * t148;
t107 = t126 * pkin(4) - t127 * pkin(11) + t122;
t111 = t155 * t123 + t151 * t124;
t109 = -pkin(11) * t184 + t111;
t99 = t154 * t107 - t150 * t109;
t100 = t150 * t107 + t154 * t109;
t133 = -t155 * pkin(4) - t151 * pkin(11) - pkin(3);
t130 = t154 * t133;
t119 = -t150 * t194 + t130;
t120 = t150 * t133 + t168;
t163 = t119 * MDP(24) - t120 * MDP(25);
t162 = t154 * MDP(24) - t150 * MDP(25);
t160 = t154 * MDP(21) - MDP(13) - t171;
t158 = t99 * MDP(24) - t100 * MDP(25) + t175 + t176 - t178;
t145 = t154 ^ 2;
t144 = t151 ^ 2;
t143 = t150 ^ 2;
t142 = t146 ^ 2;
t135 = t193 * t154;
t134 = t193 * t150;
t132 = (pkin(5) * t150 + pkin(10)) * t151;
t129 = t148 * t198 + t167;
t128 = t148 * t197 - t138;
t117 = t168 + (t133 - t188) * t150;
t112 = -t154 * t188 + t130 + (-pkin(5) - t196) * t155;
t105 = t114 * t151 - t125 * t155;
t103 = t115 * pkin(5) + t108;
t101 = -t106 * t150 + t113 * t154;
t98 = -t115 * qJ(6) + t100;
t97 = t126 * pkin(5) - t116 * qJ(6) + t99;
t1 = [MDP(1) + (t101 ^ 2 + t102 ^ 2 + t105 ^ 2) * MDP(27); (-t113 * t148 - t125 * t184) * MDP(10) + (-t114 * t148 + t125 * t185) * MDP(11) + (t105 * t184 + t113 * t126) * MDP(17) + (t106 * t184 + t113 * t127) * MDP(18) + (t101 * t126 + t105 * t115) * MDP(24) + (-t102 * t126 + t105 * t116) * MDP(25) + (-t101 * t116 - t102 * t115) * MDP(26) + (t101 * t97 + t102 * t98 + t105 * t103) * MDP(27) + (t157 * MDP(3) - t153 * MDP(4)) * t147; t142 * t152 ^ 2 * MDP(5) + t127 ^ 2 * MDP(12) + (t103 ^ 2 + t97 ^ 2 + t98 ^ 2) * MDP(27) + MDP(2) + (t181 * t203 + t183) * t148 + (-0.2e1 * t115 * MDP(20) + t177) * t116 + ((MDP(8) * t148 - t173) * t203 + (0.2e1 * MDP(6) * t152 + t180) * t142) * t156 + (0.2e1 * MDP(15) * t184 - 0.2e1 * t174 + t175 + 0.2e1 * t176 - 0.2e1 * t178) * t126 + 0.2e1 * (t128 * t148 + t142 * t197) * MDP(10) + 0.2e1 * (-t129 * t148 - t142 * t198) * MDP(11) + 0.2e1 * (-t110 * t184 + t122 * t126) * MDP(17) + 0.2e1 * (t111 * t184 + t122 * t127) * MDP(18) + (t108 * t115 + t99 * t126) * t201 + (-t100 * t126 + t108 * t116) * t200 + (-t98 * t115 - t97 * t116) * t199; -t114 * MDP(11) + (t101 * t112 + t102 * t117 + t105 * t132) * MDP(27) + (-t101 * MDP(24) + t179) * t155 + ((-t101 * t154 - t102 * t150) * MDP(26) + t161 * t105) * t151 + (-t155 * MDP(17) + t151 * MDP(18) - MDP(10)) * t113; t183 + t128 * MDP(10) - t129 * MDP(11) - t127 * t190 + (-t112 * t116 - t117 * t115) * MDP(26) + (t103 * t132 + t97 * t112 + t98 * t117) * MDP(27) + (t156 * MDP(8) + t181) * t146 + (t163 - t191) * t126 + (t174 - t122 * MDP(17) + (-MDP(15) + t189) * t184 - t158) * t155 + (t127 * MDP(12) + t122 * MDP(18) + t116 * t170 + (-t115 * t154 - t116 * t150) * MDP(20) + (pkin(10) * t115 + t187) * MDP(24) + (pkin(10) * t116 + t186) * MDP(25) + (-t150 * t98 - t154 * t97) * MDP(26) + t165 * t184 + t160 * t126) * t151; MDP(9) + (t112 ^ 2 + t117 ^ 2 + t132 ^ 2) * MDP(27) + (t169 + 0.2e1 * t191) * t155 + (t145 * MDP(19) + MDP(12) - 0.2e1 * t166) * t144 + (-t119 * t155 + t144 * t196) * t201 + (t120 * t155 + t144 * t195) * t200 + (-0.2e1 * t190 + (-t112 * t154 - t117 * t150) * t199 - 0.2e1 * t160 * t155) * t151; -t106 * MDP(18) + (-t101 * t150 + t102 * t154) * MDP(26) + (t101 * t134 - t102 * t135) * MDP(27) + (-MDP(17) - t162 + t172) * t105; t173 - t146 * t180 + t110 * MDP(17) - t111 * MDP(18) + t150 * t177 + (-t150 * t115 + t116 * t154) * MDP(20) + (-pkin(4) * t115 - t186) * MDP(24) + (-pkin(4) * t116 + t187) * MDP(25) + (t135 * t115 - t134 * t116 - t97 * t150 + t98 * t154) * MDP(26) + (t103 * t141 + t97 * t134 - t98 * t135) * MDP(27) + t202 * t126; (-t112 * t150 + t117 * t154) * MDP(26) + (t112 * t134 - t117 * t135 + t132 * t141) * MDP(27) + (-t189 - t202) * t155 + (t150 * t170 + (-t143 + t145) * MDP(20) + (-pkin(4) * t150 - t195) * MDP(24) + (-pkin(4) * t154 + t196) * MDP(25) + (-t134 * t154 + t135 * t150) * MDP(26) - t165) * t151; MDP(16) + t143 * MDP(19) + 0.2e1 * t166 + (-t134 * t150 - t135 * t154) * t199 + (t134 ^ 2 + t135 ^ 2 + t141 ^ 2) * MDP(27) + 0.2e1 * t162 * pkin(4); -t179 + (MDP(24) + t192) * t101; (-t116 * MDP(26) + t97 * MDP(27)) * pkin(5) + t158; t112 * t192 - t169 + (t164 * t154 - t171) * t151 + t163; t134 * t192 + (-MDP(25) * pkin(11) + MDP(22)) * t154 + (-MDP(24) * pkin(11) + t164) * t150; MDP(27) * pkin(5) ^ 2 + MDP(23); t105 * MDP(27); t103 * MDP(27); t132 * MDP(27); t172; 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
