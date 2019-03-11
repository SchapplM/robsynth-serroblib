% Calculate joint inertia matrix for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPPRR4_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:05
% EndTime: 2019-03-09 09:05:08
% DurationCPUTime: 0.76s
% Computational Cost: add. (1025->197), mult. (2443->286), div. (0->0), fcn. (2628->10), ass. (0->92)
t146 = sin(qJ(6));
t149 = cos(qJ(6));
t198 = MDP(29) * t146 + MDP(30) * t149;
t197 = 2 * MDP(14);
t196 = 0.2e1 * MDP(29);
t195 = 0.2e1 * MDP(30);
t194 = pkin(3) + pkin(9);
t151 = cos(qJ(2));
t193 = pkin(1) * t151;
t192 = pkin(8) + qJ(3);
t142 = sin(pkin(11));
t143 = sin(pkin(6));
t144 = cos(pkin(11));
t148 = sin(qJ(2));
t125 = (t142 * t151 + t144 * t148) * t143;
t145 = cos(pkin(6));
t131 = t145 * t193;
t185 = t143 * t148;
t118 = pkin(2) * t145 - t192 * t185 + t131;
t166 = pkin(1) * t145 * t148;
t184 = t143 * t151;
t121 = t192 * t184 + t166;
t108 = t118 * t144 - t142 * t121;
t102 = pkin(4) * t125 - t145 * t194 - t108;
t124 = t142 * t185 - t144 * t184;
t129 = (-pkin(2) * t151 - pkin(1)) * t143;
t156 = -qJ(4) * t125 + t129;
t104 = t124 * t194 + t156;
t147 = sin(qJ(5));
t150 = cos(qJ(5));
t99 = t150 * t102 - t147 * t104;
t97 = -pkin(5) * t125 - t99;
t191 = t146 * t97;
t190 = t149 * t97;
t116 = -t124 * t150 + t145 * t147;
t189 = t116 * t150;
t135 = -pkin(2) * t144 - pkin(3);
t132 = -pkin(9) + t135;
t188 = t132 * t146;
t187 = t132 * t149;
t137 = t143 ^ 2;
t186 = t137 * t148;
t183 = t146 * t147;
t182 = t147 * t149;
t109 = t142 * t118 + t144 * t121;
t180 = MDP(23) * t147;
t179 = MDP(23) * t150;
t117 = t124 * t147 + t145 * t150;
t111 = t117 * t149 + t125 * t146;
t178 = MDP(24) * t111;
t177 = MDP(24) * t149;
t176 = MDP(26) * t111;
t110 = t117 * t146 - t125 * t149;
t173 = t110 * MDP(27);
t172 = t116 * MDP(28);
t171 = t117 * MDP(17);
t170 = t117 * MDP(18);
t169 = t132 * MDP(23);
t133 = pkin(2) * t142 + qJ(4);
t168 = t133 * MDP(15);
t167 = t133 * MDP(22);
t105 = -t145 * qJ(4) - t109;
t165 = MDP(25) * t146 * t149;
t164 = t132 * MDP(22) + MDP(19);
t100 = t102 * t147 + t104 * t150;
t103 = -pkin(4) * t124 - t105;
t163 = MDP(26) * t149 - MDP(27) * t146;
t128 = pkin(5) * t147 - pkin(10) * t150 + t133;
t113 = t128 * t149 - t132 * t183;
t114 = t128 * t146 + t132 * t182;
t162 = MDP(29) * t113 - MDP(30) * t114;
t161 = MDP(29) * t149 - MDP(30) * t146;
t159 = -MDP(18) + t163;
t158 = MDP(22) + t161;
t157 = (MDP(6) * t148 + MDP(7) * t151) * t143;
t101 = pkin(5) * t116 - pkin(10) * t117 + t103;
t98 = pkin(10) * t125 + t100;
t95 = t101 * t149 - t146 * t98;
t96 = t101 * t146 + t149 * t98;
t155 = t95 * MDP(29) - t96 * MDP(30) + t172 - t173;
t154 = t146 * MDP(26) + t149 * MDP(27) - pkin(10) * t198;
t153 = -MDP(20) + t154;
t141 = t150 ^ 2;
t140 = t149 ^ 2;
t139 = t147 ^ 2;
t138 = t146 ^ 2;
t127 = pkin(8) * t184 + t166;
t126 = -pkin(8) * t185 + t131;
t112 = pkin(3) * t124 + t156;
t107 = t111 * t147;
t106 = -pkin(3) * t145 - t108;
t1 = [t125 ^ 2 * MDP(21) + (t108 ^ 2 + t109 ^ 2 + t129 ^ 2) * MDP(12) + (t105 ^ 2 + t106 ^ 2 + t112 ^ 2) * MDP(16) + MDP(1) + (MDP(4) * t148 + 0.2e1 * MDP(5) * t151) * t186 + (0.2e1 * t125 * MDP(19) + t171) * t117 + (-0.2e1 * MDP(25) * t110 + t178) * t111 + (t145 * MDP(8) + 0.2e1 * t157) * t145 + (-0.2e1 * MDP(20) * t125 - 0.2e1 * t170 + t172 - 0.2e1 * t173 + 0.2e1 * t176) * t116 + (t110 * t97 + t116 * t95) * t196 + (t111 * t97 - t116 * t96) * t195 + 0.2e1 * (t103 * t116 + t125 * t99) * MDP(22) + 0.2e1 * (t105 * t124 + t106 * t125) * MDP(13) + 0.2e1 * (-t100 * t125 + t103 * t117) * MDP(23) + 0.2e1 * (-t108 * t125 - t109 * t124) * MDP(11) + (t106 * t145 - t112 * t124) * t197 + 0.2e1 * (-t105 * t145 - t112 * t125) * MDP(15) + 0.2e1 * (-pkin(1) * t186 - t127 * t145) * MDP(10) + 0.2e1 * (t126 * t145 + t137 * t193) * MDP(9); t126 * MDP(9) - t127 * MDP(10) + (-t124 * t133 + t125 * t135) * MDP(13) - t108 * MDP(14) - t105 * MDP(15) + (-t105 * t133 + t106 * t135) * MDP(16) + t133 * t117 * MDP(23) + t107 * MDP(26) + (MDP(8) + (-pkin(3) + t135) * MDP(14) + t168) * t145 + t157 + (t162 + t167) * t116 + (-t170 + t103 * MDP(22) + (-MDP(20) - t169) * t125 + t155) * t147 + ((-t124 * t142 - t125 * t144) * MDP(11) + (t108 * t144 + t109 * t142) * MDP(12)) * pkin(2) + (t103 * MDP(23) + (-t110 * t149 - t111 * t146) * MDP(25) + (-t110 * t132 + t191) * MDP(29) + (-t111 * t132 + t190) * MDP(30) + t111 * t177 + t171 + t164 * t125 + t159 * t116) * t150; MDP(8) + (t133 ^ 2 + t135 ^ 2) * MDP(16) + 0.2e1 * t133 * t179 + t139 * MDP(28) + (t142 ^ 2 + t144 ^ 2) * MDP(12) * pkin(2) ^ 2 + (MDP(24) * t140 + MDP(17) - 0.2e1 * t165) * t141 + 0.2e1 * (t150 * t159 + t167) * t147 + t135 * t197 + 0.2e1 * t168 + (t113 * t147 - t141 * t188) * t196 + (-t114 * t147 - t141 * t187) * t195; t129 * MDP(12) - t124 * MDP(14) + t112 * MDP(16) + (t110 * t147 - t146 * t189) * MDP(29) + (-t149 * t189 + t107) * MDP(30) + (-MDP(22) * t147 - MDP(15) - t179) * t125; 0; MDP(12) + MDP(16); t145 * MDP(14) + t106 * MDP(16) + (-t110 * t150 - t116 * t183) * MDP(29) + (-t111 * t150 - t116 * t182) * MDP(30) + (MDP(22) * t150 + MDP(13) - t180) * t125; t135 * MDP(16) + MDP(14) + t198 * (-t139 - t141); 0; MDP(16); t117 * MDP(19) + t125 * MDP(21) + t99 * MDP(22) - t100 * MDP(23) + t146 * t178 + (-t110 * t146 + t111 * t149) * MDP(25) + (-pkin(5) * t110 - t190) * MDP(29) + (-pkin(5) * t111 + t191) * MDP(30) + t153 * t116; (t153 - t169) * t147 + (t146 * t177 + (-t138 + t140) * MDP(25) + (-pkin(5) * t146 + t187) * MDP(29) + (-pkin(5) * t149 - t188) * MDP(30) + t164) * t150; -t147 * t158 - t179; t150 * t158 - t180; MDP(24) * t138 + 0.2e1 * pkin(5) * t161 + MDP(21) + 0.2e1 * t165; t155 + t176; MDP(28) * t147 + t150 * t163 + t162; -t198 * t150; -t198 * t147; t154; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
