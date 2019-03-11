% Calculate joint inertia matrix for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR4_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:36:02
% EndTime: 2019-03-09 15:36:05
% DurationCPUTime: 0.82s
% Computational Cost: add. (1082->199), mult. (2066->282), div. (0->0), fcn. (2127->8), ass. (0->85)
t160 = sin(pkin(10));
t161 = cos(pkin(10));
t163 = sin(qJ(3));
t166 = cos(qJ(3));
t138 = t160 * t166 + t161 * t163;
t164 = sin(qJ(2));
t134 = t138 * t164;
t200 = t164 * t166;
t202 = t163 * t164;
t135 = -t160 * t202 + t161 * t200;
t162 = sin(qJ(6));
t165 = cos(qJ(6));
t112 = -t165 * t134 + t162 * t135;
t113 = t162 * t134 + t165 * t135;
t176 = -t113 * MDP(26) + t112 * MDP(27);
t167 = cos(qJ(2));
t145 = -t167 * pkin(2) - t164 * pkin(8) - pkin(1);
t139 = t166 * t145;
t206 = pkin(7) * t163;
t121 = -qJ(4) * t200 + t139 + (-pkin(3) - t206) * t167;
t204 = pkin(7) * t167;
t186 = t166 * t204;
t124 = t186 + (-qJ(4) * t164 + t145) * t163;
t107 = t160 * t121 + t161 * t124;
t103 = -t167 * qJ(5) + t107;
t100 = t134 * pkin(9) + t103;
t106 = t161 * t121 - t160 * t124;
t105 = t167 * pkin(4) - t106;
t99 = t167 * pkin(5) - t135 * pkin(9) + t105;
t179 = t162 * t100 - t165 * t99;
t98 = t165 * t100 + t162 * t99;
t215 = t179 * MDP(29) + t98 * MDP(30) + t176;
t203 = -qJ(4) - pkin(8);
t146 = t203 * t166;
t184 = t203 * t163;
t126 = -t160 * t146 - t161 * t184;
t128 = -t161 * t146 + t160 * t184;
t137 = t160 * t163 - t161 * t166;
t114 = t137 * pkin(9) + t128;
t119 = -t165 * t137 + t162 * t138;
t120 = t162 * t137 + t165 * t138;
t173 = -t138 * pkin(9) + t126;
t170 = t120 * MDP(26) - t119 * MDP(27) - (t162 * t114 - t165 * t173) * MDP(29) - (t165 * t114 + t162 * t173) * MDP(30);
t213 = t163 * MDP(13) + t166 * MDP(14) - t126 * MDP(20) + t128 * MDP(22) - (t163 * MDP(16) + t166 * MDP(17)) * pkin(8) - t170;
t212 = 2 * MDP(18);
t211 = 2 * MDP(21);
t210 = 0.2e1 * MDP(22);
t209 = -2 * MDP(25);
t208 = 0.2e1 * MDP(30);
t207 = -pkin(4) - pkin(5);
t205 = pkin(7) * t166;
t201 = t163 * t166;
t143 = pkin(3) * t202 + t164 * pkin(7);
t197 = t113 * MDP(24);
t153 = -t166 * pkin(3) - pkin(2);
t174 = t138 * qJ(5) - t153;
t116 = t137 * pkin(4) - t174;
t195 = t116 * MDP(23);
t194 = t119 * MDP(29);
t149 = t160 * pkin(3) + qJ(5);
t151 = t161 * pkin(3) + pkin(4);
t181 = -pkin(5) - t151;
t193 = (t162 * t149 - t165 * t181) * MDP(29);
t192 = (t165 * t149 + t162 * t181) * MDP(30);
t191 = t137 * MDP(20);
t190 = t138 * MDP(22);
t189 = t151 * MDP(20);
t188 = MDP(15) + MDP(28);
t187 = t126 ^ 2 + t128 ^ 2;
t185 = MDP(12) * t201;
t183 = t126 * t135 - t128 * t134;
t180 = t135 * qJ(5) - t143;
t178 = t166 * MDP(13) - t163 * MDP(14);
t175 = t165 * MDP(29) - t162 * MDP(30);
t172 = MDP(20) + t175;
t171 = -MDP(28) - t192 - t193;
t158 = t166 ^ 2;
t157 = t164 ^ 2;
t156 = t163 ^ 2;
t131 = t163 * t145 + t186;
t130 = -t163 * t204 + t139;
t109 = t207 * t137 + t174;
t108 = t134 * pkin(4) - t180;
t104 = t207 * t134 + t180;
t1 = [(t106 ^ 2 + t107 ^ 2 + t143 ^ 2) * MDP(19) + MDP(1) - 0.2e1 * pkin(1) * t164 * MDP(10) + (t103 ^ 2 + t105 ^ 2 + t108 ^ 2) * MDP(23) + t188 * t167 ^ 2 + (t112 * t209 + t197) * t113 + (t158 * MDP(11) + MDP(4) - 0.2e1 * t185) * t157 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t178) * t164 - t176) * t167 + (-t103 * t134 + t105 * t135) * t211 + (-t106 * t135 - t107 * t134) * t212 + 0.2e1 * (t104 * t112 - t167 * t179) * MDP(29) + (t104 * t113 - t98 * t167) * t208 + 0.2e1 * (t105 * t167 + t108 * t134) * MDP(20) + (-t103 * t167 - t108 * t135) * t210 + 0.2e1 * (-t130 * t167 + t157 * t206) * MDP(16) + 0.2e1 * (t131 * t167 + t157 * t205) * MDP(17); (-t106 * t138 - t107 * t137 + t183) * MDP(18) + (-t106 * t126 + t107 * t128 + t143 * t153) * MDP(19) + (t108 * t137 + t116 * t134) * MDP(20) + (-t103 * t137 + t105 * t138 + t183) * MDP(21) + (-t108 * t138 - t116 * t135) * MDP(22) + (t103 * t128 + t105 * t126 + t108 * t116) * MDP(23) + t120 * t197 + (-t120 * t112 - t113 * t119) * MDP(25) + (t104 * t119 + t109 * t112) * MDP(29) + (t104 * t120 + t109 * t113) * MDP(30) + (-pkin(7) * MDP(10) + MDP(7) - t213) * t167 + (MDP(6) - pkin(7) * MDP(9) + MDP(11) * t201 + (-t156 + t158) * MDP(12) + (-pkin(2) * t163 - t205) * MDP(16) + (-pkin(2) * t166 + t206) * MDP(17)) * t164; MDP(8) + t156 * MDP(11) + 0.2e1 * t185 + (t153 ^ 2 + t187) * MDP(19) + t187 * MDP(23) + 0.2e1 * t109 * t194 + (-0.2e1 * t190 + 0.2e1 * t191 + t195) * t116 + 0.2e1 * (t166 * MDP(16) - t163 * MDP(17)) * pkin(2) + (MDP(24) * t120 + t109 * t208 + t119 * t209) * t120 + (t212 + t211) * (t126 * t138 - t128 * t137); t130 * MDP(16) - t131 * MDP(17) - t105 * MDP(20) + (-t149 * t134 - t151 * t135) * MDP(21) + t107 * MDP(22) + (t103 * t149 - t105 * t151) * MDP(23) + t178 * t164 + (-MDP(15) - t189 + (-qJ(5) - t149) * MDP(22) + t171) * t167 + ((-t134 * t160 - t135 * t161) * MDP(18) + (t106 * t161 + t107 * t160) * MDP(19)) * pkin(3) + t215; (-t149 * t137 - t151 * t138) * MDP(21) + (-t126 * t151 + t128 * t149) * MDP(23) + ((-t137 * t160 - t138 * t161) * MDP(18) + (-t126 * t161 + t128 * t160) * MDP(19)) * pkin(3) + t213; (t149 ^ 2 + t151 ^ 2) * MDP(23) + (t160 ^ 2 + t161 ^ 2) * MDP(19) * pkin(3) ^ 2 + 0.2e1 * t189 + t149 * t210 + 0.2e1 * t193 + 0.2e1 * t192 + t188; t143 * MDP(19) + t134 * MDP(20) - t135 * MDP(22) + t108 * MDP(23) - t112 * MDP(29) - t113 * MDP(30); t153 * MDP(19) - t120 * MDP(30) - t190 + t191 - t194 + t195; 0; MDP(19) + MDP(23); t135 * MDP(21) + t105 * MDP(23) + t172 * t167; t138 * MDP(21) + t126 * MDP(23); -t151 * MDP(23) - t172; 0; MDP(23); t167 * MDP(28) - t215; t170; t171; 0; t175; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
