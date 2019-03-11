% Calculate joint inertia matrix for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP2_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:37:25
% EndTime: 2019-03-09 16:37:27
% DurationCPUTime: 0.80s
% Computational Cost: add. (1595->177), mult. (2858->251), div. (0->0), fcn. (3250->8), ass. (0->90)
t161 = sin(pkin(10));
t151 = pkin(3) * t161 + pkin(9);
t163 = sin(qJ(5));
t159 = t163 ^ 2;
t166 = cos(qJ(5));
t160 = t166 ^ 2;
t207 = t159 + t160;
t208 = t207 * t151;
t224 = t163 * MDP(22) + t166 * MDP(23);
t222 = t163 * MDP(26);
t168 = cos(qJ(2));
t154 = -t168 * pkin(2) - pkin(1);
t221 = 0.2e1 * t154;
t220 = 2 * MDP(28);
t219 = pkin(7) + pkin(8);
t164 = sin(qJ(3));
t218 = pkin(2) * t164;
t162 = cos(pkin(10));
t217 = pkin(3) * t162;
t165 = sin(qJ(2));
t167 = cos(qJ(3));
t140 = t164 * t168 + t165 * t167;
t180 = t164 * t165 - t167 * t168;
t123 = t140 * t161 + t162 * t180;
t216 = pkin(5) * t123;
t215 = qJ(6) * t123;
t153 = pkin(2) * t167 + pkin(3);
t135 = t161 * t153 + t162 * t218;
t133 = pkin(9) + t135;
t214 = t123 * t133;
t213 = t123 * t151;
t212 = t163 * t166;
t124 = t162 * t140 - t161 * t180;
t130 = t180 * pkin(3) + t154;
t112 = t123 * pkin(4) - t124 * pkin(9) + t130;
t145 = t219 * t165;
t146 = t219 * t168;
t181 = t164 * t145 - t167 * t146;
t119 = -t180 * qJ(4) - t181;
t182 = -t145 * t167 - t146 * t164;
t172 = -qJ(4) * t140 + t182;
t115 = t162 * t119 + t161 * t172;
t106 = t163 * t112 + t166 * t115;
t134 = t153 * t162 - t161 * t218;
t190 = -pkin(5) * t166 - qJ(6) * t163;
t175 = -pkin(4) + t190;
t127 = -t134 + t175;
t137 = t175 - t217;
t211 = -t127 - t137;
t210 = t207 * t133;
t132 = -pkin(4) - t134;
t152 = -pkin(4) - t217;
t209 = t132 + t152;
t206 = MDP(28) * t124;
t205 = MDP(30) * t127;
t204 = MDP(30) * t133;
t203 = MDP(30) * t137;
t202 = MDP(30) * t151;
t201 = MDP(30) * t163;
t113 = t119 * t161 - t162 * t172;
t189 = pkin(5) * t163 - qJ(6) * t166;
t107 = t189 * t124 + t113;
t200 = t107 * MDP(29);
t199 = t123 * MDP(24);
t198 = -MDP(26) + MDP(29);
t197 = -t189 * MDP(28) + t224;
t196 = MDP(21) * t212;
t195 = t159 * MDP(20) + MDP(15) + 0.2e1 * t196;
t194 = -MDP(30) * pkin(5) - MDP(27);
t193 = -t166 * t112 + t115 * t163;
t192 = t207 * MDP(30);
t191 = MDP(25) - t194;
t188 = MDP(30) * qJ(6) + t198;
t103 = t106 + t215;
t104 = t193 - t216;
t187 = t103 * t163 - t104 * t166;
t186 = -t124 * t127 + t214;
t185 = t124 * t132 - t214;
t184 = -t124 * t137 + t213;
t183 = t124 * t152 - t213;
t179 = t166 * MDP(22) - t163 * MDP(23);
t178 = -t193 * MDP(25) - t106 * MDP(26);
t177 = -t113 * MDP(25) - t107 * MDP(27);
t176 = -0.2e1 * t166 * MDP(27) - 0.2e1 * t163 * MDP(29);
t174 = (MDP(16) * t167 - MDP(17) * t164) * pkin(2);
t173 = -0.2e1 * t166 * MDP(25) + 0.2e1 * t222;
t171 = t113 * t222 + (t103 * t166 + t104 * t163) * MDP(28) + t182 * MDP(16) + t181 * MDP(17) - t180 * MDP(14) + t140 * MDP(13) + ((-t159 + t160) * MDP(21) + MDP(20) * t212) * t124 + t224 * t123;
t170 = -t191 * t163 + t188 * t166;
t156 = t163 * MDP(28);
t1 = [MDP(1) + t180 * MDP(16) * t221 + (t113 ^ 2 + t115 ^ 2 + t130 ^ 2) * MDP(19) + (t103 ^ 2 + t104 ^ 2 + t107 ^ 2) * MDP(30) + (MDP(20) * t160 - 0.2e1 * t196) * t124 ^ 2 + (MDP(11) * t140 - 0.2e1 * t180 * MDP(12) + MDP(17) * t221) * t140 + t199 * t123 + 0.2e1 * (-t115 * MDP(18) - t104 * MDP(27) + t103 * MDP(29) + t178) * t123 + (MDP(4) * t165 + 0.2e1 * t168 * MDP(5)) * t165 + 0.2e1 * (-MDP(10) * t165 + MDP(9) * t168) * pkin(1) + 0.2e1 * (-t187 * MDP(28) + (t163 * MDP(27) - t166 * MDP(29)) * t107 + (t163 * MDP(25) + t166 * MDP(26) + MDP(18)) * t113 + t179 * t123) * t124; (t185 * MDP(26) + t186 * MDP(29) + t103 * t204 + t177) * t166 + (-t123 * t135 - t124 * t134) * MDP(18) + (-t113 * t134 + t115 * t135) * MDP(19) + t165 * MDP(6) + t168 * MDP(7) + t107 * t205 + (-MDP(10) * t168 - MDP(9) * t165) * pkin(7) + t171 + (t185 * MDP(25) - t186 * MDP(27) + t104 * t204 - t200) * t163; MDP(8) + (t134 ^ 2 + t135 ^ 2) * MDP(19) + t210 * t220 + t132 * t173 + t133 ^ 2 * t192 + (t176 + t205) * t127 + 0.2e1 * t174 + t195; ((-t123 * t161 - t124 * t162) * MDP(18) + (-t113 * t162 + t115 * t161) * MDP(19)) * pkin(3) + (t183 * MDP(26) + t184 * MDP(29) + t103 * t202 + t177) * t166 + (t183 * MDP(25) - t184 * MDP(27) + t104 * t202 - t200) * t163 + t107 * t203 + t171; (t208 + t210) * MDP(28) + (t127 * t137 + t133 * t208) * MDP(30) + (t134 * t162 + t135 * t161) * MDP(19) * pkin(3) + t174 + (-t209 * MDP(25) + t211 * MDP(27)) * t166 + (t209 * MDP(26) + t211 * MDP(29)) * t163 + t195; t208 * t220 + (t161 ^ 2 + t162 ^ 2) * MDP(19) * pkin(3) ^ 2 + t152 * t173 + t151 ^ 2 * t192 + (t176 + t203) * t137 + t195; t130 * MDP(19) + t187 * MDP(30) - t207 * t206 + ((MDP(25) + MDP(27)) * t166 + t198 * t163) * t123; 0; 0; MDP(19) + t192; t199 + (-t193 + 0.2e1 * t216) * MDP(27) + (t106 + 0.2e1 * t215) * MDP(29) + (-pkin(5) * t104 + qJ(6) * t103) * MDP(30) + (t190 * MDP(28) + t179) * t124 + t178; t170 * t133 + t197; t170 * t151 + t197; t188 * t163 + t191 * t166; MDP(24) + 0.2e1 * pkin(5) * MDP(27) + 0.2e1 * qJ(6) * MDP(29) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(30); -t123 * MDP(27) + MDP(30) * t104 + t166 * t206; t133 * t201 + t156; t151 * t201 + t156; -t166 * MDP(30); t194; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
