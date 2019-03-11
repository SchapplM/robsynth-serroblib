% Calculate joint inertia matrix for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR1_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:23:08
% EndTime: 2019-03-09 15:23:10
% DurationCPUTime: 0.74s
% Computational Cost: add. (1609->171), mult. (2944->244), div. (0->0), fcn. (3427->10), ass. (0->91)
t185 = sin(pkin(10));
t175 = pkin(3) * t185 + qJ(5);
t184 = sin(pkin(11));
t186 = cos(pkin(11));
t220 = t184 ^ 2 + t186 ^ 2;
t221 = t220 * t175;
t188 = sin(qJ(6));
t191 = cos(qJ(6));
t163 = t184 * t188 - t186 * t191;
t164 = t184 * t191 + t186 * t188;
t222 = MDP(26) * t164 - MDP(27) * t163;
t200 = MDP(29) * t163 + MDP(30) * t164;
t236 = MDP(20) * t186 - MDP(21) * t184;
t199 = -t236 + t200;
t189 = sin(qJ(3));
t190 = sin(qJ(2));
t192 = cos(qJ(3));
t193 = cos(qJ(2));
t205 = t189 * t190 - t192 * t193;
t230 = pkin(7) + pkin(8);
t170 = t230 * t190;
t171 = t230 * t193;
t206 = t189 * t170 - t192 * t171;
t134 = -qJ(4) * t205 - t206;
t187 = cos(pkin(10));
t165 = t189 * t193 + t190 * t192;
t207 = -t170 * t192 - t171 * t189;
t196 = -qJ(4) * t165 + t207;
t124 = t134 * t185 - t187 * t196;
t235 = t124 ^ 2;
t137 = t165 * t185 + t187 * t205;
t234 = 0.2e1 * t137;
t179 = -t193 * pkin(2) - pkin(1);
t233 = 0.2e1 * t179;
t232 = 2 * MDP(22);
t231 = -2 * MDP(25);
t229 = pkin(2) * t189;
t228 = t186 * pkin(5);
t181 = t186 * pkin(9);
t227 = t124 * t186;
t138 = t165 * t187 - t185 * t205;
t226 = t138 * t184;
t147 = pkin(3) * t205 + t179;
t123 = t137 * pkin(4) - t138 * qJ(5) + t147;
t126 = t134 * t187 + t185 * t196;
t114 = t123 * t184 + t126 * t186;
t178 = pkin(2) * t192 + pkin(3);
t152 = t178 * t185 + t187 * t229;
t149 = qJ(5) + t152;
t224 = t220 * t149;
t128 = t163 * t138;
t219 = MDP(24) * t128;
t127 = t164 * t138;
t218 = t127 * MDP(27);
t217 = t137 * MDP(28);
t151 = t178 * t187 - t185 * t229;
t150 = -pkin(4) - t151;
t216 = t150 * MDP(23);
t177 = -pkin(3) * t187 - pkin(4);
t215 = t177 * MDP(23);
t213 = MDP(15) + (MDP(24) * t164 + t163 * t231) * t164;
t113 = t123 * t186 - t126 * t184;
t212 = t220 * MDP(23);
t211 = t113 * t186 + t114 * t184;
t210 = -t113 * t184 + t114 * t186;
t209 = -t137 * t149 + t138 * t150;
t208 = -t137 * t175 + t138 * t177;
t204 = 0.2e1 * t236;
t203 = MDP(20) * t184 + MDP(21) * t186;
t110 = pkin(5) * t137 - t138 * t181 + t113;
t111 = -pkin(9) * t226 + t114;
t202 = (t110 * t191 - t111 * t188) * MDP(29) - (t110 * t188 + t111 * t191) * MDP(30);
t201 = t127 * MDP(29) - t128 * MDP(30);
t198 = (MDP(16) * t192 - MDP(17) * t189) * pkin(2);
t197 = 0.2e1 * t200;
t195 = t210 * MDP(22) + (-t127 * t164 + t128 * t163) * MDP(25) - t164 * t219 + t207 * MDP(16) + t206 * MDP(17) - t205 * MDP(14) + t165 * MDP(13) + t222 * t137;
t166 = t177 - t228;
t155 = t175 * t186 + t181;
t154 = (-pkin(9) - t175) * t184;
t144 = t150 - t228;
t143 = t149 * t186 + t181;
t142 = (-pkin(9) - t149) * t184;
t136 = t154 * t188 + t155 * t191;
t135 = t154 * t191 - t155 * t188;
t130 = t142 * t188 + t143 * t191;
t129 = t142 * t191 - t143 * t188;
t122 = t124 * t184;
t117 = pkin(5) * t226 + t124;
t116 = t117 * t164;
t115 = t117 * t163;
t1 = [MDP(1) + t190 ^ 2 * MDP(4) + 0.2e1 * t190 * t193 * MDP(5) + t205 * MDP(16) * t233 + (t126 ^ 2 + t147 ^ 2 + t235) * MDP(19) + (t113 ^ 2 + t114 ^ 2 + t235) * MDP(23) + (t217 - 0.2e1 * t218) * t137 + 0.2e1 * (-MDP(10) * t190 + MDP(9) * t193) * pkin(1) + (MDP(11) * t165 - 0.2e1 * MDP(12) * t205 + MDP(17) * t233) * t165 - (MDP(26) * t234 + t127 * t231 - t219) * t128 + 0.2e1 * t201 * t117 + (-MDP(18) * t126 + MDP(20) * t113 - MDP(21) * t114 + t202) * t234 + 0.2e1 * (-t211 * MDP(22) + (MDP(18) + t203) * t124) * t138; (-MDP(10) * t193 - MDP(9) * t190) * pkin(7) + t195 + (t186 * t209 + t122) * MDP(21) + (t127 * t144 + t129 * t137 + t115) * MDP(29) + (-t128 * t144 - t130 * t137 + t116) * MDP(30) + (-t137 * t152 - t138 * t151) * MDP(18) + (-t124 * t151 + t126 * t152) * MDP(19) + t190 * MDP(6) + t193 * MDP(7) + (t124 * t150 + t149 * t210) * MDP(23) + (t184 * t209 - t227) * MDP(20); MDP(8) + (t151 ^ 2 + t152 ^ 2) * MDP(19) + t224 * t232 + t149 ^ 2 * t212 + (-t204 + t216) * t150 + t144 * t197 + 0.2e1 * t198 + t213; ((-t137 * t185 - t138 * t187) * MDP(18) + (-t124 * t187 + t126 * t185) * MDP(19)) * pkin(3) + t195 + (t186 * t208 + t122) * MDP(21) + (t184 * t208 - t227) * MDP(20) + (t127 * t166 + t135 * t137 + t115) * MDP(29) + (-t128 * t166 - t136 * t137 + t116) * MDP(30) + (t124 * t177 + t175 * t210) * MDP(23); (t221 + t224) * MDP(22) + (t149 * t221 + t150 * t177) * MDP(23) + (t151 * t187 + t152 * t185) * MDP(19) * pkin(3) + t198 + t213 + t200 * (t144 + t166) - t236 * (t150 + t177); t221 * t232 + (t185 ^ 2 + t187 ^ 2) * MDP(19) * pkin(3) ^ 2 + t175 ^ 2 * t212 + (-t204 + t215) * t177 + t166 * t197 + t213; -MDP(22) * t138 * t220 + t147 * MDP(19) + MDP(23) * t211 - t137 * t199; 0; 0; MDP(19) + t212; t124 * MDP(23) + t138 * t203 + t201; t199 + t216; t199 + t215; 0; MDP(23); -MDP(26) * t128 + t202 + t217 - t218; MDP(29) * t129 - MDP(30) * t130 + t222; MDP(29) * t135 - MDP(30) * t136 + t222; -t200; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
