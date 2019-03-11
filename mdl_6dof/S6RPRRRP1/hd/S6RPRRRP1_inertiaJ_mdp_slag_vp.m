% Calculate joint inertia matrix for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRRP1_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:57:08
% EndTime: 2019-03-09 05:57:10
% DurationCPUTime: 0.61s
% Computational Cost: add. (872->157), mult. (1501->216), div. (0->0), fcn. (1514->8), ass. (0->85)
t150 = sin(qJ(4));
t137 = pkin(3) * t150 + pkin(9);
t149 = sin(qJ(5));
t145 = t149 ^ 2;
t152 = cos(qJ(5));
t146 = t152 ^ 2;
t187 = t145 + t146;
t189 = t187 * t137;
t207 = t149 * MDP(21) + t152 * MDP(22);
t186 = MDP(25) * t149;
t206 = t186 - MDP(17);
t151 = sin(qJ(3));
t200 = cos(qJ(4));
t201 = cos(qJ(3));
t126 = t150 * t201 + t200 * t151;
t205 = 0.2e1 * t126;
t147 = sin(pkin(10));
t134 = pkin(1) * t147 + pkin(7);
t204 = pkin(8) + t134;
t148 = cos(pkin(10));
t135 = -t148 * pkin(1) - pkin(2);
t128 = -t201 * pkin(3) + t135;
t203 = 0.2e1 * t128;
t202 = 2 * MDP(27);
t125 = t150 * t151 - t200 * t201;
t199 = pkin(5) * t125;
t198 = pkin(9) * t125;
t180 = t200 * pkin(3);
t138 = -t180 - pkin(4);
t197 = pkin(4) - t138;
t104 = t125 * pkin(4) - t126 * pkin(9) + t128;
t121 = t204 * t151;
t122 = t204 * t201;
t108 = -t150 * t121 + t200 * t122;
t174 = -t152 * t104 + t108 * t149;
t96 = t174 - t199;
t196 = MDP(29) * t96;
t195 = qJ(6) * t125;
t194 = t125 * t137;
t193 = t126 * t145;
t192 = t149 * t152;
t107 = t200 * t121 + t150 * t122;
t168 = pkin(5) * t149 - qJ(6) * t152;
t99 = t168 * t126 + t107;
t191 = t99 * MDP(28);
t98 = t149 * t104 + t152 * t108;
t169 = -t152 * pkin(5) - t149 * qJ(6);
t129 = -pkin(4) + t169;
t120 = -t180 + t129;
t190 = -t120 - t129;
t188 = t187 * pkin(9);
t185 = MDP(29) * t137;
t184 = MDP(29) * t149;
t183 = t120 * MDP(29);
t119 = t126 * MDP(18);
t182 = t129 * MDP(29);
t181 = t149 * MDP(28);
t179 = -t168 * MDP(27) + t207;
t178 = MDP(20) * t192;
t177 = t145 * MDP(19) + MDP(16) + 0.2e1 * t178;
t176 = t201 * MDP(10);
t175 = -MDP(29) * pkin(5) - MDP(26);
t173 = MDP(29) * t187;
t115 = t146 * t126;
t172 = (t115 + t193) * MDP(27) - t119 + t206 * t125;
t171 = -pkin(4) * t126 - t198;
t170 = t126 * t173;
t167 = -t126 * t129 + t198;
t95 = t98 + t195;
t166 = t149 * t96 + t152 * t95;
t165 = -t174 * MDP(24) - t98 * MDP(25);
t164 = t120 * t126 - t194;
t163 = t126 * t138 - t194;
t162 = -t107 * MDP(24) - t99 * MDP(26);
t161 = t152 * MDP(21) - t149 * MDP(22);
t160 = MDP(24) * t152 - t186;
t159 = -0.2e1 * MDP(26) * t152 - 0.2e1 * t181;
t158 = (-MDP(24) - MDP(26)) * t152 - t181;
t157 = (t200 * MDP(17) - t150 * MDP(18)) * pkin(3);
t156 = -t108 * MDP(18) + t166 * MDP(27) + (t115 - t193) * MDP(20) + (MDP(19) * t192 + MDP(14)) * t126 + t206 * t107 + (-MDP(15) + t207) * t125;
t155 = (MDP(29) * qJ(6) - MDP(25) + MDP(28)) * t152 + (-MDP(24) + t175) * t149;
t140 = t149 * MDP(27);
t124 = t126 ^ 2;
t123 = t125 ^ 2;
t1 = [MDP(1) - 0.2e1 * t135 * t176 + t119 * t203 + t123 * MDP(23) + (t95 ^ 2 + t96 ^ 2 + t99 ^ 2) * MDP(29) + (t147 ^ 2 + t148 ^ 2) * MDP(4) * pkin(1) ^ 2 + (0.2e1 * t135 * MDP(11) + MDP(5) * t151 + 0.2e1 * t201 * MDP(6)) * t151 + (MDP(19) * t146 + MDP(12) - 0.2e1 * t178) * t124 + (MDP(17) * t203 + (-MDP(13) + t161) * t205) * t125 + 0.2e1 * (-MDP(26) * t96 + MDP(28) * t95 + t165) * t125 + ((-t149 * t95 + t152 * t96) * MDP(27) + (t149 * MDP(26) - t152 * MDP(28)) * t99 + (t149 * MDP(24) + t152 * MDP(25)) * t107) * t205; (t125 * t99 + t166 * t126) * MDP(29); MDP(4) + (t187 * t124 + t123) * MDP(29); t156 + (-t151 * MDP(10) - t201 * MDP(11)) * t134 + (t163 * MDP(24) + t164 * MDP(26) + t96 * t185 - t191) * t149 + t201 * MDP(8) + t99 * t183 + t151 * MDP(7) + (t163 * MDP(25) - t164 * MDP(28) + t95 * t185 + t162) * t152; t176 - t151 * MDP(11) + t137 * t170 + (t158 + t183) * t125 + t172; MDP(9) + t189 * t202 + t137 ^ 2 * t173 + (t159 + t183) * t120 + t177 - 0.2e1 * t160 * t138 + 0.2e1 * t157; t99 * t182 + (pkin(9) * t95 * MDP(29) + t171 * MDP(25) + t167 * MDP(28) + t162) * t152 + (t171 * MDP(24) - t167 * MDP(26) + pkin(9) * t196 - t191) * t149 + t156; pkin(9) * t170 + (t158 + t182) * t125 + t172; (t188 + t189) * MDP(27) + (pkin(9) * t189 + t120 * t129) * MDP(29) + t157 + (t197 * MDP(24) + t190 * MDP(26)) * t152 + (-t197 * MDP(25) + t190 * MDP(28)) * t149 + t177; t188 * t202 + pkin(9) ^ 2 * t173 + (t159 + t182) * t129 + 0.2e1 * t160 * pkin(4) + t177; t125 * MDP(23) + (-t174 + 0.2e1 * t199) * MDP(26) + (t98 + 0.2e1 * t195) * MDP(28) + (-pkin(5) * t96 + qJ(6) * t95) * MDP(29) + (t169 * MDP(27) + t161) * t126 + t165; t155 * t126; t155 * t137 + t179; t155 * pkin(9) + t179; MDP(23) + 0.2e1 * pkin(5) * MDP(26) + 0.2e1 * qJ(6) * MDP(28) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(29); MDP(27) * t126 * t152 - MDP(26) * t125 + t196; t126 * t184; t137 * t184 + t140; pkin(9) * t184 + t140; t175; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
