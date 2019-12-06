% Calculate joint inertia matrix for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRPRR6_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:36:12
% EndTime: 2019-12-05 18:36:13
% DurationCPUTime: 0.39s
% Computational Cost: add. (453->116), mult. (890->149), div. (0->0), fcn. (788->8), ass. (0->88)
t163 = sin(qJ(2));
t150 = t163 * pkin(1) + qJ(3);
t159 = sin(pkin(9));
t157 = t159 ^ 2;
t160 = cos(pkin(9));
t158 = t160 ^ 2;
t185 = t157 + t158;
t188 = t185 * t150;
t162 = sin(qJ(4));
t165 = cos(qJ(4));
t161 = sin(qJ(5));
t164 = cos(qJ(5));
t132 = -t161 * t162 + t164 * t165;
t133 = t161 * t165 + t164 * t162;
t190 = t132 * MDP(23) - t133 * MDP(24);
t214 = t165 * MDP(16) - t162 * MDP(17) + t190;
t180 = MDP(15) + MDP(22);
t212 = 2 * MDP(9);
t210 = 0.2e1 * t160;
t209 = 0.2e1 * MDP(16);
t208 = 0.2e1 * MDP(17);
t207 = 0.2e1 * MDP(23);
t206 = 0.2e1 * MDP(24);
t205 = pkin(8) * t159;
t166 = cos(qJ(2));
t204 = t166 * pkin(1);
t151 = -pkin(2) - t204;
t203 = pkin(2) - t151;
t202 = pkin(2) * MDP(10);
t123 = t132 * t159;
t197 = t159 * t162;
t148 = pkin(4) * t197;
t125 = t159 * t150 + t148;
t139 = -t160 * pkin(3) - t159 * pkin(7) - pkin(2);
t127 = t139 - t204;
t124 = t165 * t127;
t196 = t159 * t165;
t179 = pkin(8) * t196;
t101 = -t179 + t124 + (-t150 * t162 - pkin(4)) * t160;
t194 = t160 * t165;
t177 = t150 * t194;
t103 = t177 + (t127 - t205) * t162;
t96 = t161 * t101 + t164 * t103;
t201 = t125 * t123 + t96 * t160;
t131 = t159 * qJ(3) + t148;
t130 = t165 * t139;
t108 = -t179 + t130 + (-qJ(3) * t162 - pkin(4)) * t160;
t178 = qJ(3) * t194;
t113 = t178 + (t139 - t205) * t162;
t99 = t161 * t108 + t164 * t113;
t200 = t131 * t123 + t99 * t160;
t199 = t157 * t162;
t198 = t157 * t165;
t155 = t159 * MDP(8);
t195 = t160 * t162;
t95 = t164 * t101 - t161 * t103;
t193 = t95 * MDP(23);
t98 = t164 * t108 - t161 * t113;
t192 = t98 * MDP(23);
t122 = t133 * t159;
t119 = t122 * MDP(21);
t121 = t123 * MDP(20);
t191 = t121 - t119;
t112 = t162 * t127 + t177;
t189 = t112 * t160 + t150 * t198;
t118 = t162 * t139 + t178;
t187 = qJ(3) * t198 + t118 * t160;
t186 = t185 * qJ(3);
t184 = MDP(23) * t164;
t183 = t151 * MDP(10);
t176 = MDP(14) * t197;
t144 = MDP(13) * t196;
t175 = t185 * MDP(10);
t174 = -t160 * MDP(22) + t191;
t173 = MDP(7) * t210 - 0.2e1 * t155;
t172 = (t166 * MDP(5) - t163 * MDP(6)) * pkin(1);
t171 = (-MDP(24) * t161 + t184) * pkin(4);
t170 = t165 ^ 2 * t157 * MDP(11) - 0.2e1 * t162 * MDP(12) * t198 + MDP(4) - 0.2e1 * (t121 + t144) * t160 + (t119 + t176) * t210 + t180 * t158 + (MDP(18) * t123 - 0.2e1 * MDP(19) * t122) * t123;
t169 = t155 + (-MDP(7) - t214) * t160;
t168 = t144 + (-pkin(4) * t184 - t180) * t160 - t176 + t191;
t147 = t161 * pkin(4) * t160;
t145 = qJ(3) * t199;
t134 = t150 * t199;
t117 = -qJ(3) * t195 + t130;
t111 = -t150 * t195 + t124;
t109 = t131 * t122;
t105 = t125 * t122;
t1 = [(-t173 + t183) * t151 + t188 * t212 + (-t111 * t160 + t134) * t209 + t189 * t208 + (-t95 * t160 + t105) * t207 + t201 * t206 + t170 + 0.2e1 * t172 + MDP(1) + t150 ^ 2 * t175; (t187 + t189) * MDP(17) + (t186 + t188) * MDP(9) + (t200 + t201) * MDP(24) + t172 + (t134 + t145) * MDP(16) + (t105 + t109) * MDP(23) + (t203 * MDP(7) + (-t111 - t117) * MDP(16) + (-t95 - t98) * MDP(23)) * t160 + (-t151 * pkin(2) + qJ(3) * t188) * MDP(10) + t170 - t203 * t155; (t173 + t202) * pkin(2) + t186 * t212 + (-t117 * t160 + t145) * t209 + t187 * t208 + (-t98 * t160 + t109) * t207 + t200 * t206 + t170 + qJ(3) ^ 2 * t175; t169 + t183; t169 - t202; MDP(10); t111 * MDP(16) - t112 * MDP(17) + t193 + (t147 - t96) * MDP(24) + t168; t117 * MDP(16) - t118 * MDP(17) + t192 + (t147 - t99) * MDP(24) + t168; t214; 0.2e1 * t171 + t180; -t96 * MDP(24) + t174 + t193; -t99 * MDP(24) + t174 + t192; t190; MDP(22) + t171; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
