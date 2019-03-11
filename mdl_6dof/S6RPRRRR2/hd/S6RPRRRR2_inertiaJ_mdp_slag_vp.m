% Calculate joint inertia matrix for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRR2_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:58:32
% EndTime: 2019-03-09 06:58:33
% DurationCPUTime: 0.59s
% Computational Cost: add. (796->150), mult. (1471->206), div. (0->0), fcn. (1623->10), ass. (0->86)
t180 = sin(qJ(5));
t184 = cos(qJ(5));
t179 = sin(qJ(6));
t183 = cos(qJ(6));
t155 = t179 * t180 - t183 * t184;
t157 = t179 * t184 + t180 * t183;
t208 = t157 * MDP(28) - t155 * MDP(29);
t197 = t180 * MDP(21) + t184 * MDP(22) + t208;
t223 = MDP(31) * t155 + MDP(32) * t157;
t192 = MDP(24) * t184 - MDP(25) * t180;
t177 = sin(pkin(11));
t164 = pkin(1) * t177 + pkin(7);
t226 = pkin(8) + t164;
t181 = sin(qJ(4));
t168 = pkin(3) * t181 + pkin(9);
t151 = (-pkin(10) - t168) * t180;
t174 = t184 * pkin(10);
t152 = t168 * t184 + t174;
t125 = t151 * t183 - t152 * t179;
t126 = t151 * t179 + t152 * t183;
t225 = t125 * MDP(31) - t126 * MDP(32);
t161 = (-pkin(9) - pkin(10)) * t180;
t162 = pkin(9) * t184 + t174;
t135 = t161 * t183 - t162 * t179;
t136 = t161 * t179 + t162 * t183;
t224 = t135 * MDP(31) - t136 * MDP(32);
t178 = cos(pkin(11));
t165 = -t178 * pkin(1) - pkin(2);
t220 = cos(qJ(3));
t159 = -t220 * pkin(3) + t165;
t222 = 0.2e1 * t159;
t221 = -2 * MDP(27);
t219 = cos(qJ(4));
t182 = sin(qJ(3));
t156 = t181 * t182 - t219 * t220;
t218 = pkin(5) * t156;
t217 = t184 * pkin(5);
t158 = t181 * t220 + t219 * t182;
t116 = t156 * pkin(4) - t158 * pkin(9) + t159;
t149 = t226 * t182;
t150 = t226 * t220;
t124 = -t181 * t149 + t219 * t150;
t213 = t124 * t184;
t103 = t213 + (-pkin(10) * t158 + t116) * t180;
t215 = t103 * t183;
t123 = t219 * t149 + t181 * t150;
t214 = t123 * t184;
t212 = t158 * t180;
t211 = t158 * t184;
t210 = t180 * t184;
t119 = t157 * t158;
t120 = t155 * t158;
t209 = -t119 * MDP(31) + t120 * MDP(32);
t204 = MDP(26) * t120;
t112 = t119 * MDP(29);
t113 = t120 * MDP(28);
t148 = t158 * MDP(18);
t201 = MDP(23) + MDP(30);
t200 = t156 * MDP(30) - t112 - t113;
t199 = MDP(20) * t210;
t198 = t220 * MDP(10);
t104 = t184 * t116 - t124 * t180;
t101 = -pkin(10) * t211 + t104 + t218;
t98 = t183 * t101 - t103 * t179;
t169 = -t219 * pkin(3) - pkin(4);
t196 = -pkin(4) * t158 - pkin(9) * t156;
t175 = t180 ^ 2;
t195 = t175 * MDP(19) + MDP(16) + 0.2e1 * t199 + (MDP(26) * t157 + t155 * t221) * t157;
t194 = -t156 * t168 + t158 * t169;
t193 = MDP(21) * t184 - MDP(22) * t180;
t191 = -MDP(24) * t180 - MDP(25) * t184;
t190 = (MDP(31) * t183 - MDP(32) * t179) * pkin(5);
t189 = 0.2e1 * t223;
t188 = (t219 * MDP(17) - t181 * MDP(18)) * pkin(3);
t187 = -t148 + (-MDP(17) - t192 + t223) * t156;
t176 = t184 ^ 2;
t186 = (-t119 * t157 + t120 * t155) * MDP(27) - t157 * t204 - t123 * MDP(17) - t124 * MDP(18) + ((-t175 + t176) * MDP(20) + MDP(19) * t210 + MDP(14)) * t158 + (-MDP(15) + t197) * t156;
t170 = -pkin(4) - t217;
t160 = t169 - t217;
t115 = t123 * t180;
t108 = pkin(5) * t212 + t123;
t107 = t108 * t157;
t106 = t108 * t155;
t105 = t116 * t180 + t213;
t99 = t101 * t179 + t215;
t1 = [-0.2e1 * t165 * t198 + t148 * t222 + MDP(1) + (t177 ^ 2 + t178 ^ 2) * MDP(4) * pkin(1) ^ 2 + t201 * t156 ^ 2 - (t119 * t221 - t204) * t120 + (0.2e1 * t165 * MDP(11) + MDP(5) * t182 + 0.2e1 * t220 * MDP(6)) * t182 + (t176 * MDP(19) + MDP(12) - 0.2e1 * t199) * t158 ^ 2 + (MDP(17) * t222 - 0.2e1 * t113 - 0.2e1 * t112 + 0.2e1 * (-MDP(13) + t193) * t158) * t156 + 0.2e1 * (t104 * t156 + t123 * t212) * MDP(24) + 0.2e1 * (-t105 * t156 + t123 * t211) * MDP(25) + 0.2e1 * (t108 * t119 + t156 * t98) * MDP(31) + 0.2e1 * (-t108 * t120 - t156 * t99) * MDP(32); 0; MDP(4); (t194 * t184 + t115) * MDP(25) + (t194 * t180 - t214) * MDP(24) + t220 * MDP(8) + (-t182 * MDP(10) - t220 * MDP(11)) * t164 + (t119 * t160 + t125 * t156 + t106) * MDP(31) + (-t120 * t160 - t126 * t156 + t107) * MDP(32) + t182 * MDP(7) + t186; -t182 * MDP(11) + t187 + t198; t160 * t189 - 0.2e1 * t169 * t192 + MDP(9) + 0.2e1 * t188 + t195; (t196 * t180 - t214) * MDP(24) + (t119 * t170 + t135 * t156 + t106) * MDP(31) + (-t120 * t170 - t136 * t156 + t107) * MDP(32) + (t196 * t184 + t115) * MDP(25) + t186; t187; t188 + t195 + t192 * (pkin(4) - t169) + t223 * (t160 + t170); 0.2e1 * pkin(4) * t192 + t170 * t189 + t195; t156 * MDP(23) + t104 * MDP(24) - t105 * MDP(25) + (t183 * t218 + t98) * MDP(31) + (-t215 + (-t101 - t218) * t179) * MDP(32) + t193 * t158 + t200; t191 * t158 + t209; t191 * t168 + t197 + t225; t191 * pkin(9) + t197 + t224; 0.2e1 * t190 + t201; t98 * MDP(31) - t99 * MDP(32) + t200; t209; t208 + t225; t208 + t224; MDP(30) + t190; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
