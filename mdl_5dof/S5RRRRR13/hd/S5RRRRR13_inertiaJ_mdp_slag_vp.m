% Calculate joint inertia matrix for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR13_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR13_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRRR13_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:16
% EndTime: 2024-09-27 17:32:17
% DurationCPUTime: 0.38s
% Computational Cost: add. (768->159), mult. (1657->180), div. (0->0), fcn. (1533->10), ass. (0->91)
t181 = sin(pkin(5));
t179 = t181 ^ 2;
t184 = sin(qJ(4));
t214 = t179 * t184;
t188 = cos(qJ(4));
t211 = t181 * t188;
t212 = t181 * t184;
t226 = MDP(12) * t212 + MDP(13) * t211;
t225 = MDP(14) + MDP(21);
t223 = 2 * MDP(15);
t222 = 2 * MDP(16);
t221 = 2 * MDP(22);
t220 = 2 * MDP(23);
t186 = sin(qJ(2));
t219 = pkin(1) * t186;
t218 = pkin(4) * t188;
t217 = pkin(10) * t181;
t182 = cos(pkin(5));
t177 = t182 * pkin(4);
t190 = cos(qJ(2));
t171 = t190 * pkin(1) + pkin(2);
t189 = cos(qJ(3));
t161 = t189 * t171;
t185 = sin(qJ(3));
t150 = -t185 * t219 + t161;
t216 = t150 * MDP(8);
t151 = -t185 * t171 - t189 * t219;
t215 = t151 * MDP(9);
t213 = t179 * t188;
t210 = t182 * t184;
t209 = t182 * t188;
t208 = t185 * MDP(9);
t176 = t181 * pkin(9);
t141 = -t151 + t176;
t148 = pkin(3) + t150;
t119 = t188 * t141 + t148 * t210;
t167 = pkin(10) * t211;
t116 = t119 + t167;
t187 = cos(qJ(5));
t207 = t187 * t116;
t160 = t185 * pkin(2) + t176;
t178 = t189 * pkin(2);
t170 = t178 + pkin(3);
t135 = t188 * t160 + t170 * t210;
t127 = t135 + t167;
t206 = t187 * t127;
t153 = pkin(3) * t210 + pkin(9) * t211;
t140 = t153 + t167;
t205 = t187 * t140;
t204 = t190 * MDP(5);
t139 = t148 * t209;
t112 = t139 + t177 + (-t141 - t217) * t184;
t183 = sin(qJ(5));
t106 = t187 * t112 - t183 * t116;
t132 = (-t148 - t218) * t181;
t146 = t183 * t212 - t187 * t211;
t203 = t106 * t182 + t132 * t146;
t156 = t170 * t209;
t124 = t156 + t177 + (-t160 - t217) * t184;
t109 = t187 * t124 - t183 * t127;
t149 = (-t170 - t218) * t181;
t202 = t109 * t182 + t149 * t146;
t169 = pkin(3) * t209;
t133 = t169 + t177 + (-pkin(9) - pkin(10)) * t212;
t114 = t187 * t133 - t183 * t140;
t154 = (-pkin(3) - t218) * t181;
t201 = t114 * t182 + t154 * t146;
t118 = -t184 * t141 + t139;
t200 = t118 * t182 + t148 * t213;
t134 = -t184 * t160 + t156;
t199 = t134 * t182 + t170 * t213;
t152 = -pkin(9) * t212 + t169;
t198 = pkin(3) * t213 + t152 * t182;
t142 = t146 * MDP(20);
t147 = (t183 * t188 + t184 * t187) * t181;
t144 = t147 * MDP(19);
t197 = t182 * MDP(21) - t142 + t144;
t196 = MDP(16) * t214;
t195 = t182 * MDP(14) + t197 + t226;
t194 = (t189 * MDP(8) - t208) * pkin(2);
t193 = (MDP(22) * t187 - MDP(23) * t183) * pkin(4);
t192 = MDP(7) + (MDP(10) * t214 + 0.2e1 * MDP(11) * t213) * t184 + (MDP(17) * t147 - 0.2e1 * MDP(18) * t146) * t147 + (t225 * t182 - 0.2e1 * t142 + 0.2e1 * t144 + 0.2e1 * t226) * t182;
t191 = MDP(4) + t192;
t168 = t187 * t177;
t129 = t154 * t147;
t126 = t149 * t147;
t121 = t132 * t147;
t115 = t183 * t133 + t205;
t110 = t183 * t124 + t206;
t107 = t183 * t112 + t207;
t1 = [MDP(1) + t191 + 0.2e1 * t216 + 0.2e1 * t215 + (-t119 * t182 - t148 * t214) * t222 + t203 * t221 + (-t107 * t182 + t121) * t220 + t200 * t223 + 0.2e1 * (-t186 * MDP(6) + t204) * pkin(1); (-pkin(2) - t171) * t208 + t191 + (t199 + t200) * MDP(15) + (t202 + t203) * MDP(22) + (t204 + (-MDP(8) * t185 - MDP(9) * t189 - MDP(6)) * t186) * pkin(1) + (t121 + t126) * MDP(23) + ((-t119 - t135) * MDP(16) + (-t107 - t110) * MDP(23)) * t182 + (t161 + t178) * MDP(8) + (-t148 - t170) * t196; t191 + t199 * t223 + (-t135 * t182 - t170 * t214) * t222 + t202 * t221 + (-t110 * t182 + t126) * t220 + 0.2e1 * t194; t192 + (t198 + t200) * MDP(15) + t216 + t215 + (t201 + t203) * MDP(22) + ((-t119 - t153) * MDP(16) + (-t107 - t115) * MDP(23)) * t182 + (t121 + t129) * MDP(23) + (-pkin(3) - t148) * t196; t192 + (t198 + t199) * MDP(15) + (t201 + t202) * MDP(22) + (t126 + t129) * MDP(23) + ((-t135 - t153) * MDP(16) + (-t110 - t115) * MDP(23)) * t182 + t194 + (-pkin(3) - t170) * t196; t192 + t198 * t223 + (-pkin(3) * t214 - t153 * t182) * t222 + t201 * t221 + (-t115 * t182 + t129) * t220; t118 * MDP(15) - t119 * MDP(16) + (t106 + t168) * MDP(22) + (-t207 + (-t112 - t177) * t183) * MDP(23) + t195; t134 * MDP(15) - t135 * MDP(16) + (t109 + t168) * MDP(22) + (-t206 + (-t124 - t177) * t183) * MDP(23) + t195; t152 * MDP(15) - t153 * MDP(16) + (t114 + t168) * MDP(22) + (-t205 + (-t133 - t177) * t183) * MDP(23) + t195; 0.2e1 * t193 + t225; t106 * MDP(22) - t107 * MDP(23) + t197; t109 * MDP(22) - t110 * MDP(23) + t197; t114 * MDP(22) - t115 * MDP(23) + t197; MDP(21) + t193; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
