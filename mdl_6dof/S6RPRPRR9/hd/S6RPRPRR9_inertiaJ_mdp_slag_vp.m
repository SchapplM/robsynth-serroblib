% Calculate joint inertia matrix for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRPRR9_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:05:09
% EndTime: 2019-03-09 04:05:12
% DurationCPUTime: 1.05s
% Computational Cost: add. (2629->215), mult. (7200->345), div. (0->0), fcn. (8392->14), ass. (0->105)
t185 = sin(qJ(6));
t188 = cos(qJ(6));
t196 = MDP(29) * t185 + MDP(30) * t188;
t193 = t185 * MDP(26) + t188 * MDP(27) - t196 * pkin(11);
t234 = 0.2e1 * MDP(29);
t233 = 0.2e1 * MDP(30);
t180 = sin(pkin(6));
t173 = t180 ^ 2;
t232 = pkin(1) * t173;
t184 = cos(pkin(6));
t231 = pkin(1) * t184;
t178 = sin(pkin(12));
t182 = cos(pkin(12));
t222 = t180 * t182;
t163 = qJ(2) * t222 + t178 * t231;
t179 = sin(pkin(7));
t183 = cos(pkin(7));
t221 = t182 * t183;
t145 = (t179 * t184 + t180 * t221) * pkin(9) + t163;
t168 = t182 * t231;
t225 = t178 * t180;
t148 = pkin(2) * t184 + t168 + (-pkin(9) * t183 - qJ(2)) * t225;
t155 = (-pkin(9) * t178 * t179 - pkin(2) * t182 - pkin(1)) * t180;
t187 = sin(qJ(3));
t190 = cos(qJ(3));
t220 = t183 * t190;
t223 = t179 * t190;
t132 = -t145 * t187 + t148 * t220 + t155 * t223;
t224 = t179 * t187;
t147 = t184 * t224 + (t178 * t190 + t187 * t221) * t180;
t160 = t179 * t222 - t184 * t183;
t128 = -pkin(3) * t160 - qJ(4) * t147 + t132;
t133 = t145 * t190 + (t148 * t183 + t155 * t179) * t187;
t146 = -t184 * t223 + (t178 * t187 - t182 * t220) * t180;
t131 = -qJ(4) * t146 + t133;
t177 = sin(pkin(13));
t181 = cos(pkin(13));
t123 = t177 * t128 + t181 * t131;
t121 = -pkin(10) * t160 + t123;
t140 = -t148 * t179 + t183 * t155;
t134 = pkin(3) * t146 + t140;
t138 = t181 * t146 + t147 * t177;
t139 = -t146 * t177 + t147 * t181;
t124 = pkin(4) * t138 - pkin(10) * t139 + t134;
t186 = sin(qJ(5));
t189 = cos(qJ(5));
t117 = -t121 * t186 + t124 * t189;
t115 = -pkin(5) * t138 - t117;
t230 = t115 * t185;
t229 = t115 * t188;
t171 = pkin(3) * t177 + pkin(10);
t228 = t171 * t185;
t227 = t171 * t188;
t226 = t171 * t189;
t219 = t185 * t186;
t218 = t186 * t188;
t172 = -pkin(3) * t181 - pkin(4);
t217 = MDP(22) * t172;
t216 = MDP(23) * t186;
t136 = t139 * t189 - t160 * t186;
t126 = t136 * t188 + t138 * t185;
t215 = MDP(24) * t126;
t214 = MDP(24) * t188;
t213 = MDP(28) * t189;
t125 = t136 * t185 - t188 * t138;
t212 = t125 * MDP(27);
t211 = t126 * MDP(26);
t135 = t139 * t186 + t160 * t189;
t210 = t135 * MDP(28);
t209 = t136 * MDP(17);
t208 = t136 * MDP(18);
t207 = t146 * MDP(11);
t206 = t160 * MDP(12);
t205 = t172 * MDP(23);
t204 = MDP(25) * t185 * t188;
t122 = t128 * t181 - t177 * t131;
t203 = -t171 * MDP(22) + MDP(19);
t202 = -t171 * MDP(23) + MDP(20);
t118 = t121 * t189 + t124 * t186;
t201 = MDP(22) * t189 - t216;
t200 = MDP(26) * t188 - MDP(27) * t185;
t164 = -pkin(5) * t189 - pkin(11) * t186 + t172;
t149 = t164 * t188 - t185 * t226;
t150 = t164 * t185 + t188 * t226;
t198 = MDP(29) * t149 - MDP(30) * t150;
t197 = MDP(29) * t188 - MDP(30) * t185;
t120 = pkin(4) * t160 - t122;
t195 = -MDP(18) + t200;
t194 = MDP(22) + t197;
t116 = pkin(11) * t138 + t118;
t119 = pkin(5) * t135 - pkin(11) * t136 + t120;
t113 = -t116 * t185 + t119 * t188;
t114 = t116 * t188 + t119 * t185;
t192 = t113 * MDP(29) - t114 * MDP(30) + t210 + t211 - t212;
t176 = t188 ^ 2;
t175 = t186 ^ 2;
t174 = t185 ^ 2;
t162 = -qJ(2) * t225 + t168;
t159 = (t177 * t190 + t181 * t187) * t179;
t157 = t177 * t224 - t181 * t223;
t152 = t159 * t189 + t183 * t186;
t151 = t159 * t186 - t189 * t183;
t142 = t152 * t188 + t157 * t185;
t141 = -t152 * t185 + t157 * t188;
t1 = [(pkin(1) ^ 2 * t173 + t162 ^ 2 + t163 ^ 2) * MDP(7) + t138 ^ 2 * MDP(21) + (t122 ^ 2 + t123 ^ 2 + t134 ^ 2) * MDP(16) + MDP(1) + (t206 + 0.2e1 * t207) * t160 + (0.2e1 * t138 * MDP(19) + t209) * t136 + (-0.2e1 * MDP(25) * t125 + t215) * t126 + (-0.2e1 * MDP(10) * t160 + MDP(8) * t147 - 0.2e1 * MDP(9) * t146) * t147 + (-0.2e1 * MDP(20) * t138 - 0.2e1 * t208 + t210 + 0.2e1 * t211 - 0.2e1 * t212) * t135 + 0.2e1 * (-t163 * t184 - t178 * t232) * MDP(5) + 0.2e1 * (t162 * t184 + t182 * t232) * MDP(4) + 0.2e1 * (t133 * t160 + t140 * t147) * MDP(14) + 0.2e1 * (-t132 * t160 + t140 * t146) * MDP(13) + 0.2e1 * (-t122 * t139 - t123 * t138) * MDP(15) + (t113 * t135 + t115 * t125) * t234 + (-t114 * t135 + t115 * t126) * t233 + 0.2e1 * (t117 * t138 + t120 * t135) * MDP(22) + 0.2e1 * (-t118 * t138 + t120 * t136) * MDP(23) + 0.2e1 * (-t162 * t178 + t163 * t182) * MDP(6) * t180; (t146 * t183 - t160 * t223) * MDP(13) + (t147 * t183 + t160 * t224) * MDP(14) + (-t138 * t159 + t139 * t157) * MDP(15) + (-t122 * t157 + t123 * t159 + t134 * t183) * MDP(16) + (t135 * t157 - t138 * t151) * MDP(22) + (t136 * t157 - t138 * t152) * MDP(23) + (t125 * t151 + t135 * t141) * MDP(29) + (t126 * t151 - t135 * t142) * MDP(30) + (-MDP(4) * t182 + MDP(5) * t178 - MDP(7) * pkin(1)) * t180; MDP(7) + (t157 ^ 2 + t159 ^ 2 + t183 ^ 2) * MDP(16); t136 * t205 + t147 * MDP(10) - t207 - t206 + t132 * MDP(13) - t133 * MDP(14) + (t198 + t217) * t135 + ((-t138 * t177 - t139 * t181) * MDP(15) + (t122 * t181 + t123 * t177) * MDP(16)) * pkin(3) + (-t120 * MDP(22) + t202 * t138 - t192 + t208) * t189 + (t209 + t120 * MDP(23) + t126 * t214 + (-t125 * t188 - t126 * t185) * MDP(25) + (t125 * t171 + t230) * MDP(29) + (t126 * t171 + t229) * MDP(30) + t203 * t138 + t195 * t135) * t186; (-t141 * t189 + t151 * t219) * MDP(29) + (t142 * t189 + t151 * t218) * MDP(30) + (MDP(13) * t190 - MDP(14) * t187) * t179 - t201 * t157 + (-t157 * t181 + t159 * t177) * MDP(16) * pkin(3); MDP(12) + (t177 ^ 2 + t181 ^ 2) * MDP(16) * pkin(3) ^ 2 + (t213 - 0.2e1 * t217) * t189 + (MDP(24) * t176 + MDP(17) - 0.2e1 * t204) * t175 + (-t149 * t189 + t175 * t228) * t234 + (t150 * t189 + t175 * t227) * t233 + 0.2e1 * (-t189 * t195 + t205) * t186; t134 * MDP(16) + (-t125 * t189 - t135 * t219) * MDP(29) + (-t126 * t189 - t135 * t218) * MDP(30) + t201 * t138; t183 * MDP(16); 0; MDP(16); t136 * MDP(19) + t138 * MDP(21) + t117 * MDP(22) - t118 * MDP(23) + t185 * t215 + (-t125 * t185 + t126 * t188) * MDP(25) + (-pkin(5) * t125 - t229) * MDP(29) + (-pkin(5) * t126 + t230) * MDP(30) + (-MDP(20) + t193) * t135; -MDP(23) * t152 - t194 * t151; (-t193 + t202) * t189 + (t185 * t214 + (-t174 + t176) * MDP(25) + (-pkin(5) * t185 - t227) * MDP(29) + (-pkin(5) * t188 + t228) * MDP(30) + t203) * t186; t194 * t189 - t216; MDP(24) * t174 + 0.2e1 * pkin(5) * t197 + MDP(21) + 0.2e1 * t204; t192; MDP(29) * t141 - MDP(30) * t142; t200 * t186 + t198 - t213; -t196 * t186; t193; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
