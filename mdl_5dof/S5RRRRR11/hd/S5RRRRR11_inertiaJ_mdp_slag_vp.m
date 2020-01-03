% Calculate joint inertia matrix for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR11_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:43:45
% EndTime: 2019-12-31 22:43:49
% DurationCPUTime: 1.15s
% Computational Cost: add. (1073->210), mult. (2494->300), div. (0->0), fcn. (2712->10), ass. (0->104)
t174 = sin(qJ(4));
t178 = cos(qJ(4));
t173 = sin(qJ(5));
t177 = cos(qJ(5));
t155 = t173 * t174 - t177 * t178;
t156 = t173 * t178 + t174 * t177;
t228 = pkin(9) + pkin(10);
t159 = t228 * t174;
t160 = t228 * t178;
t194 = t156 * MDP(27) - t155 * MDP(28) + (-t159 * t177 - t160 * t173) * MDP(30) - (-t159 * t173 + t160 * t177) * MDP(31);
t242 = -(t174 * MDP(23) + t178 * MDP(24)) * pkin(9) + t174 * MDP(20) + t178 * MDP(21) + t194;
t175 = sin(qJ(3));
t179 = cos(qJ(3));
t158 = -pkin(3) * t179 - pkin(9) * t175 - pkin(2);
t154 = t178 * t158;
t221 = pkin(10) * t175;
t224 = pkin(8) * t174;
t128 = -t178 * t221 + t154 + (-pkin(4) - t224) * t179;
t222 = pkin(8) * t179;
t198 = t178 * t222;
t131 = t198 + (t158 - t221) * t174;
t116 = t177 * t128 - t131 * t173;
t117 = t128 * t173 + t131 * t177;
t145 = t156 * t175;
t139 = t145 * MDP(28);
t146 = t155 * t175;
t140 = t146 * MDP(27);
t241 = MDP(30) * t116 - MDP(31) * t117 - t139 - t140;
t240 = -MDP(14) + t242;
t137 = -t174 * t222 + t154;
t138 = t158 * t174 + t198;
t239 = t137 * MDP(23) - t138 * MDP(24) + t241;
t237 = -0.2e1 * t175;
t183 = (MDP(30) * t177 - MDP(31) * t173) * pkin(4);
t171 = sin(pkin(5));
t180 = cos(qJ(2));
t213 = t171 * t180;
t172 = cos(pkin(5));
t176 = sin(qJ(2));
t214 = t171 * t176;
t149 = t172 * t175 + t179 * t214;
t129 = t149 * t174 + t178 * t213;
t130 = t149 * t178 - t174 * t213;
t118 = t177 * t129 + t130 * t173;
t119 = -t129 * t173 + t130 * t177;
t236 = t119 * MDP(27) - t118 * MDP(28);
t161 = pkin(7) * t214;
t226 = pkin(1) * t180;
t141 = t161 + (-pkin(2) - t226) * t172;
t148 = -t172 * t179 + t175 * t214;
t121 = t148 * pkin(3) - t149 * pkin(9) + t141;
t199 = pkin(7) * t213;
t227 = pkin(1) * t176;
t142 = t199 + (pkin(8) + t227) * t172;
t143 = (-pkin(2) * t180 - pkin(8) * t176 - pkin(1)) * t171;
t125 = t179 * t142 + t175 * t143;
t123 = -pkin(9) * t213 + t125;
t111 = t178 * t121 - t123 * t174;
t112 = t121 * t174 + t123 * t178;
t234 = -t111 * MDP(23) + t112 * MDP(24);
t233 = 0.2e1 * MDP(23);
t232 = 0.2e1 * MDP(24);
t231 = -2 * MDP(26);
t230 = 0.2e1 * MDP(30);
t229 = 0.2e1 * MDP(31);
t225 = pkin(4) * t148;
t223 = pkin(8) * t178;
t220 = MDP(17) * pkin(2);
t219 = pkin(2) * MDP(16);
t218 = pkin(8) * MDP(17);
t110 = -pkin(10) * t129 + t112;
t217 = t110 * t177;
t124 = -t175 * t142 + t143 * t179;
t122 = pkin(3) * t213 - t124;
t216 = t122 * t174;
t215 = t122 * t178;
t212 = t172 * MDP(8);
t210 = MDP(15) * t180;
t209 = MDP(25) * t146;
t208 = MDP(25) * t156;
t203 = t130 * MDP(18);
t202 = t149 * MDP(13);
t201 = t178 * MDP(18);
t200 = MDP(22) + MDP(29);
t197 = t148 * MDP(29) + t236;
t196 = MDP(19) * t174 * t178;
t195 = pkin(8) * MDP(16) - MDP(13);
t109 = -pkin(10) * t130 + t111 + t225;
t106 = t177 * t109 - t110 * t173;
t193 = t130 * MDP(20) - t129 * MDP(21);
t192 = MDP(20) * t178 - MDP(21) * t174;
t107 = t109 * t173 + t217;
t188 = -t106 * MDP(30) + t107 * MDP(31);
t185 = -MDP(12) + t192;
t182 = t149 * MDP(12) - t193 - t236;
t169 = t178 ^ 2;
t167 = t174 ^ 2;
t166 = t171 ^ 2;
t165 = -pkin(4) * t178 - pkin(3);
t157 = (pkin(4) * t174 + pkin(8)) * t175;
t151 = t172 * t227 + t199;
t150 = t172 * t226 - t161;
t113 = pkin(4) * t129 + t122;
t1 = [t166 * t176 ^ 2 * MDP(4) + t149 ^ 2 * MDP(11) + MDP(1) + (0.2e1 * MDP(6) * t214 + t212) * t172 + t200 * t148 ^ 2 + (-0.2e1 * t129 * MDP(19) + t203) * t130 + (MDP(25) * t119 + t118 * t231) * t119 + (0.2e1 * MDP(5) * t176 + t210) * t166 * t180 + 0.2e1 * (-t124 * t213 + t141 * t148) * MDP(16) + 0.2e1 * (t125 * t213 + t141 * t149) * MDP(17) + 0.2e1 * (t150 * t172 + t166 * t226) * MDP(9) + 0.2e1 * (-t151 * t172 - t166 * t227) * MDP(10) + (t111 * t148 + t122 * t129) * t233 + (t106 * t148 + t113 * t118) * t230 + (-t107 * t148 + t113 * t119) * t229 + (-t112 * t148 + t122 * t130) * t232 + 0.2e1 * (MDP(7) * t172 - t202) * t213 + 0.2e1 * (MDP(14) * t213 - t182) * t148; (-t113 * t146 + t119 * t157) * MDP(31) + (t113 * t145 + t118 * t157) * MDP(30) - t149 * t220 + t212 + t150 * MDP(9) - t151 * MDP(10) + (t118 * t146 - t119 * t145) * MDP(26) - t119 * t209 + (MDP(6) * t176 + MDP(7) * t180) * t171 + (-t219 + t239) * t148 + (-t141 * MDP(16) + (-MDP(14) + t218) * t213 - t200 * t148 + t182 + t188 + t234) * t179 + ((pkin(8) * t129 + t216) * MDP(23) + (pkin(8) * t130 + t215) * MDP(24) + t141 * MDP(17) + t130 * t201 + (-t129 * t178 - t130 * t174) * MDP(19) + t149 * MDP(11) + t195 * t213 + t185 * t148) * t175; t145 * t157 * t230 + t220 * t237 + MDP(8) - (t145 * t231 + t157 * t229 - t209) * t146 + (MDP(18) * t169 + t223 * t232 + t224 * t233 + MDP(11) - 0.2e1 * t196) * t175 ^ 2 + (-t116 * t230 + t117 * t229 - t137 * t233 + t138 * t232 + t200 * t179 + t185 * t237 + 0.2e1 * t139 + 0.2e1 * t140 + 0.2e1 * t219) * t179; t202 - t171 * t210 + t124 * MDP(16) - t125 * MDP(17) + t174 * t203 + (-t129 * t174 + t130 * t178) * MDP(19) + (-pkin(3) * t129 - t215) * MDP(23) + (-pkin(3) * t130 + t216) * MDP(24) + t119 * t208 + (-t118 * t156 - t119 * t155) * MDP(26) + (t113 * t155 + t118 * t165) * MDP(30) + (t113 * t156 + t119 * t165) * MDP(31) + t240 * t148; -t146 * t208 + (-t145 * t156 + t146 * t155) * MDP(26) + (t145 * t165 + t155 * t157) * MDP(30) + (-t146 * t165 + t156 * t157) * MDP(31) + (-t218 - t240) * t179 + (t174 * t201 + (-t167 + t169) * MDP(19) + (-pkin(3) * t174 - t223) * MDP(23) + (-pkin(3) * t178 + t224) * MDP(24) - t195) * t175; 0.2e1 * t196 + t155 * t165 * t230 + MDP(18) * t167 + MDP(15) + 0.2e1 * (MDP(23) * t178 - MDP(24) * t174) * pkin(3) + (t155 * t231 + t165 * t229 + t208) * t156; t148 * MDP(22) + (t177 * t225 + t106) * MDP(30) + (-t217 + (-t109 - t225) * t173) * MDP(31) + t193 + t197 - t234; (-t200 - t183) * t179 + t192 * t175 + t239; t242; 0.2e1 * t183 + t200; -t188 + t197; -t179 * MDP(29) + t241; t194; MDP(29) + t183; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
