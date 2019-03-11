% Calculate joint inertia matrix for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR8_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:05:33
% EndTime: 2019-03-09 14:05:36
% DurationCPUTime: 0.95s
% Computational Cost: add. (1768->196), mult. (3623->280), div. (0->0), fcn. (4119->10), ass. (0->99)
t182 = sin(pkin(11));
t183 = cos(pkin(11));
t186 = sin(qJ(4));
t190 = cos(qJ(4));
t166 = t182 * t190 + t183 * t186;
t187 = sin(qJ(2));
t157 = t166 * t187;
t165 = t182 * t186 - t190 * t183;
t158 = t165 * t187;
t185 = sin(qJ(5));
t189 = cos(qJ(5));
t138 = t189 * t157 - t158 * t185;
t139 = -t157 * t185 - t158 * t189;
t244 = t139 * MDP(24) - t138 * MDP(25);
t227 = pkin(8) + qJ(3);
t169 = t227 * t182;
t170 = t227 * t183;
t151 = -t190 * t169 - t170 * t186;
t152 = -t169 * t186 + t170 * t190;
t140 = -pkin(9) * t166 + t151;
t141 = -pkin(9) * t165 + t152;
t123 = t189 * t140 - t141 * t185;
t124 = t140 * t185 + t141 * t189;
t146 = t189 * t165 + t166 * t185;
t147 = -t165 * t185 + t166 * t189;
t113 = -pkin(10) * t147 + t123;
t114 = -pkin(10) * t146 + t124;
t184 = sin(qJ(6));
t188 = cos(qJ(6));
t129 = t188 * t146 + t147 * t184;
t130 = -t146 * t184 + t147 * t188;
t204 = t130 * MDP(31) - t129 * MDP(32) + (t113 * t188 - t114 * t184) * MDP(34) - (t113 * t184 + t114 * t188) * MDP(35);
t195 = t147 * MDP(24) - t146 * MDP(25) + t123 * MDP(27) - t124 * MDP(28) + t204;
t243 = t166 * MDP(17) - t165 * MDP(18) + t151 * MDP(20) - t152 * MDP(21) + t195;
t214 = MDP(35) * t184;
t197 = (MDP(34) * t188 - t214) * pkin(5);
t207 = t189 * MDP(27);
t242 = (-MDP(28) * t185 + t207) * pkin(4);
t241 = (MDP(14) * qJ(3));
t191 = cos(qJ(2));
t168 = -pkin(2) * t191 - qJ(3) * t187 - pkin(1);
t163 = t183 * t168;
t229 = pkin(7) * t182;
t148 = -pkin(8) * t183 * t187 + t163 + (-pkin(3) - t229) * t191;
t228 = pkin(7) * t191;
t154 = t182 * t168 + t183 * t228;
t225 = t182 * t187;
t150 = -pkin(8) * t225 + t154;
t131 = t190 * t148 - t150 * t186;
t122 = -pkin(4) * t191 + pkin(9) * t158 + t131;
t132 = t148 * t186 + t150 * t190;
t125 = -pkin(9) * t157 + t132;
t111 = t189 * t122 - t125 * t185;
t112 = t122 * t185 + t125 * t189;
t109 = -pkin(5) * t191 - pkin(10) * t139 + t111;
t110 = -pkin(10) * t138 + t112;
t102 = t188 * t109 - t110 * t184;
t103 = t109 * t184 + t110 * t188;
t118 = t188 * t138 + t139 * t184;
t119 = -t138 * t184 + t139 * t188;
t224 = t119 * MDP(31) - t118 * MDP(32);
t239 = t102 * MDP(34) - t103 * MDP(35) + t224;
t240 = t111 * MDP(27) - t112 * MDP(28) + t239 + t244;
t237 = 2 * MDP(13);
t236 = -2 * MDP(16);
t235 = 0.2e1 * MDP(21);
t234 = -2 * MDP(23);
t233 = 0.2e1 * MDP(28);
t232 = -2 * MDP(30);
t231 = 0.2e1 * MDP(35);
t230 = pkin(4) * t185;
t226 = pkin(2) * MDP(14);
t167 = pkin(3) * t225 + t187 * pkin(7);
t222 = MDP(11) * t183;
t221 = MDP(12) * t182;
t220 = MDP(15) * t166;
t219 = MDP(20) * t165;
t218 = MDP(22) * t147;
t217 = MDP(27) * t146;
t216 = MDP(29) * t130;
t215 = MDP(34) * t129;
t176 = pkin(4) * t189 + pkin(5);
t171 = t188 * t176;
t209 = (-t184 * t230 + t171) * MDP(34);
t208 = (t176 * t184 + t188 * t230) * MDP(35);
t206 = MDP(26) + MDP(33);
t175 = -pkin(3) * t183 - pkin(2);
t205 = MDP(19) + t206;
t149 = pkin(4) * t157 + t167;
t156 = pkin(4) * t165 + t175;
t202 = MDP(11) * t182 + MDP(12) * t183;
t201 = -t158 * MDP(17) - t157 * MDP(18);
t198 = MDP(33) - t208 + t209;
t196 = t221 - t222 - t226;
t180 = t187 ^ 2;
t153 = -t182 * t228 + t163;
t133 = pkin(5) * t146 + t156;
t126 = pkin(5) * t138 + t149;
t1 = [(pkin(7) ^ 2 * t180 + t153 ^ 2 + t154 ^ 2) * MDP(14) + t180 * MDP(4) - 0.2e1 * pkin(1) * t187 * MDP(10) + MDP(1) - (-MDP(15) * t158 + t157 * t236) * t158 + (MDP(22) * t139 + t138 * t234) * t139 + (MDP(29) * t119 + t118 * t232) * t119 + t205 * t191 ^ 2 + 0.2e1 * (t187 * MDP(5) + pkin(1) * MDP(9) - t201 - t224 - t244) * t191 + 0.2e1 * (-t153 * t191 + t180 * t229) * MDP(11) + 0.2e1 * (pkin(7) * t180 * t183 + t154 * t191) * MDP(12) + 0.2e1 * (-t131 * t191 + t157 * t167) * MDP(20) + (t132 * t191 - t158 * t167) * t235 + 0.2e1 * (-t111 * t191 + t138 * t149) * MDP(27) + (t112 * t191 + t139 * t149) * t233 + 0.2e1 * (-t102 * t191 + t118 * t126) * MDP(34) + (t103 * t191 + t119 * t126) * t231 + (-t153 * t183 - t154 * t182) * t187 * t237; (-t118 * t130 - t119 * t129) * MDP(30) + (-t138 * t147 - t139 * t146) * MDP(23) + (-t157 * t166 + t158 * t165) * MDP(16) + (t118 * t133 + t126 * t129) * MDP(34) + (t119 * t133 + t126 * t130) * MDP(35) + (t138 * t156 + t146 * t149) * MDP(27) + (t139 * t156 + t147 * t149) * MDP(28) + (t157 * t175 + t165 * t167) * MDP(20) + (-t158 * t175 + t166 * t167) * MDP(21) - t158 * t220 + t139 * t218 + t119 * t216 + (MDP(6) - t202 * pkin(2) + (-MDP(9) + t196) * pkin(7)) * t187 + (-pkin(7) * MDP(10) + t202 * qJ(3) + MDP(7) - t243) * t191 + (MDP(13) + t241) * (-t153 * t182 + t154 * t183); 0.2e1 * t175 * t219 + 0.2e1 * t156 * t217 + 0.2e1 * t133 * t215 + MDP(8) + (-0.2e1 * t221 + 0.2e1 * t222 + t226) * pkin(2) + (t165 * t236 + t175 * t235 + t220) * t166 + (t146 * t234 + t156 * t233 + t218) * t147 + (t129 * t232 + t133 * t231 + t216) * t130 + (t237 + t241) * (t182 ^ 2 + t183 ^ 2) * qJ(3); t157 * MDP(20) - t158 * MDP(21) + t138 * MDP(27) + t139 * MDP(28) + t118 * MDP(34) + t119 * MDP(35) + (MDP(14) * pkin(7) + t202) * t187; MDP(21) * t166 + MDP(28) * t147 + MDP(35) * t130 + t196 + t215 + t217 + t219; MDP(14); t131 * MDP(20) - t132 * MDP(21) + (-MDP(19) - MDP(26) - t198 - t242) * t191 + t201 + t240; t243; 0; t205 - 0.2e1 * t208 + 0.2e1 * t209 + 0.2e1 * t242; (-t206 - t197) * t191 + t240; t195; 0; (pkin(5) * t188 + t171) * MDP(34) + (-pkin(5) - t176) * t214 + (t207 + (-MDP(34) * t184 - MDP(35) * t188 - MDP(28)) * t185) * pkin(4) + t206; 0.2e1 * t197 + t206; -t191 * MDP(33) + t239; t204; 0; t198; MDP(33) + t197; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
