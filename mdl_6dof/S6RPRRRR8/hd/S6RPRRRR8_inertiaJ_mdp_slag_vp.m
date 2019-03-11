% Calculate joint inertia matrix for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% MDP [34x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(34,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [34 1]), ...
  'S6RPRRRR8_inertiaJ_mdp_slag_vp: MDP has to be [34x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:21:15
% EndTime: 2019-03-09 07:21:17
% DurationCPUTime: 0.61s
% Computational Cost: add. (785->157), mult. (1397->212), div. (0->0), fcn. (1549->8), ass. (0->84)
t182 = sin(qJ(5));
t186 = cos(qJ(5));
t181 = sin(qJ(6));
t185 = cos(qJ(6));
t157 = t181 * t182 - t185 * t186;
t160 = t181 * t186 + t182 * t185;
t212 = MDP(30) * t160 - MDP(31) * t157;
t200 = MDP(23) * t182 + MDP(24) * t186 + t212;
t233 = MDP(26) * t182 + MDP(27) * t186;
t195 = MDP(26) * t186 - MDP(27) * t182;
t228 = MDP(33) * t157 + MDP(34) * t160;
t183 = sin(qJ(4));
t184 = sin(qJ(3));
t187 = cos(qJ(3));
t225 = cos(qJ(4));
t158 = t183 * t184 - t187 * t225;
t117 = t160 * t158;
t119 = t157 * t158;
t232 = t119 * MDP(30) + t117 * MDP(31);
t171 = pkin(3) * t183 + pkin(9);
t153 = (-pkin(10) - t171) * t182;
t178 = t186 * pkin(10);
t154 = t171 * t186 + t178;
t123 = t153 * t185 - t154 * t181;
t124 = t153 * t181 + t154 * t185;
t231 = MDP(33) * t123 - MDP(34) * t124;
t165 = (-pkin(9) - pkin(10)) * t182;
t166 = pkin(9) * t186 + t178;
t139 = t165 * t185 - t166 * t181;
t140 = t165 * t181 + t166 * t185;
t230 = MDP(33) * t139 - MDP(34) * t140;
t227 = t158 ^ 2;
t159 = t183 * t187 + t184 * t225;
t155 = t159 ^ 2;
t226 = -2 * MDP(29);
t224 = (pkin(1) * MDP(6));
t223 = t159 * pkin(5);
t222 = t186 * pkin(5);
t188 = -pkin(1) - pkin(7);
t220 = -pkin(8) + t188;
t162 = t220 * t184;
t163 = t220 * t187;
t130 = t162 * t183 - t163 * t225;
t219 = t130 * t186;
t218 = t158 * t182;
t217 = t158 * t186;
t216 = t182 * t186;
t168 = pkin(3) * t184 + qJ(2);
t126 = pkin(4) * t159 + pkin(9) * t158 + t168;
t131 = t162 * t225 + t163 * t183;
t214 = t186 * t131;
t105 = t214 + (pkin(10) * t158 + t126) * t182;
t215 = t185 * t105;
t116 = t160 * t159;
t118 = t157 * t159;
t213 = -MDP(33) * t116 + MDP(34) * t118;
t210 = MDP(28) * t119;
t203 = MDP(25) + MDP(32);
t202 = MDP(32) * t159 + t232;
t201 = MDP(22) * t216;
t106 = t126 * t186 - t131 * t182;
t104 = pkin(10) * t217 + t106 + t223;
t100 = t185 * t104 - t105 * t181;
t172 = -pkin(3) * t225 - pkin(4);
t199 = pkin(4) * t158 - pkin(9) * t159;
t179 = t182 ^ 2;
t198 = t179 * MDP(21) + MDP(18) + 0.2e1 * t201 + (MDP(28) * t160 + t157 * t226) * t160;
t197 = -t158 * t172 - t159 * t171;
t196 = -MDP(23) * t186 + MDP(24) * t182;
t193 = (MDP(33) * t185 - MDP(34) * t181) * pkin(5);
t192 = 0.2e1 * t228;
t191 = (MDP(19) * t225 - MDP(20) * t183) * pkin(3);
t190 = -t159 * MDP(20) + (-MDP(19) + t228 - t195) * t158;
t180 = t186 ^ 2;
t189 = (t117 * t160 - t119 * t157) * MDP(29) + t160 * t210 - t130 * MDP(19) - t131 * MDP(20) + ((t179 - t180) * MDP(22) - MDP(21) * t216 - MDP(16)) * t158 + (-MDP(17) + t200) * t159;
t173 = -pkin(4) - t222;
t164 = t172 - t222;
t127 = t130 * t182;
t111 = -pkin(5) * t218 + t130;
t109 = t111 * t160;
t108 = t111 * t157;
t107 = t126 * t182 + t214;
t101 = t104 * t181 + t215;
t1 = [-0.2e1 * t168 * t158 * MDP(20) + MDP(1) + (MDP(7) * t187 - 0.2e1 * MDP(8) * t184) * t187 + t203 * t155 + (-t117 * t226 + t210) * t119 + ((-2 * MDP(4) + t224) * pkin(1)) + (MDP(21) * t180 + MDP(14) - 0.2e1 * t201) * t227 + (0.2e1 * MDP(12) * t184 + 0.2e1 * t187 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * (t168 * MDP(19) + (MDP(15) + t196) * t158 + t232) * t159 + 0.2e1 * (t100 * t159 - t111 * t117) * MDP(33) + 0.2e1 * (-t101 * t159 + t111 * t119) * MDP(34) + 0.2e1 * (t106 * t159 - t130 * t218) * MDP(26) + 0.2e1 * (-t107 * t159 - t130 * t217) * MDP(27); MDP(4) - t224 + (-t116 * t159 - t158 * t117) * MDP(33) + (t118 * t159 + t119 * t158) * MDP(34) + t233 * (-t155 - t227); MDP(6); (-t117 * t164 + t123 * t159 + t108) * MDP(33) + (t119 * t164 - t124 * t159 + t109) * MDP(34) + t189 + (t186 * t197 + t127) * MDP(27) + (t182 * t197 - t219) * MDP(26) + (-t188 * MDP(13) - MDP(10)) * t184 + (MDP(12) * t188 + MDP(9)) * t187; MDP(12) * t187 - t184 * MDP(13) + t190; t164 * t192 - 0.2e1 * t172 * t195 + MDP(11) + 0.2e1 * t191 + t198; (-t117 * t173 + t139 * t159 + t108) * MDP(33) + (t119 * t173 - t140 * t159 + t109) * MDP(34) + t189 + (t186 * t199 + t127) * MDP(27) + (t182 * t199 - t219) * MDP(26); t190; t191 + t198 + t195 * (pkin(4) - t172) + t228 * (t164 + t173); 0.2e1 * pkin(4) * t195 + t173 * t192 + t198; t159 * MDP(25) + t106 * MDP(26) - t107 * MDP(27) + (t185 * t223 + t100) * MDP(33) + (-t215 + (-t104 - t223) * t181) * MDP(34) + t196 * t158 + t202; -t159 * t233 + t213; -t171 * t233 + t200 + t231; -pkin(9) * t233 + t200 + t230; 0.2e1 * t193 + t203; t100 * MDP(33) - t101 * MDP(34) + t202; t213; t212 + t231; t212 + t230; MDP(32) + t193; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
