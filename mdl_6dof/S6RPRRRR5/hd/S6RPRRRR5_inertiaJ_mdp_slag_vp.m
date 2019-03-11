% Calculate joint inertia matrix for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR5_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:09:52
% EndTime: 2019-03-09 07:09:54
% DurationCPUTime: 0.67s
% Computational Cost: add. (1303->161), mult. (2467->217), div. (0->0), fcn. (2982->10), ass. (0->94)
t188 = sin(qJ(5));
t192 = cos(qJ(5));
t187 = sin(qJ(6));
t191 = cos(qJ(6));
t163 = t187 * t188 - t191 * t192;
t164 = t187 * t192 + t188 * t191;
t220 = t164 * MDP(31) - t163 * MDP(32);
t208 = t188 * MDP(24) + t192 * MDP(25) + t220;
t201 = MDP(27) * t192 - MDP(28) * t188;
t239 = t163 * MDP(34) + t164 * MDP(35);
t242 = t201 - t239;
t189 = sin(qJ(4));
t174 = pkin(3) * t189 + pkin(9);
t159 = (-pkin(10) - t174) * t188;
t180 = t192 * pkin(10);
t160 = t174 * t192 + t180;
t141 = t159 * t191 - t160 * t187;
t142 = t159 * t187 + t160 * t191;
t241 = t141 * MDP(34) - t142 * MDP(35);
t168 = (-pkin(9) - pkin(10)) * t188;
t169 = pkin(9) * t192 + t180;
t150 = t168 * t191 - t169 * t187;
t151 = t168 * t187 + t169 * t191;
t240 = t150 * MDP(34) - t151 * MDP(35);
t186 = cos(pkin(11));
t172 = -t186 * pkin(2) - pkin(1);
t185 = sin(pkin(11));
t190 = sin(qJ(3));
t193 = cos(qJ(3));
t209 = -t190 * t185 + t193 * t186;
t152 = -pkin(3) * t209 + t172;
t238 = 0.2e1 * t152;
t237 = 0.2e1 * t172;
t236 = -2 * MDP(30);
t235 = cos(qJ(4));
t234 = pkin(1) * MDP(7);
t162 = t185 * t193 + t186 * t190;
t145 = t162 * t189 - t209 * t235;
t233 = pkin(5) * t145;
t232 = t192 * pkin(5);
t230 = pkin(7) + qJ(2);
t146 = t162 * t235 + t189 * t209;
t123 = t145 * pkin(4) - t146 * pkin(9) + t152;
t165 = t230 * t185;
t166 = t230 * t186;
t210 = -t193 * t165 - t166 * t190;
t131 = -pkin(8) * t162 + t210;
t204 = t190 * t165 - t193 * t166;
t132 = pkin(8) * t209 - t204;
t122 = t189 * t131 + t132 * t235;
t227 = t122 * t192;
t109 = t227 + (-pkin(10) * t146 + t123) * t188;
t229 = t109 * t191;
t121 = -t131 * t235 + t189 * t132;
t228 = t121 * t192;
t226 = t146 * t188;
t225 = t146 * t192;
t224 = t185 * MDP(5);
t223 = t186 * MDP(4);
t222 = t188 * t192;
t127 = t163 * t146;
t215 = MDP(29) * t127;
t126 = t164 * t146;
t124 = t126 * MDP(32);
t125 = t127 * MDP(31);
t214 = t146 * MDP(21);
t213 = MDP(26) + MDP(33);
t212 = t145 * MDP(33) - t124 - t125;
t211 = MDP(23) * t222;
t110 = -t122 * t188 + t192 * t123;
t108 = -pkin(10) * t225 + t110 + t233;
t105 = t191 * t108 - t109 * t187;
t175 = -pkin(3) * t235 - pkin(4);
t207 = -pkin(4) * t146 - pkin(9) * t145;
t183 = t188 ^ 2;
t206 = t183 * MDP(22) + MDP(19) + 0.2e1 * t211 + (MDP(29) * t164 + t163 * t236) * t164;
t205 = -t145 * t174 + t146 * t175;
t203 = t209 * MDP(13);
t202 = MDP(24) * t192 - MDP(25) * t188;
t200 = -MDP(27) * t188 - MDP(28) * t192;
t198 = (MDP(34) * t191 - MDP(35) * t187) * pkin(5);
t197 = 0.2e1 * t239;
t196 = (MDP(20) * t235 - t189 * MDP(21)) * pkin(3);
t184 = t192 ^ 2;
t195 = (-t126 * t164 + t127 * t163) * MDP(30) - t164 * t215 - t121 * MDP(20) - t122 * MDP(21) + ((-t183 + t184) * MDP(23) + MDP(22) * t222 + MDP(17)) * t146 + (-MDP(18) + t208) * t145;
t176 = -pkin(4) - t232;
t167 = t175 - t232;
t116 = t121 * t188;
t115 = pkin(5) * t226 + t121;
t114 = t115 * t164;
t113 = t115 * t163;
t111 = t123 * t188 + t227;
t106 = t108 * t187 + t229;
t1 = [t214 * t238 - t203 * t237 + MDP(1) + t213 * t145 ^ 2 - (t126 * t236 - t215) * t127 + (0.2e1 * t223 - 0.2e1 * t224 + t234) * pkin(1) + (MDP(14) * t237 + MDP(8) * t162 + 0.2e1 * MDP(9) * t209) * t162 + (t184 * MDP(22) + MDP(15) - 0.2e1 * t211) * t146 ^ 2 + (MDP(20) * t238 - 0.2e1 * t125 - 0.2e1 * t124 + 0.2e1 * (-MDP(16) + t202) * t146) * t145 + 0.2e1 * (t105 * t145 + t115 * t126) * MDP(34) + 0.2e1 * (-t106 * t145 - t115 * t127) * MDP(35) + 0.2e1 * (t110 * t145 + t121 * t226) * MDP(27) + 0.2e1 * (-t111 * t145 + t121 * t225) * MDP(28) + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t185 ^ 2 + t186 ^ 2) * qJ(2); -t223 + t224 - t234 - t203 + t162 * MDP(14) + t214 + (MDP(20) + t242) * t145; MDP(7); (t188 * t205 - t228) * MDP(27) + t209 * MDP(11) + t204 * MDP(14) + t210 * MDP(13) + t162 * MDP(10) + (t126 * t167 + t141 * t145 + t113) * MDP(34) + (-t127 * t167 - t142 * t145 + t114) * MDP(35) + (t192 * t205 + t116) * MDP(28) + t195; 0; t167 * t197 - 0.2e1 * t175 * t201 + MDP(12) + 0.2e1 * t196 + t206; (t188 * t207 - t228) * MDP(27) + (t126 * t176 + t145 * t150 + t113) * MDP(34) + (-t127 * t176 - t145 * t151 + t114) * MDP(35) + (t192 * t207 + t116) * MDP(28) + t195; 0; t196 + t206 + t201 * (pkin(4) - t175) + t239 * (t167 + t176); 0.2e1 * pkin(4) * t201 + t176 * t197 + t206; t145 * MDP(26) + t110 * MDP(27) - t111 * MDP(28) + (t191 * t233 + t105) * MDP(34) + (-t229 + (-t108 - t233) * t187) * MDP(35) + t202 * t146 + t212; t242; t174 * t200 + t208 + t241; pkin(9) * t200 + t208 + t240; 0.2e1 * t198 + t213; t105 * MDP(34) - t106 * MDP(35) + t212; -t239; t220 + t241; t220 + t240; MDP(33) + t198; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
