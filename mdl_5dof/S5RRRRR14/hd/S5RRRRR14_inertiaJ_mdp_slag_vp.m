% Calculate joint inertia matrix for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR14_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR14_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR14_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:38
% EndTime: 2024-09-27 18:43:39
% DurationCPUTime: 0.64s
% Computational Cost: add. (976->163), mult. (2253->196), div. (0->0), fcn. (2284->10), ass. (0->105)
t228 = MDP(18) + MDP(25);
t264 = MDP(11) + t228;
t211 = sin(pkin(5));
t209 = t211 ^ 2;
t215 = sin(qJ(3));
t248 = t209 * t215;
t219 = cos(qJ(3));
t245 = t211 * t219;
t246 = t211 * t215;
t263 = MDP(10) * t245 + MDP(9) * t246;
t260 = 2 * MDP(12);
t259 = 2 * MDP(13);
t258 = 2 * MDP(19);
t257 = 2 * MDP(20);
t256 = 2 * MDP(26);
t255 = 2 * MDP(27);
t214 = sin(qJ(4));
t254 = pkin(3) * t214;
t253 = pkin(3) * t219;
t218 = cos(qJ(4));
t175 = t214 * t246 - t218 * t245;
t252 = t175 * pkin(4);
t212 = cos(pkin(5));
t208 = t212 * pkin(3);
t207 = t212 * pkin(4);
t217 = cos(qJ(5));
t251 = t217 * pkin(4);
t250 = t218 * pkin(3);
t199 = pkin(4) + t250;
t213 = sin(qJ(5));
t179 = t213 * t199 + t217 * t254;
t249 = t179 * t212;
t247 = t209 * t219;
t244 = t212 * t215;
t243 = t212 * t219;
t220 = cos(qJ(2));
t200 = t220 * pkin(1) + pkin(2);
t184 = t200 * t243;
t216 = sin(qJ(2));
t188 = t216 * pkin(1) + t211 * pkin(8);
t153 = t184 + t208 + (-pkin(9) * t211 - t188) * t215;
t164 = t219 * t188 + t200 * t244;
t195 = pkin(9) * t245;
t156 = t164 + t195;
t240 = t218 * t156;
t138 = t214 * t153 + t240;
t173 = t175 * pkin(10);
t128 = t138 - t173;
t242 = t217 * t128;
t198 = pkin(2) * t243;
t162 = t198 + t208 + (-pkin(8) - pkin(9)) * t246;
t181 = pkin(2) * t244 + pkin(8) * t245;
t167 = t181 + t195;
t239 = t218 * t167;
t142 = t214 * t162 + t239;
t139 = t142 - t173;
t241 = t217 * t139;
t137 = t218 * t153 - t214 * t156;
t176 = (t214 * t219 + t215 * t218) * t211;
t226 = -t176 * pkin(10) + t207;
t127 = t137 + t226;
t121 = t217 * t127 - t213 * t128;
t148 = t217 * t175 + t213 * t176;
t177 = (-t200 - t253) * t211;
t152 = t177 + t252;
t238 = t121 * t212 + t152 * t148;
t141 = t218 * t162 - t214 * t167;
t133 = t141 + t226;
t124 = t217 * t133 - t213 * t139;
t182 = (-pkin(2) - t253) * t211;
t157 = t182 + t252;
t237 = t124 * t212 + t157 * t148;
t236 = t137 * t212 + t177 * t175;
t235 = t141 * t212 + t182 * t175;
t163 = -t215 * t188 + t184;
t234 = t163 * t212 + t200 * t247;
t180 = -pkin(8) * t246 + t198;
t233 = pkin(2) * t247 + t180 * t212;
t145 = t148 * MDP(24);
t149 = -t213 * t175 + t217 * t176;
t147 = t149 * MDP(23);
t169 = t175 * MDP(17);
t171 = t176 * MDP(16);
t189 = t217 * t199;
t178 = -t213 * t254 + t189;
t232 = t178 * MDP(26);
t231 = t179 * MDP(27);
t230 = t213 * MDP(27);
t229 = t218 * MDP(19);
t227 = t212 * MDP(25) - t145 + t147;
t122 = t213 * t127 + t242;
t125 = t213 * t133 + t241;
t225 = t212 * MDP(18) - t169 + t171 + t227;
t224 = (t220 * MDP(5) - t216 * MDP(6)) * pkin(1);
t223 = (t217 * MDP(26) - t230) * pkin(4);
t222 = t212 * MDP(11) + t225 + t263;
t221 = MDP(4) - 0.2e1 * (t145 + t169) * t212 + (MDP(7) * t248 + 0.2e1 * MDP(8) * t247) * t215 + (MDP(14) * t176 - 0.2e1 * MDP(15) * t175) * t176 + (MDP(21) * t149 - 0.2e1 * MDP(22) * t148) * t149 + t264 * t212 ^ 2 + 0.2e1 * (t147 + t171 + t263) * t212;
t197 = t212 * t250;
t196 = t212 * t251;
t172 = t178 * t212;
t159 = t182 * t176;
t155 = t177 * t176;
t136 = t157 * t149;
t132 = t152 * t149;
t1 = [MDP(1) + t221 + t238 * t256 + t236 * t258 + t234 * t260 + (-t122 * t212 + t132) * t255 + (-t138 * t212 + t155) * t257 + (-t164 * t212 - t200 * t248) * t259 + 0.2e1 * t224; t221 + t224 + (t155 + t159) * MDP(20) + (t132 + t136) * MDP(27) + ((-t164 - t181) * MDP(13) + (-t138 - t142) * MDP(20) + (-t122 - t125) * MDP(27)) * t212 + (t233 + t234) * MDP(12) + (t235 + t236) * MDP(19) + (t237 + t238) * MDP(26) + (-pkin(2) - t200) * MDP(13) * t248; t233 * t260 + (-pkin(2) * t248 - t181 * t212) * t259 + t235 * t258 + (-t142 * t212 + t159) * t257 + t237 * t256 + (-t125 * t212 + t136) * t255 + t221; t222 + (-t240 + (-t153 - t208) * t214) * MDP(20) + (-t122 - t249) * MDP(27) + (t121 + t172) * MDP(26) + (t137 + t197) * MDP(19) + t163 * MDP(12) - t164 * MDP(13); t222 + (-t239 + (-t162 - t208) * t214) * MDP(20) + (-t125 - t249) * MDP(27) + (t124 + t172) * MDP(26) + (t141 + t197) * MDP(19) + t180 * MDP(12) - t181 * MDP(13); 0.2e1 * (-t214 * MDP(20) + t229) * pkin(3) + 0.2e1 * t232 - 0.2e1 * t231 + t264; t137 * MDP(19) - t138 * MDP(20) + (t121 + t196) * MDP(26) + (-t242 + (-t127 - t207) * t213) * MDP(27) + t225; t141 * MDP(19) - t142 * MDP(20) + (t124 + t196) * MDP(26) + (-t241 + (-t133 - t207) * t213) * MDP(27) + t225; (t189 + t251) * MDP(26) + (-pkin(4) - t199) * t230 + (t229 + (-MDP(26) * t213 - MDP(27) * t217 - MDP(20)) * t214) * pkin(3) + t228; 0.2e1 * t223 + t228; t121 * MDP(26) - t122 * MDP(27) + t227; t124 * MDP(26) - t125 * MDP(27) + t227; MDP(25) - t231 + t232; MDP(25) + t223; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
