% Calculate joint inertia matrix for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR8_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:49
% EndTime: 2019-12-31 22:25:52
% DurationCPUTime: 0.56s
% Computational Cost: add. (637->139), mult. (1221->192), div. (0->0), fcn. (1325->8), ass. (0->79)
t161 = sin(qJ(4));
t165 = cos(qJ(4));
t160 = sin(qJ(5));
t164 = cos(qJ(5));
t138 = t160 * t161 - t164 * t165;
t140 = t160 * t165 + t164 * t161;
t187 = t140 * MDP(27) - t138 * MDP(28);
t176 = t161 * MDP(20) + t165 * MDP(21) + t187;
t204 = pkin(6) + pkin(7);
t162 = sin(qJ(3));
t150 = t162 * pkin(2) + pkin(8);
t134 = (-pkin(9) - t150) * t161;
t157 = t165 * pkin(9);
t135 = t165 * t150 + t157;
t109 = t164 * t134 - t160 * t135;
t110 = t160 * t134 + t164 * t135;
t203 = t109 * MDP(30) - t110 * MDP(31);
t143 = (-pkin(8) - pkin(9)) * t161;
t145 = t165 * pkin(8) + t157;
t121 = t164 * t143 - t160 * t145;
t123 = t160 * t143 + t164 * t145;
t202 = t121 * MDP(30) - t123 * MDP(31);
t171 = t165 * MDP(23) - t161 * MDP(24);
t201 = t138 * MDP(30) + t140 * MDP(31);
t198 = cos(qJ(2));
t153 = -t198 * pkin(2) - pkin(1);
t200 = 0.2e1 * t153;
t199 = -2 * MDP(26);
t197 = cos(qJ(3));
t163 = sin(qJ(2));
t139 = t162 * t163 - t197 * t198;
t196 = t139 * pkin(4);
t195 = t165 * pkin(4);
t141 = t162 * t198 + t197 * t163;
t112 = t139 * pkin(3) - t141 * pkin(8) + t153;
t144 = t204 * t163;
t146 = t204 * t198;
t124 = -t162 * t144 + t197 * t146;
t188 = t165 * t124;
t95 = t188 + (-pkin(9) * t141 + t112) * t161;
t193 = t164 * t95;
t122 = t197 * t144 + t162 * t146;
t192 = t122 * t165;
t191 = t141 * t161;
t190 = t141 * t165;
t189 = t161 * t165;
t105 = t138 * t141;
t185 = MDP(25) * t105;
t104 = t140 * t141;
t102 = t104 * MDP(28);
t103 = t105 * MDP(27);
t180 = 0.2e1 * t198;
t179 = MDP(22) + MDP(29);
t178 = t139 * MDP(29) - t102 - t103;
t177 = MDP(19) * t189;
t96 = t165 * t112 - t161 * t124;
t94 = -pkin(9) * t190 + t196 + t96;
t90 = -t160 * t95 + t164 * t94;
t151 = -t197 * pkin(2) - pkin(3);
t175 = -pkin(3) * t141 - pkin(8) * t139;
t158 = t161 ^ 2;
t174 = t158 * MDP(18) + MDP(15) + 0.2e1 * t177 + (MDP(25) * t140 + t138 * t199) * t140;
t173 = -t139 * t150 + t141 * t151;
t172 = t165 * MDP(20) - t161 * MDP(21);
t170 = -MDP(23) * t161 - MDP(24) * t165;
t169 = (MDP(30) * t164 - MDP(31) * t160) * pkin(4);
t168 = 0.2e1 * t201;
t167 = (t197 * MDP(16) - t162 * MDP(17)) * pkin(2);
t159 = t165 ^ 2;
t166 = -t140 * t185 + (-t140 * t104 + t105 * t138) * MDP(26) - t122 * MDP(16) - t124 * MDP(17) + ((-t158 + t159) * MDP(19) + MDP(18) * t189 + MDP(13)) * t141 + (-MDP(14) + t176) * t139;
t152 = -pkin(3) - t195;
t142 = t151 - t195;
t113 = t122 * t161;
t101 = pkin(4) * t191 + t122;
t99 = t101 * t140;
t98 = t101 * t138;
t97 = t161 * t112 + t188;
t91 = t160 * t94 + t193;
t1 = [pkin(1) * MDP(9) * t180 + MDP(1) + t179 * t139 ^ 2 - (t104 * t199 - t185) * t105 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t163 + MDP(5) * t180) * t163 + (MDP(16) * t200 - 0.2e1 * t102 - 0.2e1 * t103) * t139 + 0.2e1 * (t122 * t191 + t96 * t139) * MDP(23) + 0.2e1 * (t122 * t190 - t97 * t139) * MDP(24) + 0.2e1 * (t101 * t104 + t90 * t139) * MDP(30) + 0.2e1 * (-t101 * t105 - t91 * t139) * MDP(31) + (MDP(17) * t200 + 0.2e1 * (-MDP(12) + t172) * t139 + (t159 * MDP(18) + MDP(11) - 0.2e1 * t177) * t141) * t141; t166 + (t173 * t165 + t113) * MDP(24) + (t173 * t161 - t192) * MDP(23) + (t142 * t104 + t109 * t139 + t98) * MDP(30) + (-t142 * t105 - t110 * t139 + t99) * MDP(31) + t163 * MDP(6) + (-t198 * MDP(10) - t163 * MDP(9)) * pkin(6) + t198 * MDP(7); t142 * t168 - 0.2e1 * t151 * t171 + MDP(8) + 0.2e1 * t167 + t174; t166 + (t175 * t165 + t113) * MDP(24) + (t175 * t161 - t192) * MDP(23) + (t152 * t104 + t121 * t139 + t98) * MDP(30) + (-t152 * t105 - t123 * t139 + t99) * MDP(31); t167 + t174 + t171 * (pkin(3) - t151) + t201 * (t142 + t152); 0.2e1 * pkin(3) * t171 + t152 * t168 + t174; t139 * MDP(22) + t96 * MDP(23) - t97 * MDP(24) + (t164 * t196 + t90) * MDP(30) + (-t193 + (-t94 - t196) * t160) * MDP(31) + t172 * t141 + t178; t170 * t150 + t176 + t203; t170 * pkin(8) + t176 + t202; 0.2e1 * t169 + t179; t90 * MDP(30) - t91 * MDP(31) + t178; t187 + t203; t187 + t202; MDP(29) + t169; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
