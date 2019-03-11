% Calculate joint inertia matrix for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR5_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:20:49
% EndTime: 2019-03-08 22:20:51
% DurationCPUTime: 0.60s
% Computational Cost: add. (789->160), mult. (1731->243), div. (0->0), fcn. (1950->12), ass. (0->85)
t153 = sin(qJ(6));
t157 = cos(qJ(6));
t167 = (MDP(28) * t157 - MDP(29) * t153) * pkin(5);
t149 = sin(pkin(12));
t193 = pkin(9) + qJ(4);
t138 = t193 * t149;
t151 = cos(pkin(12));
t139 = t193 * t151;
t154 = sin(qJ(5));
t158 = cos(qJ(5));
t120 = -t138 * t158 - t139 * t154;
t121 = -t138 * t154 + t139 * t158;
t134 = t149 * t154 - t151 * t158;
t135 = t149 * t158 + t151 * t154;
t108 = -pkin(10) * t135 + t120;
t109 = -pkin(10) * t134 + t121;
t113 = t134 * t157 + t135 * t153;
t114 = -t134 * t153 + t135 * t157;
t174 = t114 * MDP(25) - t113 * MDP(26) + (t108 * t157 - t109 * t153) * MDP(28) - (t108 * t153 + t109 * t157) * MDP(29);
t204 = t135 * MDP(18) - t134 * MDP(19) + t120 * MDP(21) - t121 * MDP(22) + t174;
t152 = cos(pkin(6));
t155 = sin(qJ(3));
t159 = cos(qJ(3));
t150 = sin(pkin(6));
t156 = sin(qJ(2));
t189 = t150 * t156;
t130 = t152 * t155 + t159 * t189;
t160 = cos(qJ(2));
t188 = t150 * t160;
t117 = -t130 * t149 - t151 * t188;
t118 = t130 * t151 - t149 * t188;
t102 = t117 * t158 - t118 * t154;
t103 = t117 * t154 + t118 * t158;
t192 = (t102 * t157 - t103 * t153) * MDP(28) - (t102 * t153 + t103 * t157) * MDP(29);
t203 = MDP(21) * t102 - MDP(22) * t103 + t192;
t183 = (MDP(15) * qJ(4));
t202 = MDP(14) + t183;
t126 = t135 * t155;
t127 = t134 * t155;
t106 = t126 * t157 - t127 * t153;
t107 = -t126 * t153 - t127 * t157;
t185 = MDP(25) * t107 - MDP(26) * t106;
t137 = -pkin(3) * t159 - qJ(4) * t155 - pkin(2);
t132 = t151 * t137;
t195 = pkin(8) * t149;
t115 = -pkin(9) * t151 * t155 + t132 + (-pkin(4) - t195) * t159;
t194 = pkin(8) * t159;
t123 = t137 * t149 + t151 * t194;
t190 = t149 * t155;
t119 = -pkin(9) * t190 + t123;
t100 = t115 * t158 - t119 * t154;
t96 = -pkin(5) * t159 + pkin(10) * t127 + t100;
t101 = t115 * t154 + t119 * t158;
t99 = -pkin(10) * t126 + t101;
t87 = -t153 * t99 + t157 * t96;
t88 = t153 * t96 + t157 * t99;
t201 = t87 * MDP(28) - MDP(29) * t88 + t185;
t200 = 2 * MDP(14);
t199 = -2 * MDP(17);
t198 = 0.2e1 * MDP(22);
t197 = -2 * MDP(24);
t196 = 0.2e1 * MDP(29);
t191 = pkin(3) * MDP(15);
t136 = pkin(4) * t190 + pkin(8) * t155;
t182 = t113 * MDP(28);
t181 = t114 * MDP(23);
t180 = t134 * MDP(21);
t179 = t135 * MDP(16);
t178 = t149 * MDP(13);
t177 = t151 * MDP(12);
t176 = t155 * MDP(11);
t175 = MDP(20) + MDP(27);
t143 = -pkin(4) * t151 - pkin(3);
t171 = MDP(12) * t149 + MDP(13) * t151;
t170 = -t127 * MDP(18) - t126 * MDP(19);
t166 = -t177 + t178 - t191;
t165 = MDP(15) * pkin(8) + t171;
t163 = t126 * MDP(21) - t127 * MDP(22) + t106 * MDP(28) + t107 * MDP(29);
t162 = MDP(22) * t135 + MDP(29) * t114 + t166 + t180 + t182;
t147 = t155 ^ 2;
t129 = -t152 * t159 + t155 * t189;
t125 = pkin(5) * t134 + t143;
t122 = -t149 * t194 + t132;
t116 = pkin(5) * t126 + t136;
t1 = [MDP(1) + (t117 ^ 2 + t118 ^ 2 + t129 ^ 2) * MDP(15); (t117 * t122 + t118 * t123) * MDP(15) + t163 * t129 + (-MDP(12) * t117 + MDP(13) * t118 - t203) * t159 + ((-t117 * t151 - t118 * t149) * MDP(14) + t165 * t129) * t155 + (-t156 * MDP(4) + (MDP(10) * t159 + MDP(3) - t176) * t160) * t150; MDP(2) + t147 * MDP(5) - 0.2e1 * pkin(2) * t176 + (pkin(8) ^ 2 * t147 + t122 ^ 2 + t123 ^ 2) * MDP(15) + t175 * t159 ^ 2 - (-MDP(16) * t127 + t126 * t199) * t127 + (MDP(23) * t107 + t106 * t197) * t107 + 0.2e1 * (MDP(10) * pkin(2) + MDP(6) * t155 - t170 - t185) * t159 + 0.2e1 * (-t122 * t159 + t147 * t195) * MDP(12) + 0.2e1 * (pkin(8) * t147 * t151 + t123 * t159) * MDP(13) + 0.2e1 * (-t100 * t159 + t126 * t136) * MDP(21) + (t101 * t159 - t127 * t136) * t198 + 0.2e1 * (t106 * t116 - t159 * t87) * MDP(28) + (t107 * t116 + t159 * t88) * t196 + (-t122 * t151 - t123 * t149) * t155 * t200; -t130 * MDP(11) + (-MDP(10) + t162) * t129 + t202 * (-t117 * t149 + t118 * t151); -t127 * t179 + (-t126 * t135 + t127 * t134) * MDP(17) + (t126 * t143 + t134 * t136) * MDP(21) + (-t127 * t143 + t135 * t136) * MDP(22) + t107 * t181 + (-t106 * t114 - t107 * t113) * MDP(24) + (t106 * t125 + t113 * t116) * MDP(28) + (t107 * t125 + t114 * t116) * MDP(29) + (MDP(7) - t171 * pkin(3) + (-MDP(10) + t166) * pkin(8)) * t155 + (-pkin(8) * MDP(11) + qJ(4) * t171 + MDP(8) - t204) * t159 + t202 * (-t122 * t149 + t123 * t151); 0.2e1 * t143 * t180 + 0.2e1 * t125 * t182 + MDP(9) + (0.2e1 * t177 - 0.2e1 * t178 + t191) * pkin(3) + (t134 * t199 + t143 * t198 + t179) * t135 + (t113 * t197 + t125 * t196 + t181) * t114 + (t200 + t183) * (t149 ^ 2 + t151 ^ 2) * qJ(4); t129 * MDP(15); t155 * t165 + t163; t162; MDP(15); t203; t100 * MDP(21) - t101 * MDP(22) + (-t175 - t167) * t159 + t170 + t201; t204; 0; 0.2e1 * t167 + t175; t192; -t159 * MDP(27) + t201; t174; 0; MDP(27) + t167; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
