% Calculate joint inertia matrix for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPP5_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:43
% EndTime: 2019-03-09 04:44:45
% DurationCPUTime: 0.52s
% Computational Cost: add. (906->174), mult. (1611->216), div. (0->0), fcn. (1682->6), ass. (0->67)
t123 = sin(pkin(9));
t124 = cos(pkin(9));
t126 = sin(qJ(3));
t165 = cos(qJ(3));
t106 = t165 * t123 + t126 * t124;
t169 = 0.2e1 * t106;
t105 = t126 * t123 - t165 * t124;
t125 = sin(qJ(4));
t127 = cos(qJ(4));
t114 = -t124 * pkin(2) - pkin(1);
t92 = t105 * pkin(3) - t106 * pkin(8) + t114;
t163 = pkin(7) + qJ(2);
t108 = t163 * t123;
t109 = t163 * t124;
t97 = -t126 * t108 + t165 * t109;
t87 = -t125 * t97 + t127 * t92;
t85 = -t105 * pkin(4) - t87;
t88 = t125 * t92 + t127 * t97;
t168 = t87 * MDP(20) - t88 * MDP(21) - t85 * MDP(22);
t128 = pkin(4) + pkin(5);
t156 = t127 * qJ(5);
t167 = t128 * t125 - t156;
t166 = 0.2e1 * t114;
t164 = pkin(1) * MDP(7);
t162 = pkin(8) - qJ(6);
t161 = pkin(4) * MDP(22);
t160 = pkin(4) * MDP(25);
t159 = t106 * t127;
t158 = t123 * MDP(5);
t157 = t124 * MDP(4);
t154 = t127 * MDP(26) + t125 * MDP(27);
t117 = t125 * qJ(5);
t153 = t127 * pkin(4) + t117;
t121 = t125 ^ 2;
t122 = t127 ^ 2;
t151 = t121 + t122;
t150 = t106 * MDP(14);
t149 = t127 * MDP(16);
t148 = MDP(20) + MDP(22);
t147 = -MDP(21) + MDP(24);
t146 = MDP(22) + MDP(26);
t145 = MDP(23) - MDP(28);
t144 = MDP(24) + MDP(27);
t143 = MDP(25) + MDP(29);
t100 = t105 * qJ(5);
t142 = 0.2e1 * t100 + t88;
t107 = -pkin(3) - t153;
t84 = t100 + t88;
t141 = MDP(25) * pkin(8) + MDP(23);
t96 = t165 * t108 + t126 * t109;
t140 = t128 * MDP(29) + MDP(22);
t139 = pkin(4) * t125 - t156;
t98 = t125 * t106 * qJ(6);
t83 = t98 + t84;
t86 = -t167 * t106 - t96;
t89 = t139 * t106 + t96;
t137 = -t96 * MDP(20) - t89 * MDP(22) + t86 * MDP(26) - t83 * MDP(28);
t82 = -t105 * pkin(5) - qJ(6) * t159 + t85;
t136 = t96 * MDP(21) - t89 * MDP(24) + t86 * MDP(27) - t82 * MDP(28);
t102 = t127 * pkin(5) - t107;
t111 = t162 * t127;
t135 = pkin(3) * MDP(20) - t107 * MDP(22) + t102 * MDP(26) - t111 * MDP(28);
t110 = t162 * t125;
t134 = -pkin(3) * MDP(21) - t107 * MDP(24) + t102 * MDP(27) - t110 * MDP(28);
t133 = t125 * MDP(17) + t127 * MDP(18) - t110 * MDP(26) + t111 * MDP(27);
t130 = qJ(5) ^ 2;
t1 = [t150 * t166 + (t84 ^ 2 + t85 ^ 2 + t89 ^ 2) * MDP(25) + (t82 ^ 2 + t83 ^ 2 + t86 ^ 2) * MDP(29) + MDP(1) + (0.2e1 * t157 - 0.2e1 * t158 + t164) * pkin(1) + (t122 * MDP(15) - 0.2e1 * t125 * t149 + MDP(8)) * t106 ^ 2 + (MDP(13) * t166 + t105 * MDP(19) + (t127 * MDP(17) - t125 * MDP(18) - MDP(9)) * t169) * t105 + 0.2e1 * (t84 * MDP(24) - t82 * MDP(26) + t83 * MDP(27) + t168) * t105 + ((t85 * MDP(23) + t136) * t127 + (-t84 * MDP(23) - t137) * t125) * t169 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t123 ^ 2 + t124 ^ 2) * qJ(2); -t157 + t158 - t164 + t150 + (t84 * t125 - t85 * t127) * MDP(25) + (t83 * t125 - t82 * t127) * MDP(29) - t145 * t151 * t106 + (MDP(13) + (MDP(20) + t146) * t127 + (-MDP(21) + t144) * t125) * t105; t143 * t151 + MDP(7); -t96 * MDP(13) - t97 * MDP(14) + t89 * t107 * MDP(25) + (t86 * t102 + t82 * t110 + t83 * t111) * MDP(29) + (t141 * t84 + t137) * t127 + (t141 * t85 + t136) * t125 + (-MDP(11) + (-t148 * t125 + t147 * t127) * pkin(8) + t133) * t105 + (MDP(10) + (-t121 + t122) * MDP(16) + t134 * t127 + (t127 * MDP(15) - t135) * t125) * t106; (-t127 * t110 + t125 * t111) * MDP(29); MDP(12) + t121 * MDP(15) + (t151 * pkin(8) ^ 2 + t107 ^ 2) * MDP(25) + (t102 ^ 2 + t110 ^ 2 + t111 ^ 2) * MDP(29) + 0.2e1 * t151 * MDP(23) * pkin(8) + 0.2e1 * t135 * t127 + 0.2e1 * (t134 + t149) * t125; t142 * MDP(24) + (-t85 * pkin(4) + t84 * qJ(5)) * MDP(25) - t85 * MDP(26) + (t98 + t142) * MDP(27) + (t83 * qJ(5) - t82 * t128) * MDP(29) + (MDP(19) + t161 + (pkin(5) + t128) * MDP(26)) * t105 + ((-t145 * qJ(5) - MDP(18)) * t125 + (-pkin(4) * MDP(23) + qJ(6) * MDP(26) + t128 * MDP(28) + MDP(17)) * t127) * t106 + t168; t153 * MDP(25) + t117 * MDP(29) + (MDP(20) + t140) * t127 + t147 * t125 + t154; -t139 * MDP(23) + t167 * MDP(28) + (t111 * qJ(5) - t110 * t128) * MDP(29) + ((qJ(5) * MDP(25) + t147) * t127 + (-t148 - t160) * t125) * pkin(8) + t133; MDP(19) + 0.2e1 * t161 + (pkin(4) ^ 2 + t130) * MDP(25) + 0.2e1 * t128 * MDP(26) + (t128 ^ 2 + t130) * MDP(29) + 0.2e1 * t144 * qJ(5); t85 * MDP(25) + t82 * MDP(29) - t146 * t105 + t145 * t159; -t143 * t127; t110 * MDP(29) + (-MDP(28) + t141) * t125; -MDP(26) - t140 - t160; t143; t86 * MDP(29) + (-t125 * MDP(26) + t127 * MDP(27)) * t106; 0; t102 * MDP(29) + t154; 0; 0; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
