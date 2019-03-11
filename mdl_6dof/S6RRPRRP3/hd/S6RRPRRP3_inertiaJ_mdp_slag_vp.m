% Calculate joint inertia matrix for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRRP3_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:50:36
% EndTime: 2019-03-09 11:50:38
% DurationCPUTime: 0.67s
% Computational Cost: add. (1124->158), mult. (2082->244), div. (0->0), fcn. (2335->8), ass. (0->74)
t139 = sin(pkin(10));
t132 = pkin(2) * t139 + pkin(8);
t142 = sin(qJ(4));
t144 = cos(qJ(4));
t148 = t142 * MDP(18) + t144 * MDP(19);
t170 = pkin(9) + t132;
t117 = t170 * t142;
t118 = t170 * t144;
t141 = sin(qJ(5));
t173 = cos(qJ(5));
t102 = -t173 * t117 - t118 * t141;
t103 = -t141 * t117 + t118 * t173;
t127 = t141 * t144 + t142 * t173;
t178 = -t141 * t142 + t173 * t144;
t150 = t127 * MDP(22) + MDP(23) * t178 + t102 * MDP(25) - t103 * MDP(26);
t182 = t142 * MDP(15) + t144 * MDP(16) - t132 * t148 + t150;
t149 = t144 * MDP(18) - t142 * MDP(19);
t161 = MDP(25) * t178 - t127 * MDP(26);
t181 = t161 + t149;
t180 = -qJ(3) - pkin(7);
t140 = cos(pkin(10));
t143 = sin(qJ(2));
t174 = cos(qJ(2));
t122 = t139 * t174 + t140 * t143;
t97 = t127 * t122;
t98 = t178 * t122;
t179 = t98 * MDP(22) - t97 * MDP(23);
t177 = (MDP(15) * t144 - MDP(16) * t142) * t122;
t176 = -2 * MDP(21);
t175 = 2 * MDP(27);
t121 = t139 * t143 - t140 * t174;
t172 = pkin(4) * t121;
t171 = pkin(4) * t141;
t169 = MDP(28) * pkin(5);
t168 = t178 * t98;
t129 = t180 * t143;
t130 = t180 * t174;
t107 = t129 * t139 - t130 * t140;
t167 = t107 * t144;
t166 = t122 * t142;
t165 = t122 * t144;
t163 = t142 * t144;
t162 = t98 * MDP(20);
t158 = 0.2e1 * t174;
t157 = MDP(17) + MDP(24);
t156 = t121 * MDP(24) + t179;
t155 = t173 * pkin(4);
t136 = -pkin(2) * t174 - pkin(1);
t104 = t121 * pkin(3) - t122 * pkin(8) + t136;
t87 = t167 + (-pkin(9) * t122 + t104) * t142;
t154 = t173 * t87;
t133 = -pkin(2) * t140 - pkin(3);
t152 = MDP(14) * t163;
t151 = MDP(25) * t173;
t89 = t144 * t104 - t107 * t142;
t86 = -pkin(9) * t165 + t172 + t89;
t83 = -t141 * t87 + t173 * t86;
t105 = -t140 * t129 - t130 * t139;
t94 = pkin(4) * t166 + t105;
t128 = -pkin(4) * t144 + t133;
t84 = t141 * t86 + t154;
t138 = t144 ^ 2;
t137 = t142 ^ 2;
t135 = t155 + pkin(5);
t123 = t127 ^ 2;
t108 = -pkin(5) * t178 + t128;
t93 = qJ(6) * t178 + t103;
t92 = -qJ(6) * t127 + t102;
t91 = t127 * t97;
t90 = t104 * t142 + t167;
t88 = pkin(5) * t97 + t94;
t82 = -t97 * qJ(6) + t84;
t81 = pkin(5) * t121 - qJ(6) * t98 + t83;
t1 = [MDP(1) + pkin(1) * MDP(9) * t158 + (t105 ^ 2 + t107 ^ 2 + t136 ^ 2) * MDP(12) + (t81 ^ 2 + t82 ^ 2 + t88 ^ 2) * MDP(28) + (t176 * t97 + t162) * t98 + (t138 * MDP(13) - 0.2e1 * t152) * t122 ^ 2 + t157 * t121 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t143 + MDP(5) * t158) * t143 + 0.2e1 * (t177 + t179) * t121 + 0.2e1 * (t105 * t122 - t107 * t121) * MDP(11) + 0.2e1 * (t105 * t166 + t121 * t89) * MDP(18) + 0.2e1 * (t105 * t165 - t121 * t90) * MDP(19) + 0.2e1 * (t121 * t83 + t94 * t97) * MDP(25) + 0.2e1 * (-t121 * t84 + t94 * t98) * MDP(26) + (-t81 * t98 - t82 * t97) * t175; t143 * MDP(6) + t174 * MDP(7) + t127 * t162 + (-t91 + t168) * MDP(21) + (t128 * t97 - t178 * t94) * MDP(25) + (t127 * t94 + t128 * t98) * MDP(26) + (-t127 * t81 + t178 * t82 - t92 * t98 - t93 * t97) * MDP(27) + (t108 * t88 + t81 * t92 + t82 * t93) * MDP(28) - t149 * t105 + (-MDP(10) * t174 - t143 * MDP(9)) * pkin(7) + (MDP(13) * t163 + (-t137 + t138) * MDP(14) + t148 * t133) * t122 + t182 * t121 + ((-t121 * t139 - t122 * t140) * MDP(11) + (-t105 * t140 + t107 * t139) * MDP(12)) * pkin(2); MDP(8) + t137 * MDP(13) + 0.2e1 * t152 + t123 * MDP(20) - t127 * t178 * t176 + (-t127 * t92 + t178 * t93) * t175 + (t108 ^ 2 + t92 ^ 2 + t93 ^ 2) * MDP(28) + (t139 ^ 2 + t140 ^ 2) * MDP(12) * pkin(2) ^ 2 - 0.2e1 * t128 * t161 - 0.2e1 * t133 * t149; t136 * MDP(12) + (-t91 - t168) * MDP(27) + (t127 * t82 + t178 * t81) * MDP(28) + t181 * t121; (t127 * t93 + t178 * t92) * MDP(28); MDP(12) + (t178 ^ 2 + t123) * MDP(28); t121 * MDP(17) + t89 * MDP(18) - t90 * MDP(19) + (t121 * t155 + t83) * MDP(25) + (-t154 + (-t86 - t172) * t141) * MDP(26) + (-t135 * t98 - t171 * t97) * MDP(27) + (t135 * t81 + t171 * t82) * MDP(28) + t156 + t177; (-t127 * t135 + t171 * t178) * MDP(27) + (t135 * t92 + t171 * t93) * MDP(28) + t182; (t127 * t171 + t135 * t178) * MDP(28) + t181; t135 ^ 2 * MDP(28) + (0.2e1 * t151 + (MDP(28) * t171 - 0.2e1 * MDP(26)) * t141) * pkin(4) + t157; t83 * MDP(25) - t84 * MDP(26) + (-t98 * MDP(27) + t81 * MDP(28)) * pkin(5) + t156; (-t127 * MDP(27) + t92 * MDP(28)) * pkin(5) + t150; t169 * t178 + t161; t135 * t169 + MDP(24) + (-MDP(26) * t141 + t151) * pkin(4); MDP(28) * pkin(5) ^ 2 + MDP(24); t88 * MDP(28); t108 * MDP(28); 0; 0; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
