% Calculate joint inertia matrix for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RRPRPP1_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:47:53
% EndTime: 2019-03-09 09:47:55
% DurationCPUTime: 0.68s
% Computational Cost: add. (1474->180), mult. (2649->272), div. (0->0), fcn. (2956->8), ass. (0->66)
t180 = MDP(21) + MDP(25);
t147 = cos(qJ(4));
t142 = sin(pkin(9));
t163 = pkin(2) * t142 + pkin(8);
t156 = qJ(5) + t163;
t120 = t156 * t147;
t141 = sin(pkin(10));
t143 = cos(pkin(10));
t145 = sin(qJ(4));
t153 = t156 * t145;
t108 = t120 * t141 + t143 * t153;
t110 = t143 * t120 - t141 * t153;
t179 = -t108 * MDP(22) + t110 * MDP(24);
t178 = 2 * MDP(20);
t177 = 0.2e1 * MDP(22);
t176 = 2 * MDP(23);
t175 = cos(qJ(2));
t144 = cos(pkin(9));
t146 = sin(qJ(2));
t125 = t142 * t146 - t144 * t175;
t127 = t142 * t175 + t144 * t146;
t172 = t147 * t127;
t138 = -t175 * pkin(2) - pkin(1);
t111 = t125 * pkin(3) - t127 * pkin(8) + t138;
t130 = (-qJ(3) - pkin(7)) * t146;
t164 = t175 * pkin(7);
t131 = t175 * qJ(3) + t164;
t115 = t130 * t142 + t131 * t144;
t97 = t147 * t111 - t115 * t145;
t93 = pkin(4) * t125 - qJ(5) * t172 + t97;
t174 = t115 * t147;
t95 = t174 + (-qJ(5) * t127 + t111) * t145;
t89 = t141 * t93 + t143 * t95;
t173 = t145 * t127;
t123 = t141 * t145 - t143 * t147;
t126 = t141 * t147 + t143 * t145;
t137 = -pkin(2) * t144 - pkin(3);
t129 = -pkin(4) * t147 + t137;
t103 = pkin(5) * t123 - qJ(6) * t126 + t129;
t171 = MDP(25) * t103;
t170 = t123 * MDP(22);
t169 = t126 * MDP(24);
t132 = pkin(4) * t141 + qJ(6);
t168 = t132 * MDP(24);
t167 = 0.2e1 * t175;
t166 = t108 ^ 2 + t110 ^ 2;
t86 = t125 * qJ(6) + t89;
t162 = t145 * t147 * MDP(14);
t88 = -t141 * t95 + t143 * t93;
t104 = t126 * t127;
t105 = -t141 * t173 + t143 * t172;
t161 = -t110 * t104 + t105 * t108;
t113 = -t144 * t130 + t131 * t142;
t102 = pkin(4) * t173 + t113;
t155 = -t147 * MDP(18) + t145 * MDP(19);
t154 = t169 - t170;
t152 = (t147 * MDP(15) - t145 * MDP(16)) * t127;
t151 = -t125 * t163 + t137 * t127;
t150 = t154 - t155;
t140 = t147 ^ 2;
t139 = t145 ^ 2;
t135 = pkin(4) * t143 + pkin(5);
t98 = t111 * t145 + t174;
t90 = pkin(5) * t104 - qJ(6) * t105 + t102;
t87 = -pkin(5) * t125 - t88;
t1 = [MDP(1) + pkin(1) * MDP(9) * t167 + (t113 ^ 2 + t115 ^ 2 + t138 ^ 2) * MDP(12) + (t102 ^ 2 + t88 ^ 2 + t89 ^ 2) * MDP(21) + (t86 ^ 2 + t87 ^ 2 + t90 ^ 2) * MDP(25) + (t140 * MDP(13) - 0.2e1 * t162) * t127 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t146 + MDP(5) * t167) * t146 + (t125 * MDP(17) + 0.2e1 * t152) * t125 + 0.2e1 * (t113 * t127 - t115 * t125) * MDP(11) + 0.2e1 * (t113 * t173 + t125 * t97) * MDP(18) + 0.2e1 * (t113 * t172 - t125 * t98) * MDP(19) + (-t104 * t89 - t105 * t88) * t178 + (t104 * t90 - t125 * t87) * t177 + (-t104 * t86 + t105 * t87) * t176 + 0.2e1 * (-t105 * t90 + t125 * t86) * MDP(24); t175 * MDP(7) - MDP(10) * t164 + t145 * MDP(13) * t172 + (-t139 + t140) * t127 * MDP(14) + (-t113 * t147 + t145 * t151) * MDP(18) + (t113 * t145 + t147 * t151) * MDP(19) + (-t123 * t89 - t126 * t88 + t161) * MDP(20) + (t102 * t129 - t108 * t88 + t110 * t89) * MDP(21) + (t103 * t104 + t123 * t90) * MDP(22) + (-t123 * t86 + t126 * t87 + t161) * MDP(23) + (-t103 * t105 - t126 * t90) * MDP(24) + (t103 * t90 + t108 * t87 + t110 * t86) * MDP(25) + (-pkin(7) * MDP(9) + MDP(6)) * t146 + (t145 * MDP(15) + t147 * MDP(16) + t179) * t125 + ((-t125 * t142 - t127 * t144) * MDP(11) + (-t113 * t144 + t115 * t142) * MDP(12)) * pkin(2); MDP(8) + t139 * MDP(13) + 0.2e1 * t162 + (t129 ^ 2 + t166) * MDP(21) + t166 * MDP(25) + (t142 ^ 2 + t144 ^ 2) * MDP(12) * pkin(2) ^ 2 + 0.2e1 * t155 * t137 + (-0.2e1 * t169 + 0.2e1 * t170 + t171) * t103 + (t178 + t176) * (t108 * t126 - t110 * t123); t138 * MDP(12) + (-t123 * t88 + t126 * t89) * MDP(21) + (t123 * t87 + t126 * t86) * MDP(25) + t150 * t125 + (MDP(20) + MDP(23)) * (-t126 * t104 + t105 * t123); t180 * (t108 * t123 + t110 * t126); MDP(12) + t180 * (t123 ^ 2 + t126 ^ 2); t97 * MDP(18) - t98 * MDP(19) + t88 * MDP(22) + (-t104 * t132 - t105 * t135) * MDP(23) + t86 * MDP(24) + (t132 * t86 - t135 * t87) * MDP(25) + t152 + (MDP(17) + (pkin(5) + t135) * MDP(22) + t168) * t125 + ((-t104 * t141 - t105 * t143) * MDP(20) + (t141 * t89 + t143 * t88) * MDP(21)) * pkin(4); (-t123 * t132 - t126 * t135) * MDP(23) + (-t108 * t135 + t110 * t132) * MDP(25) + (-MDP(19) * t163 + MDP(16)) * t147 + (-MDP(18) * t163 + MDP(15)) * t145 + ((-t123 * t141 - t126 * t143) * MDP(20) + (-t108 * t143 + t110 * t141) * MDP(21)) * pkin(4) + t179; (-t123 * t135 + t126 * t132) * MDP(25) + (-t123 * t143 + t126 * t141) * MDP(21) * pkin(4) + t150; MDP(17) + (t132 ^ 2 + t135 ^ 2) * MDP(25) + (t141 ^ 2 + t143 ^ 2) * MDP(21) * pkin(4) ^ 2 + t135 * t177 + 0.2e1 * t168; MDP(21) * t102 + t104 * MDP(22) - t105 * MDP(24) + t90 * MDP(25); MDP(21) * t129 - t154 + t171; 0; 0; t180; -t125 * MDP(22) + t105 * MDP(23) + t87 * MDP(25); t126 * MDP(23) + t108 * MDP(25); t123 * MDP(25); -t135 * MDP(25) - MDP(22); 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
