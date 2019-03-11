% Calculate joint inertia matrix for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRPR2_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:02:35
% EndTime: 2019-03-09 05:02:37
% DurationCPUTime: 0.57s
% Computational Cost: add. (858->154), mult. (1653->233), div. (0->0), fcn. (1735->10), ass. (0->77)
t150 = sin(qJ(4));
t153 = cos(qJ(4));
t159 = t150 * MDP(17) + t153 * MDP(18);
t182 = -qJ(5) - pkin(8);
t132 = t182 * t150;
t133 = t182 * t153;
t145 = sin(pkin(11));
t147 = cos(pkin(11));
t115 = t147 * t132 + t145 * t133;
t128 = t145 * t153 + t147 * t150;
t104 = -t128 * pkin(9) + t115;
t116 = t145 * t132 - t147 * t133;
t127 = -t145 * t150 + t147 * t153;
t105 = t127 * pkin(9) + t116;
t149 = sin(qJ(6));
t152 = cos(qJ(6));
t111 = -t152 * t127 + t149 * t128;
t112 = t149 * t127 + t152 * t128;
t163 = t112 * MDP(23) - t111 * MDP(24) + (t152 * t104 - t149 * t105) * MDP(26) - (t149 * t104 + t152 * t105) * MDP(27);
t189 = t150 * MDP(14) + t153 * MDP(15) - t159 * pkin(8) + t163;
t151 = sin(qJ(3));
t119 = t128 * t151;
t174 = t151 * t153;
t176 = t150 * t151;
t120 = -t145 * t176 + t147 * t174;
t101 = t152 * t119 + t149 * t120;
t102 = -t149 * t119 + t152 * t120;
t180 = t102 * MDP(23) - t101 * MDP(24);
t154 = cos(qJ(3));
t148 = cos(pkin(10));
t139 = -t148 * pkin(1) - pkin(2);
t126 = -t154 * pkin(3) - t151 * pkin(8) + t139;
t121 = t153 * t126;
t146 = sin(pkin(10));
t137 = t146 * pkin(1) + pkin(7);
t179 = t137 * t150;
t103 = -qJ(5) * t174 + t121 + (-pkin(4) - t179) * t154;
t177 = t137 * t154;
t165 = t153 * t177;
t107 = t165 + (-qJ(5) * t151 + t126) * t150;
t93 = t147 * t103 - t145 * t107;
t87 = -t154 * pkin(5) - t120 * pkin(9) + t93;
t94 = t145 * t103 + t147 * t107;
t88 = -t119 * pkin(9) + t94;
t84 = -t149 * t88 + t152 * t87;
t85 = t149 * t87 + t152 * t88;
t188 = t84 * MDP(26) - t85 * MDP(27) + t180;
t186 = 2 * MDP(19);
t185 = -2 * MDP(22);
t184 = 0.2e1 * MDP(27);
t183 = pkin(4) * t145;
t181 = -t101 * MDP(26) - t102 * MDP(27);
t178 = t137 * t153;
t175 = t150 * t153;
t124 = pkin(4) * t176 + t151 * t137;
t171 = t102 * MDP(21);
t170 = t111 * MDP(26);
t138 = t147 * pkin(4) + pkin(5);
t169 = (t152 * t138 - t149 * t183) * MDP(26);
t168 = (t149 * t138 + t152 * t183) * MDP(27);
t167 = t151 * MDP(11);
t166 = MDP(16) + MDP(25);
t140 = -t153 * pkin(4) - pkin(3);
t164 = MDP(13) * t175;
t162 = t153 * MDP(14) - t150 * MDP(15);
t160 = t153 * MDP(17) - t150 * MDP(18);
t158 = MDP(25) - t168 + t169;
t157 = t140 * MDP(20) + t112 * MDP(27) + t170;
t144 = t154 ^ 2;
t143 = t153 ^ 2;
t142 = t151 ^ 2;
t141 = t150 ^ 2;
t118 = -t127 * pkin(5) + t140;
t114 = t150 * t126 + t165;
t113 = -t150 * t177 + t121;
t108 = t119 * pkin(5) + t124;
t1 = [MDP(1) + 0.2e1 * t139 * t167 + (t124 ^ 2 + t93 ^ 2 + t94 ^ 2) * MDP(20) + (t146 ^ 2 + t148 ^ 2) * MDP(4) * pkin(1) ^ 2 + t166 * t144 + (t101 * t185 + t171) * t102 + (t143 * MDP(12) + MDP(5) - 0.2e1 * t164) * t142 + 0.2e1 * (-t139 * MDP(10) + (MDP(6) - t162) * t151 - t180) * t154 + 0.2e1 * (-t113 * t154 + t142 * t179) * MDP(17) + 0.2e1 * (t114 * t154 + t142 * t178) * MDP(18) + (-t94 * t119 - t93 * t120) * t186 + 0.2e1 * (t108 * t101 - t84 * t154) * MDP(26) + (t108 * t102 + t85 * t154) * t184; (-t93 * t119 + t94 * t120 - t124 * t154) * MDP(20); MDP(4) + (t119 ^ 2 + t120 ^ 2 + t144) * MDP(20); (-t115 * t120 - t116 * t119 + t94 * t127 - t93 * t128) * MDP(19) + (t93 * t115 + t94 * t116 + t124 * t140) * MDP(20) + t112 * t171 + (-t112 * t101 - t102 * t111) * MDP(22) + (t118 * t101 + t108 * t111) * MDP(26) + (t118 * t102 + t108 * t112) * MDP(27) + (-t137 * MDP(11) + MDP(8) - t189) * t154 + (MDP(7) - t137 * MDP(10) + MDP(12) * t175 + (-t141 + t143) * MDP(13) + (-pkin(3) * t150 - t178) * MDP(17) + (-pkin(3) * t153 + t179) * MDP(18)) * t151; -t167 + (t119 * t128 + t120 * t127) * MDP(19) + (-t119 * t115 + t120 * t116) * MDP(20) + (MDP(10) - t157 + t160) * t154; MDP(9) + t141 * MDP(12) + 0.2e1 * t164 + (-t115 * t128 + t116 * t127) * t186 + (t115 ^ 2 + t116 ^ 2 + t140 ^ 2) * MDP(20) + 0.2e1 * t118 * t170 + 0.2e1 * t160 * pkin(3) + (MDP(21) * t112 + t111 * t185 + t118 * t184) * t112; t113 * MDP(17) - t114 * MDP(18) + (-MDP(16) - t158) * t154 + t162 * t151 + ((-t119 * t145 - t120 * t147) * MDP(19) + (t145 * t94 + t147 * t93) * MDP(20)) * pkin(4) + t188; -t159 * t151 + (-t119 * t147 + t120 * t145) * MDP(20) * pkin(4) + t181; ((t127 * t145 - t128 * t147) * MDP(19) + (t115 * t147 + t116 * t145) * MDP(20)) * pkin(4) + t189; (t145 ^ 2 + t147 ^ 2) * MDP(20) * pkin(4) ^ 2 + 0.2e1 * t169 - 0.2e1 * t168 + t166; t124 * MDP(20) - t181; -t154 * MDP(20); t157; 0; MDP(20); -t154 * MDP(25) + t188; t181; t163; t158; 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
