% Calculate joint inertia matrix for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRRP1_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:59:19
% EndTime: 2019-03-08 23:59:21
% DurationCPUTime: 0.50s
% Computational Cost: add. (681->156), mult. (1353->220), div. (0->0), fcn. (1483->10), ass. (0->75)
t133 = sin(qJ(5));
t137 = cos(qJ(5));
t159 = t133 * MDP(21) + t137 * MDP(22);
t153 = t133 * MDP(25);
t173 = t153 - MDP(17);
t134 = sin(qJ(4));
t135 = sin(qJ(3));
t138 = cos(qJ(4));
t139 = cos(qJ(3));
t113 = t134 * t139 + t138 * t135;
t172 = 0.2e1 * t113;
t124 = -t139 * pkin(3) - pkin(2);
t171 = 0.2e1 * t124;
t170 = 2 * MDP(26);
t169 = -pkin(9) - pkin(8);
t168 = t138 * pkin(3);
t122 = -pkin(4) - t168;
t167 = pkin(4) - t122;
t166 = MDP(27) * pkin(5);
t117 = t169 * t135;
t118 = t169 * t139;
t99 = t134 * t117 - t138 * t118;
t165 = t137 * t99;
t131 = sin(pkin(6));
t136 = sin(qJ(2));
t164 = t131 * t136;
t140 = cos(qJ(2));
t163 = t131 * t140;
t162 = t133 * t137;
t128 = t137 * qJ(6);
t132 = cos(pkin(6));
t104 = t132 * t139 - t135 * t164;
t105 = t132 * t135 + t139 * t164;
t91 = t134 * t104 + t138 * t105;
t85 = -t133 * t163 + t137 * t91;
t161 = t85 * MDP(25);
t98 = -t138 * t117 - t134 * t118;
t160 = t98 * MDP(24);
t158 = MDP(10) * t139;
t157 = MDP(18) * t113;
t156 = MDP(22) * t133;
t155 = MDP(26) * t113;
t112 = t134 * t135 - t138 * t139;
t154 = t112 * MDP(23);
t152 = t133 * MDP(26);
t151 = t137 * MDP(24);
t123 = -t137 * pkin(5) - pkin(4);
t150 = MDP(20) * t162;
t129 = t133 ^ 2;
t149 = t129 * MDP(19) + MDP(16) + 0.2e1 * t150;
t94 = t112 * pkin(4) - t113 * pkin(10) + t124;
t82 = -t133 * t99 + t137 * t94;
t148 = -pkin(4) * t113 - pkin(10) * t112;
t147 = t82 * MDP(24) - (t133 * t94 + t165) * MDP(25);
t121 = t134 * pkin(3) + pkin(10);
t146 = -t112 * t121 + t113 * t122;
t145 = t151 - t153;
t144 = t133 * MDP(24) + t137 * MDP(25);
t143 = (t138 * MDP(17) - t134 * MDP(18)) * pkin(3);
t84 = -t133 * t91 - t137 * t163;
t90 = -t138 * t104 + t134 * t105;
t142 = (-t84 * t133 + t85 * t137) * MDP(26) - t91 * MDP(18) + (-t151 + t173) * t90;
t130 = t137 ^ 2;
t81 = t165 + (-qJ(6) * t113 + t94) * t133;
t141 = t81 * t137 * MDP(26) - t99 * MDP(18) + t173 * t98 + (MDP(19) * t162 + MDP(14) + (-t129 + t130) * MDP(20)) * t113 + (-MDP(15) + t159) * t112;
t116 = t137 * pkin(10) + t128;
t115 = (-qJ(6) - pkin(10)) * t133;
t114 = t123 - t168;
t111 = t116 * t137;
t109 = t137 * t121 + t128;
t108 = (-qJ(6) - t121) * t133;
t101 = t109 * t137;
t87 = t133 * t113 * pkin(5) + t98;
t79 = t112 * pkin(5) - t113 * t128 + t82;
t1 = [MDP(1) + (t84 ^ 2 + t85 ^ 2 + t90 ^ 2) * MDP(27); (t84 * t79 + t85 * t81 + t90 * t87) * MDP(27) + (t84 * MDP(24) - t161) * t112 + ((-t133 * t85 - t137 * t84) * MDP(26) + t144 * t90) * t113 + (-t136 * MDP(4) + (-MDP(11) * t135 - MDP(17) * t112 + MDP(3) - t157 + t158) * t140) * t131; MDP(2) + 0.2e1 * pkin(2) * t158 + t157 * t171 + (t79 ^ 2 + t81 ^ 2 + t87 ^ 2) * MDP(27) + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t135 + 0.2e1 * t139 * MDP(6)) * t135 + (t130 * MDP(19) + MDP(12) - 0.2e1 * t150) * t113 ^ 2 + ((-t133 * t81 - t137 * t79) * MDP(26) + t144 * t98) * t172 + (MDP(17) * t171 + t154 + (MDP(21) * t137 - MDP(13) - t156) * t172 + 0.2e1 * t147) * t112; t104 * MDP(10) - t105 * MDP(11) + (t84 * t108 + t85 * t109 + t90 * t114) * MDP(27) + t142; (t79 * t108 + t81 * t109 + t87 * t114) * MDP(27) + t135 * MDP(7) + t139 * MDP(8) + (t146 * MDP(25) - t108 * t155 - t160) * t137 + (t146 * MDP(24) + (-t109 * t113 - t79) * MDP(26)) * t133 + (-t135 * MDP(10) - t139 * MDP(11)) * pkin(8) + t141; MDP(9) + (-t108 * t133 + t101) * t170 + (t108 ^ 2 + t109 ^ 2 + t114 ^ 2) * MDP(27) + t149 - 0.2e1 * t145 * t122 + 0.2e1 * t143; (t84 * t115 + t85 * t116 + t90 * t123) * MDP(27) + t142; (t79 * t115 + t81 * t116 + t87 * t123) * MDP(27) + (t148 * MDP(25) - t115 * t155 - t160) * t137 + (t148 * MDP(24) + (-t113 * t116 - t79) * MDP(26)) * t133 + t141; (t101 + t111) * MDP(26) + (t108 * t115 + t109 * t116 + t114 * t123) * MDP(27) + t167 * t151 + t143 + (-t167 * MDP(25) + (-t108 - t115) * MDP(26)) * t133 + t149; (-t115 * t133 + t111) * t170 + (t115 ^ 2 + t116 ^ 2 + t123 ^ 2) * MDP(27) + 0.2e1 * t145 * pkin(4) + t149; -t161 + (MDP(24) + t166) * t84; t79 * t166 + t154 + (-t156 + (-MDP(26) * pkin(5) + MDP(21)) * t137) * t113 + t147; -t144 * t121 + (t108 * MDP(27) - t152) * pkin(5) + t159; -t144 * pkin(10) + (MDP(27) * t115 - t152) * pkin(5) + t159; MDP(27) * pkin(5) ^ 2 + MDP(23); t90 * MDP(27); t87 * MDP(27); t114 * MDP(27); t123 * MDP(27); 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
