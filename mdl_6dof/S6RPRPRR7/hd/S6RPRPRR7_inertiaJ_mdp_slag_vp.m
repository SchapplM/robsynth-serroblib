% Calculate joint inertia matrix for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR7_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:05
% EndTime: 2019-03-09 03:56:06
% DurationCPUTime: 0.39s
% Computational Cost: add. (669->110), mult. (1140->156), div. (0->0), fcn. (1309->8), ass. (0->60)
t133 = sin(qJ(6));
t136 = cos(qJ(6));
t144 = t136 * MDP(28) - t133 * MDP(29);
t173 = MDP(21) + t144;
t131 = sin(pkin(10));
t132 = cos(pkin(10));
t135 = sin(qJ(3));
t138 = cos(qJ(3));
t114 = -t131 * t135 + t132 * t138;
t115 = -t131 * t138 - t132 * t135;
t134 = sin(qJ(5));
t137 = cos(qJ(5));
t100 = t137 * t114 + t134 * t115;
t172 = t100 ^ 2;
t162 = t133 * MDP(25) + t136 * MDP(26);
t171 = MDP(28) * t133 + MDP(29) * t136;
t124 = t135 * pkin(3) + qJ(2);
t104 = -t115 * pkin(4) + t124;
t170 = 0.2e1 * t104;
t169 = (pkin(1) * MDP(6));
t168 = pkin(3) * t131;
t139 = -pkin(1) - pkin(7);
t163 = -qJ(4) + t139;
t117 = t163 * t135;
t118 = t163 * t138;
t102 = -t131 * t117 + t132 * t118;
t87 = -t114 * pkin(8) + t102;
t103 = t132 * t117 + t131 * t118;
t88 = t115 * pkin(8) + t103;
t83 = t134 * t88 - t137 * t87;
t80 = t83 * t133;
t166 = t83 * t136;
t165 = t133 * t136;
t146 = t134 * t114 - t137 * t115;
t164 = t146 * MDP(27);
t159 = t100 * MDP(22);
t123 = t132 * pkin(3) + pkin(4);
t107 = t137 * t123 - t134 * t168;
t158 = t107 * MDP(21);
t108 = -t134 * t123 - t137 * t168;
t157 = t108 * MDP(22);
t154 = MDP(24) * t165;
t129 = t133 ^ 2;
t153 = t129 * MDP(23) + MDP(20) + 0.2e1 * t154;
t152 = t114 ^ 2 + t115 ^ 2;
t151 = -pkin(5) * t100 - pkin(9) * t146;
t105 = -pkin(5) - t107;
t106 = pkin(9) - t108;
t149 = t100 * t105 - t106 * t146;
t148 = t102 * t114 - t103 * t115;
t147 = t114 * t132 - t115 * t131;
t145 = MDP(25) * t136 - MDP(26) * t133;
t142 = -MDP(22) * t146 + t173 * t100;
t130 = t136 ^ 2;
t84 = t134 * t87 + t137 * t88;
t141 = -t83 * MDP(21) - t84 * MDP(22) + (-MDP(19) + t162) * t146 + ((-t129 + t130) * MDP(24) + MDP(23) * t165 + MDP(18)) * t100;
t85 = pkin(5) * t146 - t100 * pkin(9) + t104;
t79 = t133 * t85 + t136 * t84;
t78 = -t133 * t84 + t136 * t85;
t1 = [MDP(1) + (t102 ^ 2 + t103 ^ 2 + t124 ^ 2) * MDP(15) + t159 * t170 + (MDP(7) * t138 - 0.2e1 * t135 * MDP(8)) * t138 + ((-2 * MDP(4) + t169) * pkin(1)) + (t130 * MDP(23) + MDP(16) - 0.2e1 * t154) * t172 + (0.2e1 * t135 * MDP(12) + 0.2e1 * t138 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + (MDP(21) * t170 + t164 + 0.2e1 * (-MDP(17) + t145) * t100) * t146 - 0.2e1 * t148 * MDP(14) + 0.2e1 * (t100 * t80 + t146 * t78) * MDP(28) + 0.2e1 * (t100 * t166 - t146 * t79) * MDP(29); -t152 * MDP(14) + t148 * MDP(15) + MDP(4) - t169 + t171 * (-t146 ^ 2 - t172); t152 * MDP(15) + MDP(6); (t149 * t133 - t166) * MDP(28) + (t149 * t136 + t80) * MDP(29) + (t139 * MDP(12) + MDP(9)) * t138 + (-t139 * MDP(13) - MDP(10)) * t135 + (-t147 * MDP(14) + (t102 * t132 + t103 * t131) * MDP(15)) * pkin(3) + t141; t147 * MDP(15) * pkin(3) + t138 * MDP(12) - t135 * MDP(13) + t142; MDP(11) + (t131 ^ 2 + t132 ^ 2) * MDP(15) * pkin(3) ^ 2 - 0.2e1 * t144 * t105 + 0.2e1 * t158 + 0.2e1 * t157 + t153; t124 * MDP(15) + t173 * t146 + t159; 0; 0; MDP(15); (t151 * t133 - t166) * MDP(28) + (t151 * t136 + t80) * MDP(29) + t141; t142; t153 + t157 + t158 + t144 * (pkin(5) - t105); 0; 0.2e1 * pkin(5) * t144 + t153; t78 * MDP(28) - t79 * MDP(29) + t145 * t100 + t164; -t171 * t146; -t106 * t171 + t162; t144; -pkin(9) * t171 + t162; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
