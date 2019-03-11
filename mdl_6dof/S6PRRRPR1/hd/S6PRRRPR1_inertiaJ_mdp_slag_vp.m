% Calculate joint inertia matrix for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRPR1_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:03:24
% EndTime: 2019-03-08 23:03:26
% DurationCPUTime: 0.50s
% Computational Cost: add. (789->138), mult. (1604->215), div. (0->0), fcn. (1888->12), ass. (0->76)
t137 = sin(qJ(6));
t141 = cos(qJ(6));
t172 = t137 * MDP(23) + t141 * MDP(24);
t185 = -MDP(26) * t141 + MDP(27) * t137;
t143 = cos(qJ(3));
t127 = -t143 * pkin(3) - pkin(2);
t184 = 0.2e1 * t127;
t183 = pkin(8) + pkin(9);
t138 = sin(qJ(4));
t182 = pkin(3) * t138;
t181 = MDP(20) * pkin(4);
t139 = sin(qJ(3));
t125 = t138 * t139;
t166 = t183 * t139;
t161 = t138 * t166;
t142 = cos(qJ(4));
t174 = t142 * t143;
t102 = -t125 * qJ(5) + (qJ(5) + t183) * t174 - t161;
t133 = sin(pkin(12));
t135 = cos(pkin(12));
t119 = t138 * t143 + t139 * t142;
t165 = t183 * t143;
t148 = -t138 * t165 - t142 * t166;
t147 = -t119 * qJ(5) + t148;
t90 = t102 * t133 - t135 * t147;
t180 = t141 * t90;
t162 = -t125 + t174;
t108 = t135 * t119 + t133 * t162;
t179 = t108 * t137;
t178 = t108 * t141;
t134 = sin(pkin(6));
t140 = sin(qJ(2));
t177 = t134 * t140;
t144 = cos(qJ(2));
t176 = t134 * t144;
t175 = t137 * t141;
t126 = pkin(3) * t142 + pkin(4);
t115 = t133 * t126 + t135 * t182;
t107 = t119 * t133 - t135 * t162;
t171 = MDP(25) * t107;
t111 = -t162 * pkin(4) + t127;
t168 = t111 * MDP(20);
t167 = t143 * MDP(10);
t164 = MDP(22) * t175;
t131 = t137 ^ 2;
t163 = t131 * MDP(21) + MDP(16) + 0.2e1 * t164;
t114 = t126 * t135 - t133 * t182;
t112 = -pkin(5) - t114;
t113 = pkin(10) + t115;
t160 = -t107 * t113 + t108 * t112;
t123 = pkin(4) * t133 + pkin(10);
t124 = -pkin(4) * t135 - pkin(5);
t159 = -t107 * t123 + t108 * t124;
t158 = t162 * MDP(17);
t156 = -MDP(26) * t137 - MDP(27) * t141;
t136 = cos(pkin(6));
t155 = t136 * t139 + t143 * t177;
t154 = t136 * t143 - t139 * t177;
t103 = t138 * t154 + t142 * t155;
t146 = -t138 * t155 + t142 * t154;
t93 = t103 * t133 - t135 * t146;
t153 = t146 * MDP(17) - t103 * MDP(18) + t185 * t93;
t152 = (MDP(17) * t142 - MDP(18) * t138) * pkin(3);
t151 = (MDP(23) * t141 - MDP(24) * t137) * t108;
t150 = -0.2e1 * t185;
t132 = t141 ^ 2;
t149 = t148 * MDP(17) + (-t142 * t165 + t161) * MDP(18) + t162 * MDP(15) + t119 * MDP(14) + (MDP(21) * t175 + (-t131 + t132) * MDP(22)) * t108 + t172 * t107;
t95 = t135 * t103 + t133 * t146;
t92 = t135 * t102 + t133 * t147;
t89 = t107 * pkin(5) - t108 * pkin(10) + t111;
t87 = t90 * t137;
t86 = -t137 * t176 + t141 * t95;
t85 = -t137 * t95 - t141 * t176;
t84 = t137 * t89 + t141 * t92;
t83 = -t137 * t92 + t141 * t89;
t1 = [MDP(1) + (t134 ^ 2 * t144 ^ 2 + t93 ^ 2 + t95 ^ 2) * MDP(20); (-t107 * t95 + t108 * t93) * MDP(19) + (t90 * t93 + t92 * t95) * MDP(20) + (t107 * t85 + t93 * t179) * MDP(26) + (-t107 * t86 + t93 * t178) * MDP(27) + (-t140 * MDP(4) + (-t139 * MDP(11) - t119 * MDP(18) + MDP(3) + t158 + t167 - t168) * t144) * t134; MDP(2) + 0.2e1 * pkin(2) * t167 - t158 * t184 + (t111 ^ 2 + t90 ^ 2 + t92 ^ 2) * MDP(20) + (MDP(21) * t132 - 0.2e1 * t164) * t108 ^ 2 + (-0.2e1 * MDP(11) * pkin(2) + MDP(5) * t139 + 0.2e1 * MDP(6) * t143) * t139 + (MDP(12) * t119 + 0.2e1 * t162 * MDP(13) + MDP(18) * t184) * t119 + (0.2e1 * t151 + t171) * t107 + 0.2e1 * (-t107 * t92 + t108 * t90) * MDP(19) + 0.2e1 * (t107 * t83 + t90 * t179) * MDP(26) + 0.2e1 * (-t107 * t84 + t90 * t178) * MDP(27); t154 * MDP(10) - t155 * MDP(11) + (-t114 * t93 + t115 * t95) * MDP(20) + t153; t139 * MDP(7) + t143 * MDP(8) + (-t107 * t115 - t108 * t114) * MDP(19) + (-t114 * t90 + t115 * t92) * MDP(20) + (t160 * t137 - t180) * MDP(26) + (t160 * t141 + t87) * MDP(27) + (-MDP(10) * t139 - MDP(11) * t143) * pkin(8) + t149; MDP(9) + (t114 ^ 2 + t115 ^ 2) * MDP(20) - t112 * t150 + 0.2e1 * t152 + t163; (t133 * t95 - t135 * t93) * t181 + t153; (t159 * t137 - t180) * MDP(26) + (t159 * t141 + t87) * MDP(27) + ((-t107 * t133 - t108 * t135) * MDP(19) + (t133 * t92 - t135 * t90) * MDP(20)) * pkin(4) + t149; (t114 * t135 + t115 * t133) * t181 + t152 + t163 + t185 * (t112 + t124); (t133 ^ 2 + t135 ^ 2) * MDP(20) * pkin(4) ^ 2 - t124 * t150 + t163; -MDP(20) * t176; -t107 * t185 + t168; 0; 0; MDP(20); MDP(26) * t85 - MDP(27) * t86; t83 * MDP(26) - t84 * MDP(27) + t151 + t171; t156 * t113 + t172; t156 * t123 + t172; -t185; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
