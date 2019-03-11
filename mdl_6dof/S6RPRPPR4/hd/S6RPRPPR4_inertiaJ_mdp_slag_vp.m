% Calculate joint inertia matrix for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPPR4_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:25
% EndTime: 2019-03-09 02:48:26
% DurationCPUTime: 0.53s
% Computational Cost: add. (788->141), mult. (1469->188), div. (0->0), fcn. (1591->8), ass. (0->76)
t151 = MDP(18) + MDP(22);
t131 = sin(pkin(10));
t133 = cos(pkin(10));
t163 = t131 ^ 2 + t133 ^ 2;
t187 = t151 * t163;
t186 = (MDP(17) + MDP(20)) * t163;
t135 = sin(qJ(6));
t137 = cos(qJ(6));
t185 = -t131 * t137 + t133 * t135;
t153 = MDP(16) - MDP(21);
t154 = MDP(15) + MDP(19);
t184 = t154 * t131 + t153 * t133;
t111 = t131 * t135 + t133 * t137;
t107 = t111 * MDP(28);
t165 = MDP(29) * t185 - t107;
t183 = -t153 * t131 + t154 * t133 - t165;
t150 = qJ(5) * t131 + pkin(3);
t179 = pkin(4) + pkin(5);
t106 = t179 * t133 + t150;
t182 = 0.2e1 * t106;
t134 = cos(pkin(9));
t124 = -pkin(2) * t134 - pkin(1);
t181 = 0.2e1 * t124;
t180 = -2 * MDP(24);
t178 = cos(qJ(3));
t177 = pkin(1) * MDP(7);
t132 = sin(pkin(9));
t136 = sin(qJ(3));
t115 = t178 * t132 + t136 * t134;
t176 = pkin(8) * t115;
t175 = pkin(7) + qJ(2);
t174 = -pkin(8) + qJ(4);
t173 = pkin(3) * MDP(18);
t119 = t175 * t132;
t121 = t175 * t134;
t104 = -t136 * t119 + t178 * t121;
t112 = t132 * t136 - t178 * t134;
t97 = pkin(3) * t112 - qJ(4) * t115 + t124;
t92 = t133 * t104 + t131 * t97;
t172 = qJ(5) * t133;
t170 = t132 * MDP(5);
t168 = t134 * MDP(4);
t94 = t185 * t115;
t167 = t94 * MDP(26);
t95 = t111 * t115;
t166 = t95 * MDP(25);
t161 = MDP(27) * t112;
t160 = t185 * MDP(23);
t116 = -pkin(4) * t133 - t150;
t159 = t116 * MDP(22);
t158 = t131 * MDP(16);
t157 = t131 * MDP(21);
t156 = t133 * MDP(15);
t155 = t133 * MDP(19);
t88 = t112 * qJ(5) + t92;
t99 = t131 * t104;
t91 = t133 * t97 - t99;
t102 = t178 * t119 + t121 * t136;
t89 = -pkin(4) * t112 - t91;
t149 = t131 * t89 + t133 * t88;
t148 = t131 * t88 - t133 * t89;
t147 = t131 * t92 + t133 * t91;
t146 = -t131 * t91 + t133 * t92;
t86 = t99 + (-t97 - t176) * t133 - t179 * t112;
t87 = t131 * t176 + t88;
t145 = (-t135 * t87 + t137 * t86) * MDP(28) - (t135 * t86 + t137 * t87) * MDP(29);
t144 = t94 * MDP(28) + t95 * MDP(29);
t143 = MDP(15) * t131 + MDP(16) * t133;
t142 = MDP(19) * t131 - MDP(21) * t133;
t141 = MDP(28) * t137 - t135 * MDP(29);
t118 = t174 * t131;
t120 = t174 * t133;
t140 = -MDP(25) * t185 - MDP(26) * t111 + MDP(28) * (t118 * t137 - t120 * t135) - MDP(29) * (t118 * t135 + t120 * t137);
t93 = (pkin(4) * t131 - t172) * t115 + t102;
t90 = (-t179 * t131 + t172) * t115 - t102;
t1 = [(t102 ^ 2 + t91 ^ 2 + t92 ^ 2) * MDP(18) + (t88 ^ 2 + t89 ^ 2 + t93 ^ 2) * MDP(22) + MDP(1) + (MDP(23) * t95 + t180 * t94) * t95 + (MDP(14) * t181 + MDP(8) * t115) * t115 + (0.2e1 * t168 - 0.2e1 * t170 + t177) * pkin(1) + (MDP(13) * t181 - 0.2e1 * MDP(9) * t115 + t161 - 0.2e1 * t166 + 0.2e1 * t167) * t112 + 0.2e1 * t144 * t90 + 0.2e1 * (MDP(15) * t91 - MDP(16) * t92 - MDP(19) * t89 + MDP(21) * t88 - t145) * t112 + 0.2e1 * (-MDP(17) * t147 - MDP(20) * t148 + t102 * t143 + t142 * t93) * t115 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t132 ^ 2 + t134 ^ 2) * qJ(2); -t168 + t170 - t177 + t147 * MDP(18) + t148 * MDP(22) + (MDP(13) + t183) * t112 + (MDP(14) - t186) * t115; MDP(7) + t187; -t104 * MDP(14) + t146 * MDP(17) + t149 * MDP(20) - t95 * t160 + (-t111 * t95 + t185 * t94) * MDP(24) + (t106 * t94 + t111 * t90) * MDP(28) + (t106 * t95 - t185 * t90) * MDP(29) + (-t155 - t157 + t159) * t93 + (-pkin(3) * t143 + t116 * t142 + MDP(10)) * t115 + (-MDP(13) - t156 + t158 - t173) * t102 + (t146 * MDP(18) + t149 * MDP(22)) * qJ(4) + (-qJ(4) * t184 - MDP(11) - t140) * t112; 0; MDP(12) + t107 * t182 + (-0.2e1 * t155 - 0.2e1 * t157 + t159) * t116 + (0.2e1 * t156 - 0.2e1 * t158 + t173) * pkin(3) - (MDP(29) * t182 + t111 * t180 - t160) * t185 + (qJ(4) * t187 + 0.2e1 * t186) * qJ(4); t102 * MDP(18) + t93 * MDP(22) + t115 * t184 - t144; 0; t159 - t173 - t183; t151; t115 * t133 * MDP(20) + t89 * MDP(22) + (-MDP(19) - t141) * t112; -t133 * MDP(22); (MDP(22) * qJ(4) + MDP(20)) * t131; 0; MDP(22); t145 - t161 + t166 - t167; t165; t140; 0; t141; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
