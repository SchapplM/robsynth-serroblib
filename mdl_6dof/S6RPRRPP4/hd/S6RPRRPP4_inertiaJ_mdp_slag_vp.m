% Calculate joint inertia matrix for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRPP4_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:40:48
% EndTime: 2019-03-09 04:40:50
% DurationCPUTime: 0.66s
% Computational Cost: add. (1391->171), mult. (2556->246), div. (0->0), fcn. (2869->8), ass. (0->69)
t142 = sin(pkin(10));
t144 = cos(pkin(10));
t146 = sin(qJ(4));
t148 = cos(qJ(4));
t189 = -t142 * t146 + t144 * t148;
t188 = MDP(23) + MDP(27);
t145 = cos(pkin(9));
t136 = -pkin(2) * t145 - pkin(1);
t187 = 0.2e1 * t136;
t186 = 2 * MDP(22);
t185 = 2 * MDP(24);
t184 = 2 * MDP(25);
t183 = cos(qJ(3));
t182 = pkin(1) * MDP(7);
t181 = pkin(7) + qJ(2);
t180 = -qJ(5) - pkin(8);
t143 = sin(pkin(9));
t147 = sin(qJ(3));
t124 = t143 * t147 - t145 * t183;
t126 = t143 * t183 + t147 * t145;
t177 = t126 * t148;
t108 = pkin(3) * t124 - pkin(8) * t126 + t136;
t128 = t181 * t143;
t129 = t181 * t145;
t110 = -t147 * t128 + t129 * t183;
t97 = t148 * t108 - t110 * t146;
t94 = pkin(4) * t124 - qJ(5) * t177 + t97;
t179 = t110 * t148;
t96 = t179 + (-qJ(5) * t126 + t108) * t146;
t90 = t142 * t94 + t144 * t96;
t178 = t126 * t146;
t175 = t143 * MDP(5);
t173 = t145 * MDP(4);
t172 = t146 * t148;
t125 = t142 * t148 + t144 * t146;
t137 = -pkin(4) * t148 - pkin(3);
t107 = -pkin(5) * t189 - qJ(6) * t125 + t137;
t170 = MDP(27) * t107;
t169 = t189 * MDP(24);
t168 = t125 * MDP(26);
t167 = t126 * MDP(14);
t131 = pkin(4) * t142 + qJ(6);
t166 = t131 * MDP(26);
t130 = t180 * t148;
t162 = t180 * t146;
t113 = -t130 * t142 - t144 * t162;
t115 = -t144 * t130 + t142 * t162;
t165 = t113 ^ 2 + t115 ^ 2;
t87 = t124 * qJ(6) + t90;
t163 = MDP(16) * t172;
t89 = -t142 * t96 + t144 * t94;
t104 = t125 * t126;
t105 = t189 * t126;
t161 = -t115 * t104 + t105 * t113;
t109 = t128 * t183 + t147 * t129;
t156 = t148 * MDP(17) - t146 * MDP(18);
t155 = t148 * MDP(20) - t146 * MDP(21);
t154 = -MDP(20) * t146 - MDP(21) * t148;
t153 = -t168 - t169;
t101 = pkin(4) * t178 + t109;
t152 = -t153 + t155;
t151 = t146 * MDP(17) + t148 * MDP(18) - t113 * MDP(24) + t115 * MDP(26) + pkin(8) * t154;
t141 = t148 ^ 2;
t140 = t146 ^ 2;
t134 = pkin(4) * t144 + pkin(5);
t98 = t108 * t146 + t179;
t91 = t104 * pkin(5) - t105 * qJ(6) + t101;
t88 = -pkin(5) * t124 - t89;
t1 = [MDP(1) + t167 * t187 + (t101 ^ 2 + t89 ^ 2 + t90 ^ 2) * MDP(23) + (t87 ^ 2 + t88 ^ 2 + t91 ^ 2) * MDP(27) + (0.2e1 * t173 - 0.2e1 * t175 + t182) * pkin(1) + (t141 * MDP(15) + MDP(8) - 0.2e1 * t163) * t126 ^ 2 + (MDP(13) * t187 + t124 * MDP(19) + 0.2e1 * (-MDP(9) + t156) * t126) * t124 + 0.2e1 * (t109 * t178 + t124 * t97) * MDP(20) + 0.2e1 * (t109 * t177 - t124 * t98) * MDP(21) + (-t104 * t90 - t105 * t89) * t186 + (t104 * t91 - t124 * t88) * t185 + (-t104 * t87 + t105 * t88) * t184 + 0.2e1 * (-t105 * t91 + t124 * t87) * MDP(26) + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t143 ^ 2 + t145 ^ 2) * qJ(2); -t173 + t175 - t182 + t167 + (t125 * t90 + t189 * t89) * MDP(23) + (t125 * t87 - t189 * t88) * MDP(27) + (MDP(13) + t152) * t124 + (MDP(22) + MDP(25)) * (-t125 * t104 - t105 * t189); MDP(7) + t188 * (t125 ^ 2 + t189 ^ 2); -t110 * MDP(14) + (-t125 * t89 + t189 * t90 + t161) * MDP(22) + (t101 * t137 - t113 * t89 + t115 * t90) * MDP(23) + (t104 * t107 - t189 * t91) * MDP(24) + (t125 * t88 + t189 * t87 + t161) * MDP(25) + (-t105 * t107 - t125 * t91) * MDP(26) + (t107 * t91 + t113 * t88 + t115 * t87) * MDP(27) + (-MDP(13) - t155) * t109 + (-MDP(11) + t151) * t124 + (MDP(10) + MDP(15) * t172 + (-t140 + t141) * MDP(16) + t154 * pkin(3)) * t126; t188 * (-t113 * t189 + t125 * t115); MDP(12) + t140 * MDP(15) + 0.2e1 * t163 + (t137 ^ 2 + t165) * MDP(23) + t165 * MDP(27) + (-0.2e1 * t168 - 0.2e1 * t169 + t170) * t107 + 0.2e1 * t155 * pkin(3) + (t186 + t184) * (t113 * t125 + t115 * t189); t97 * MDP(20) - t98 * MDP(21) + t89 * MDP(24) + (-t104 * t131 - t105 * t134) * MDP(25) + t87 * MDP(26) + (t131 * t87 - t134 * t88) * MDP(27) + t156 * t126 + (MDP(19) + (pkin(5) + t134) * MDP(24) + t166) * t124 + ((-t104 * t142 - t105 * t144) * MDP(22) + (t142 * t90 + t144 * t89) * MDP(23)) * pkin(4); (t125 * t131 + t134 * t189) * MDP(27) + (t125 * t142 + t144 * t189) * MDP(23) * pkin(4) + t152; (-t125 * t134 + t131 * t189) * MDP(25) + (-t113 * t134 + t115 * t131) * MDP(27) + ((-t125 * t144 + t142 * t189) * MDP(22) + (-t113 * t144 + t115 * t142) * MDP(23)) * pkin(4) + t151; MDP(19) + (t131 ^ 2 + t134 ^ 2) * MDP(27) + (t142 ^ 2 + t144 ^ 2) * MDP(23) * pkin(4) ^ 2 + t134 * t185 + 0.2e1 * t166; MDP(23) * t101 + t104 * MDP(24) - t105 * MDP(26) + t91 * MDP(27); 0; MDP(23) * t137 + t153 + t170; 0; t188; -t124 * MDP(24) + t105 * MDP(25) + t88 * MDP(27); -t189 * MDP(27); t125 * MDP(25) + t113 * MDP(27); -MDP(27) * t134 - MDP(24); 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
