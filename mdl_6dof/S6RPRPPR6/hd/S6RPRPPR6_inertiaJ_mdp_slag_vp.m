% Calculate joint inertia matrix for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPPR6_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:54:39
% EndTime: 2019-03-09 02:54:41
% DurationCPUTime: 0.42s
% Computational Cost: add. (718->129), mult. (1216->187), div. (0->0), fcn. (1331->8), ass. (0->66)
t123 = sin(pkin(9));
t112 = pkin(3) * t123 + qJ(5);
t122 = sin(pkin(10));
t124 = cos(pkin(10));
t153 = t122 ^ 2 + t124 ^ 2;
t142 = t153 * MDP(19);
t168 = t142 * t112;
t129 = -pkin(1) - pkin(7);
t167 = -qJ(4) + t129;
t126 = sin(qJ(6));
t127 = cos(qJ(6));
t125 = cos(pkin(9));
t128 = cos(qJ(3));
t160 = sin(qJ(3));
t104 = -t123 * t128 - t125 * t160;
t103 = -t123 * t160 + t125 * t128;
t154 = t103 * t124;
t118 = t160 * pkin(3) + qJ(2);
t94 = -pkin(4) * t104 - qJ(5) * t103 + t118;
t109 = t167 * t160;
t144 = t167 * t128;
t98 = t125 * t109 + t123 * t144;
t83 = -t122 * t98 + t124 * t94;
t81 = -pkin(5) * t104 - pkin(8) * t154 + t83;
t155 = t103 * t122;
t84 = t122 * t94 + t124 * t98;
t82 = -pkin(8) * t155 + t84;
t108 = t122 * t127 + t124 * t126;
t86 = t108 * t103;
t166 = MDP(25) * (-t126 * t82 + t127 * t81) - MDP(26) * (t126 * t81 + t127 * t82) - MDP(23) * t86;
t143 = t153 * MDP(18);
t96 = t109 * t123 - t125 * t144;
t165 = t96 ^ 2;
t164 = t103 ^ 2;
t117 = -pkin(3) * t125 - pkin(4);
t110 = -pkin(5) * t124 + t117;
t162 = 0.2e1 * t110;
t161 = -2 * MDP(21);
t159 = (pkin(1) * MDP(6));
t158 = pkin(8) + t112;
t157 = MDP(15) * pkin(3);
t152 = MDP(16) * t122;
t151 = MDP(16) * t124;
t150 = MDP(17) * t122;
t149 = MDP(17) * t124;
t148 = MDP(19) * t117;
t147 = MDP(20) * t108;
t105 = t122 * t126 - t127 * t124;
t146 = MDP(25) * t105;
t145 = t160 * MDP(13);
t141 = t122 * t84 + t124 * t83;
t140 = t122 * t83 - t124 * t84;
t88 = t105 * t103;
t138 = t86 * MDP(25) - t88 * MDP(26);
t137 = (MDP(25) * t108 - MDP(26) * t105) * t104;
t136 = t149 + t152;
t135 = -t108 * MDP(26) - t146;
t134 = MDP(14) + t136;
t100 = t158 * t124;
t99 = t158 * t122;
t133 = t108 * MDP(22) - t105 * MDP(23) + (-t100 * t126 - t127 * t99) * MDP(25) - (t100 * t127 - t126 * t99) * MDP(26);
t132 = -t135 + t150 - t151;
t131 = t132 + t148;
t102 = t104 ^ 2;
t85 = pkin(5) * t155 + t96;
t1 = [MDP(1) + (t118 ^ 2 + t98 ^ 2 + t165) * MDP(15) + (t83 ^ 2 + t84 ^ 2 + t165) * MDP(19) + t102 * MDP(24) + (MDP(7) * t128 - 0.2e1 * t160 * MDP(8)) * t128 + ((-2 * MDP(4) + t159) * pkin(1)) - (-MDP(20) * t88 - 0.2e1 * MDP(22) * t104 + t86 * t161) * t88 + (0.2e1 * t160 * MDP(12) + 0.2e1 * t128 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * t138 * t85 + 0.2e1 * (-t141 * MDP(18) + t134 * t96) * t103 + 0.2e1 * (MDP(14) * t98 - MDP(16) * t83 + MDP(17) * t84 - t166) * t104; -t102 * t152 - t159 + MDP(4) + (-t98 * MDP(15) + t140 * MDP(19) + (-MDP(14) - t149) * t104 - t137) * t104 + (-t134 * t103 + (-MDP(15) - MDP(19)) * t96 - t138) * t103; MDP(6) + (t102 + t164) * MDP(15) + (t153 * t102 + t164) * MDP(19); -t160 * MDP(10) - t129 * t145 + (t117 * t155 - t124 * t96) * MDP(16) + (t117 * t154 + t122 * t96) * MDP(17) - t140 * MDP(18) + (-t140 * t112 + t117 * t96) * MDP(19) - t88 * t147 + (t105 * t88 - t108 * t86) * MDP(21) + (t105 * t85 + t110 * t86) * MDP(25) + (t108 * t85 - t110 * t88) * MDP(26) + (MDP(12) * t129 + MDP(9)) * t128 + (t136 * t112 - t133) * t104 + ((-t103 * t125 + t104 * t123) * MDP(14) + (t123 * t98 - t125 * t96) * MDP(15)) * pkin(3); t128 * MDP(12) - t145 + (-t123 * t157 - t143 - t168) * t104 + (t125 * t157 - t131) * t103; t146 * t162 + MDP(11) + (t123 ^ 2 + t125 ^ 2) * MDP(15) * pkin(3) ^ 2 + (t148 + 0.2e1 * t150 - 0.2e1 * t151) * t117 + (MDP(26) * t162 + t105 * t161 + t147) * t108 + (0.2e1 * t143 + t168) * t112; t118 * MDP(15) + t141 * MDP(19) - t103 * t143 + t132 * t104; 0; 0; MDP(15) + t142; t96 * MDP(19) + t136 * t103 + t138; -t103 * MDP(19); t131; 0; MDP(19); -MDP(22) * t88 - MDP(24) * t104 + t166; t137; t133; t135; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
