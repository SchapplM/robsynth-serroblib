% Calculate joint inertia matrix for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPRP8_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:05
% EndTime: 2019-03-09 03:26:07
% DurationCPUTime: 0.60s
% Computational Cost: add. (711->141), mult. (1154->203), div. (0->0), fcn. (1197->6), ass. (0->61)
t126 = MDP(22) - MDP(25);
t104 = sin(pkin(9));
t105 = cos(pkin(9));
t108 = cos(qJ(3));
t140 = sin(qJ(3));
t89 = -t104 * t140 + t105 * t108;
t87 = t89 ^ 2;
t90 = -t104 * t108 - t105 * t140;
t88 = t90 ^ 2;
t146 = t88 + t87;
t141 = pkin(5) * t90;
t106 = sin(qJ(5));
t107 = cos(qJ(5));
t100 = t140 * pkin(3) + qJ(2);
t79 = -pkin(4) * t90 - pkin(8) * t89 + t100;
t109 = -pkin(1) - pkin(7);
t145 = -qJ(4) + t109;
t122 = t145 * t108;
t93 = t145 * t140;
t83 = t104 * t122 + t105 * t93;
t74 = -t106 * t83 + t107 * t79;
t73 = -t74 + t141;
t148 = t73 * MDP(26);
t147 = MDP(21) + MDP(23);
t102 = t106 ^ 2;
t103 = t107 ^ 2;
t130 = t102 + t103;
t120 = t130 * MDP(24);
t112 = MDP(26) * qJ(6) - t126;
t121 = -MDP(26) * pkin(5) - MDP(23);
t116 = MDP(21) - t121;
t143 = t116 * t106 - t112 * t107;
t142 = t126 * t106 - t147 * t107;
t139 = (pkin(1) * MDP(6));
t98 = pkin(3) * t104 + pkin(8);
t138 = t90 * t98;
t75 = t106 * t79 + t107 * t83;
t137 = MDP(15) * pkin(3);
t136 = qJ(6) * t90;
t135 = t107 * t89;
t134 = MDP(26) * t90;
t115 = -pkin(5) * t107 - qJ(6) * t106;
t99 = -pkin(3) * t105 - pkin(4);
t85 = t115 + t99;
t133 = t85 * MDP(26);
t132 = t90 * MDP(19);
t131 = t98 * MDP(26);
t129 = MDP(17) * t107;
t128 = t107 * MDP(18);
t125 = 0.2e1 * t106;
t123 = t140 * MDP(13);
t81 = t104 * t93 - t105 * t122;
t119 = MDP(24) + t131;
t118 = t85 * t89 + t138;
t117 = t89 * t99 + t138;
t114 = pkin(5) * t106 - qJ(6) * t107;
t72 = -t136 + t75;
t113 = t106 * t72 - t107 * t73;
t111 = t74 * MDP(21) - t75 * MDP(22);
t76 = t114 * t89 + t81;
t1 = [MDP(1) + (t100 ^ 2 + t81 ^ 2 + t83 ^ 2) * MDP(15) + t103 * t87 * MDP(16) - 0.2e1 * t89 * t90 * t128 + t88 * MDP(20) + (t72 ^ 2 + t73 ^ 2 + t76 ^ 2) * MDP(26) + (MDP(7) * t108 - 0.2e1 * t140 * MDP(8)) * t108 + (-t87 * t129 + t89 * t132) * t125 + ((-2 * MDP(4) + t139) * pkin(1)) + (0.2e1 * t140 * MDP(12) + 0.2e1 * t108 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * (t83 * MDP(14) + t73 * MDP(23) - t72 * MDP(25) - t111) * t90 + 0.2e1 * (-t113 * MDP(24) + (t106 * MDP(23) - t107 * MDP(25)) * t76 + (t106 * MDP(21) + t107 * MDP(22) + MDP(14)) * t81) * t89; MDP(4) - t139 + (-t81 * t89 - t83 * t90) * MDP(15) - t76 * t89 * MDP(26) - t72 * t134 * t107 - t90 * t148 * t106 + (-t147 * t106 - t126 * t107 - MDP(14)) * t146; MDP(6) + t146 * MDP(15) + (t130 * t88 + t87) * MDP(26); -t109 * t123 + t76 * t133 - t140 * MDP(10) + (-t102 + t103) * MDP(17) * t89 + (t109 * MDP(12) + MDP(9)) * t108 + (-t81 * MDP(21) + t117 * MDP(22) - t76 * MDP(23) - t118 * MDP(25) + t119 * t72 - t132) * t107 + (MDP(16) * t135 - t90 * MDP(18) + t117 * MDP(21) + t81 * MDP(22) + t118 * MDP(23) - t76 * MDP(25) + t119 * t73) * t106 + ((t104 * t90 - t105 * t89) * MDP(14) + (t104 * t83 - t105 * t81) * MDP(15)) * pkin(3); t108 * MDP(12) - t123 + (-t104 * t137 - t130 * t131 - t120) * t90 + (t105 * t137 - t133 - t142) * t89; MDP(11) + t102 * MDP(16) + (t130 * t98 ^ 2 + t85 ^ 2) * MDP(26) + 0.2e1 * t98 * t120 + (t104 ^ 2 + t105 ^ 2) * MDP(15) * pkin(3) ^ 2 + 0.2e1 * (-t99 * MDP(21) - t85 * MDP(23)) * t107 + (t99 * MDP(22) - t85 * MDP(25) + t129) * t125; t100 * MDP(15) + t113 * MDP(26) - t89 * t120 + t142 * t90; 0; 0; t130 * MDP(26) + MDP(15); -t90 * MDP(20) + (t74 - 0.2e1 * t141) * MDP(23) + (-0.2e1 * t136 + t75) * MDP(25) + (-pkin(5) * t73 + qJ(6) * t72) * MDP(26) + (-t106 * MDP(19) + t115 * MDP(24) + t128) * t89 + t111; t143 * t90; t106 * MDP(18) + t107 * MDP(19) - t114 * MDP(24) - t143 * t98; t112 * t106 + t116 * t107; MDP(20) + 0.2e1 * pkin(5) * MDP(23) + 0.2e1 * qJ(6) * MDP(25) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(26); t90 * MDP(23) + MDP(24) * t135 + t148; -t106 * t134; t119 * t106; -t107 * MDP(26); t121; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
