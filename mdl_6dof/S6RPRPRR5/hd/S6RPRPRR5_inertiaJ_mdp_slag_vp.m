% Calculate joint inertia matrix for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR5_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:49:36
% EndTime: 2019-03-09 03:49:37
% DurationCPUTime: 0.41s
% Computational Cost: add. (619->126), mult. (1127->172), div. (0->0), fcn. (1262->8), ass. (0->63)
t107 = sin(qJ(6));
t110 = cos(qJ(6));
t149 = MDP(31) * t107 + MDP(32) * t110;
t108 = sin(qJ(5));
t111 = cos(qJ(5));
t129 = t110 * MDP(31);
t119 = -t107 * MDP(32) + t129;
t117 = MDP(24) + t119;
t148 = -t108 * MDP(25) + t111 * t117;
t105 = sin(pkin(10));
t106 = cos(pkin(10));
t109 = sin(qJ(3));
t145 = cos(qJ(3));
t88 = t109 * t105 - t145 * t106;
t89 = t145 * t105 + t109 * t106;
t98 = -t106 * pkin(2) - pkin(1);
t78 = t88 * pkin(3) - t89 * qJ(4) + t98;
t75 = -t88 * pkin(4) - t78;
t147 = 0.2e1 * t75;
t112 = -pkin(3) - pkin(4);
t92 = t108 * qJ(4) - t111 * t112;
t90 = pkin(5) + t92;
t146 = pkin(5) + t90;
t144 = pkin(1) * MDP(7);
t143 = pkin(7) + qJ(2);
t81 = t108 * t88 + t111 * t89;
t142 = t110 * t81;
t141 = t105 * MDP(5);
t140 = t106 * MDP(4);
t80 = t108 * t89 - t111 * t88;
t139 = t80 * MDP(30);
t138 = t81 * MDP(25);
t137 = t92 * MDP(24);
t93 = t111 * qJ(4) + t108 * t112;
t136 = t93 * MDP(25);
t103 = t107 ^ 2;
t135 = t103 * MDP(26) + MDP(23);
t130 = t110 * MDP(27);
t128 = MDP(14) - MDP(17);
t126 = t107 * t130;
t127 = 0.2e1 * t126 + t135;
t94 = t143 * t105;
t95 = t143 * t106;
t82 = t109 * t95 + t145 * t94;
t125 = -pkin(3) * MDP(18) - MDP(15);
t124 = -pkin(5) * t81 - pkin(9) * t80;
t91 = -pkin(9) + t93;
t123 = -t80 * t91 + t81 * t90;
t76 = -t89 * pkin(8) + t82;
t83 = -t109 * t94 + t145 * t95;
t77 = t88 * pkin(8) + t83;
t72 = t108 * t77 - t111 * t76;
t121 = -t80 * MDP(29) + t72 * MDP(31);
t120 = MDP(28) * t110 - MDP(29) * t107;
t116 = 0.2e1 * t119;
t115 = -MDP(26) * t142 - t80 * MDP(28) - t72 * MDP(32);
t104 = t110 ^ 2;
t73 = t108 * t76 + t111 * t77;
t114 = -t80 * MDP(22) - t72 * MDP(24) - t73 * MDP(25) + (MDP(21) - (t103 - t104) * MDP(27)) * t81;
t71 = t80 * pkin(5) - t81 * pkin(9) + t75;
t70 = t107 * t71 + t110 * t73;
t69 = -t107 * t73 + t110 * t71;
t1 = [(t78 ^ 2 + t82 ^ 2 + t83 ^ 2) * MDP(18) + t138 * t147 + MDP(1) + (0.2e1 * t140 - 0.2e1 * t141 + t144) * pkin(1) + (t104 * MDP(26) + MDP(19) - 0.2e1 * t126) * t81 ^ 2 + (0.2e1 * t98 * MDP(14) - 0.2e1 * t78 * MDP(17) + MDP(8) * t89 - 0.2e1 * t88 * MDP(9)) * t89 + (MDP(24) * t147 + t139) * t80 + 0.2e1 * (t82 * t89 - t83 * t88) * MDP(16) + 0.2e1 * (t72 * t107 * t81 + t69 * t80) * MDP(31) + 0.2e1 * (t72 * t142 - t70 * t80) * MDP(32) + 0.2e1 * (t98 * MDP(13) + t78 * MDP(15)) * t88 + 0.2e1 * (-MDP(20) + t120) * t81 * t80 + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t105 ^ 2 + t106 ^ 2) * qJ(2); t78 * MDP(18) - t138 - t140 + t141 - t144 + t128 * t89 + (MDP(13) + MDP(15)) * t88 - t117 * t80; MDP(7) + MDP(18); t89 * MDP(10) - t88 * MDP(11) + (-pkin(3) * t89 - t88 * qJ(4)) * MDP(16) + (MDP(18) * qJ(4) - t128) * t83 + (-MDP(13) + t125) * t82 + (t123 * MDP(32) + t121) * t110 + (t123 * MDP(31) + t115) * t107 - t114; 0; MDP(12) + 0.2e1 * pkin(3) * MDP(15) + 0.2e1 * qJ(4) * MDP(17) + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(18) + t90 * t116 + 0.2e1 * t137 + 0.2e1 * t136 + t127; t89 * MDP(16) + t82 * MDP(18) + t149 * (-t108 * t80 - t111 * t81); 0; t125 - t148; MDP(18); (t124 * MDP(32) - t121) * t110 + (t124 * MDP(31) - t115) * t107 + t114; 0; -t137 - t136 - t146 * t129 + (t146 * MDP(32) - 0.2e1 * t130) * t107 - t135; t148; pkin(5) * t116 + t127; t69 * MDP(31) - t70 * MDP(32) + t120 * t81 + t139; -t119; (-t91 * MDP(32) - MDP(29)) * t110 + (-t91 * MDP(31) - MDP(28)) * t107; -t149 * t108; t107 * MDP(28) + t110 * MDP(29) - pkin(9) * t149; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
