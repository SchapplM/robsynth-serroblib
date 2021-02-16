% Calculate joint inertia matrix for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR2_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:21:20
% EndTime: 2021-01-16 02:21:22
% DurationCPUTime: 0.41s
% Computational Cost: add. (463->125), mult. (939->182), div. (0->0), fcn. (1024->10), ass. (0->63)
t102 = sin(pkin(11));
t104 = cos(pkin(11));
t106 = sin(qJ(3));
t108 = cos(qJ(3));
t90 = t102 * t108 + t104 * t106;
t99 = -t108 * pkin(3) - pkin(2);
t114 = -t90 * qJ(5) + t99;
t89 = t102 * t106 - t104 * t108;
t81 = pkin(4) * t89 + t114;
t148 = -0.2e1 * t81;
t147 = 2 * MDP(12);
t142 = MDP(15) * pkin(3);
t146 = MDP(14) + MDP(16);
t129 = MDP(15) + MDP(19);
t145 = t99 * MDP(15) + t81 * MDP(19) + (MDP(12) - MDP(17)) * t89;
t144 = sin(qJ(2));
t143 = qJ(4) + pkin(8);
t105 = sin(qJ(6));
t141 = t105 * t89;
t107 = cos(qJ(6));
t140 = t107 * t89;
t139 = cos(pkin(6));
t103 = sin(pkin(6));
t109 = cos(qJ(2));
t138 = t103 * t109;
t137 = t90 * MDP(24);
t94 = pkin(3) * t102 + qJ(5);
t136 = t94 * MDP(19);
t98 = -pkin(3) * t104 - pkin(4);
t135 = t98 * MDP(19);
t134 = MDP(10) * t108;
t133 = t105 * MDP(25);
t132 = t107 * MDP(26);
t130 = -MDP(13) + MDP(18);
t125 = t143 * t106;
t92 = t143 * t108;
t82 = t102 * t92 + t104 * t125;
t84 = -t102 * t125 + t104 * t92;
t128 = t82 ^ 2 + t84 ^ 2;
t127 = t103 * t144;
t126 = t107 * t105 * MDP(21);
t124 = MDP(17) + t135;
t93 = -pkin(9) + t98;
t119 = t89 * t94 - t90 * t93;
t118 = -MDP(12) + t124;
t116 = MDP(25) * t107 - MDP(26) * t105;
t115 = -t132 - t133;
t113 = (MDP(22) * t105 + MDP(23) * t107) * t89;
t112 = t115 - t130;
t111 = -t106 * t127 + t108 * t139;
t101 = t107 ^ 2;
t100 = t105 ^ 2;
t87 = t106 * t139 + t108 * t127;
t79 = t102 * t111 + t104 * t87;
t77 = t102 * t87 - t104 * t111;
t76 = -t89 * pkin(5) + t84;
t75 = pkin(5) * t90 + t82;
t74 = (pkin(4) + pkin(9)) * t89 + t114;
t73 = -t105 * t77 + t107 * t138;
t72 = t105 * t138 + t107 * t77;
t71 = t105 * t75 + t107 * t74;
t70 = -t105 * t74 + t107 * t75;
t1 = [MDP(1) + t129 * (t103 ^ 2 * t109 ^ 2 + t77 ^ 2 + t79 ^ 2); (-t140 * t79 + t72 * t90) * MDP(25) + (t141 * t79 + t73 * t90) * MDP(26) + (-t144 * MDP(4) + (-MDP(11) * t106 + t130 * t90 + MDP(3) + t134 - t145) * t109) * t103 + t129 * (t77 * t82 + t79 * t84) + t146 * (t77 * t90 - t79 * t89); MDP(2) + 0.2e1 * pkin(2) * t134 + (t99 ^ 2 + t128) * MDP(15) + (t81 ^ 2 + t128) * MDP(19) + (-0.2e1 * MDP(11) * pkin(2) + MDP(5) * t106 + 0.2e1 * MDP(6) * t108) * t106 + (0.2e1 * t99 * MDP(13) + MDP(18) * t148 + 0.2e1 * t113 + t137) * t90 + 0.2e1 * (-t140 * t76 + t70 * t90) * MDP(25) + 0.2e1 * (t141 * t76 - t71 * t90) * MDP(26) + (t99 * t147 + MDP(17) * t148 + (MDP(20) * t100 + 0.2e1 * t126) * t89) * t89 + 0.2e1 * t146 * (t82 * t90 - t84 * t89); t111 * MDP(10) - t87 * MDP(11) + (-t104 * t142 + t118) * t77 + (t102 * t142 - t112 + t136) * t79; t98 * t90 * MDP(16) + t106 * MDP(7) + t108 * MDP(8) + (-t94 * MDP(16) + (-t100 + t101) * MDP(21)) * t89 + (t130 + t136) * t84 + t118 * t82 + (-MDP(10) * t106 - MDP(11) * t108) * pkin(8) + (t90 * MDP(22) - MDP(25) * t119 + t76 * MDP(26)) * t107 + (MDP(20) * t140 - t90 * MDP(23) + t76 * MDP(25) + MDP(26) * t119) * t105 + ((-t102 * t89 - t104 * t90) * MDP(14) + (t102 * t84 - t104 * t82) * MDP(15)) * pkin(3); -0.2e1 * t126 + t101 * MDP(20) + MDP(9) + ((2 * MDP(17)) + t135) * t98 + (0.2e1 * MDP(18) + 0.2e1 * t132 + 0.2e1 * t133 + t136) * t94 + (t104 * t147 - 0.2e1 * t102 * MDP(13) + (t102 ^ 2 + t104 ^ 2) * t142) * pkin(3); -t129 * t138; t112 * t90 + t145; 0; t129; t77 * MDP(19); t82 * MDP(19) + (MDP(16) + t116) * t90; t124; 0; MDP(19); MDP(25) * t72 + MDP(26) * t73; t70 * MDP(25) - t71 * MDP(26) + t113 + t137; (MDP(25) * t93 + MDP(22)) * t107 + (-MDP(26) * t93 - MDP(23)) * t105; t115; t116; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11), t1(16); t1(2), t1(3), t1(5), t1(8), t1(12), t1(17); t1(4), t1(5), t1(6), t1(9), t1(13), t1(18); t1(7), t1(8), t1(9), t1(10), t1(14), t1(19); t1(11), t1(12), t1(13), t1(14), t1(15), t1(20); t1(16), t1(17), t1(18), t1(19), t1(20), t1(21);];
Mq = res;
