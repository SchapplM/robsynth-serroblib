% Calculate joint inertia matrix for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRRP2_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:01:17
% EndTime: 2019-03-09 06:01:19
% DurationCPUTime: 0.52s
% Computational Cost: add. (675->150), mult. (1267->218), div. (0->0), fcn. (1254->8), ass. (0->72)
t123 = sin(qJ(4));
t126 = cos(qJ(4));
t131 = t123 * MDP(17) + t126 * MDP(18);
t122 = sin(qJ(5));
t125 = cos(qJ(5));
t145 = t125 * t126;
t103 = t122 * t123 - t145;
t104 = t122 * t126 + t125 * t123;
t157 = pkin(8) + pkin(9);
t107 = t157 * t123;
t108 = t157 * t126;
t87 = -t125 * t107 - t122 * t108;
t88 = -t122 * t107 + t125 * t108;
t135 = t104 * MDP(21) - t103 * MDP(22) + t87 * MDP(24) - t88 * MDP(25);
t164 = t123 * MDP(14) + t126 * MDP(15) - t131 * pkin(8) + t135;
t163 = -2 * MDP(20);
t139 = t125 * MDP(24);
t162 = (-t122 * MDP(25) + t139) * pkin(4);
t124 = sin(qJ(3));
t95 = t104 * t124;
t147 = t123 * t124;
t97 = -t122 * t147 + t124 * t145;
t153 = t97 * MDP(21) - t95 * MDP(22);
t161 = t104 * MDP(19);
t127 = cos(qJ(3));
t120 = sin(pkin(10));
t111 = t120 * pkin(1) + pkin(7);
t150 = t111 * t123;
t155 = pkin(9) * t124;
t121 = cos(pkin(10));
t112 = -t121 * pkin(1) - pkin(2);
t102 = -t127 * pkin(3) - t124 * pkin(8) + t112;
t98 = t126 * t102;
t77 = -t126 * t155 + t98 + (-pkin(4) - t150) * t127;
t148 = t111 * t127;
t137 = t126 * t148;
t80 = t137 + (t102 - t155) * t123;
t74 = -t122 * t80 + t125 * t77;
t75 = t122 * t77 + t125 * t80;
t160 = t74 * MDP(24) - t75 * MDP(25) + t153;
t158 = 2 * MDP(26);
t156 = pkin(4) * t122;
t154 = -t95 * MDP(24) - t97 * MDP(25);
t152 = MDP(27) * pkin(5);
t151 = t104 * t95;
t149 = t111 * t126;
t146 = t123 * t126;
t115 = -t126 * pkin(4) - pkin(3);
t89 = t103 * pkin(5) + t115;
t142 = t89 * MDP(27);
t99 = pkin(4) * t147 + t124 * t111;
t140 = t124 * MDP(11);
t138 = MDP(16) + MDP(23);
t136 = MDP(13) * t146;
t134 = t126 * MDP(14) - t123 * MDP(15);
t132 = t126 * MDP(17) - t123 * MDP(18);
t130 = t103 * MDP(24) + t104 * MDP(25);
t119 = t127 ^ 2;
t118 = t126 ^ 2;
t117 = t124 ^ 2;
t116 = t123 ^ 2;
t114 = t125 * pkin(4) + pkin(5);
t94 = t97 ^ 2;
t84 = t123 * t102 + t137;
t83 = -t123 * t148 + t98;
t82 = t97 * t103;
t81 = t95 * pkin(5) + t99;
t79 = -t103 * qJ(6) + t88;
t78 = -t104 * qJ(6) + t87;
t73 = -t95 * qJ(6) + t75;
t72 = -t127 * pkin(5) - t97 * qJ(6) + t74;
t1 = [MDP(1) + 0.2e1 * t112 * t140 + t94 * MDP(19) + t97 * t95 * t163 + (t72 ^ 2 + t73 ^ 2 + t81 ^ 2) * MDP(27) + (t120 ^ 2 + t121 ^ 2) * MDP(4) * pkin(1) ^ 2 + t138 * t119 + (t118 * MDP(12) + MDP(5) - 0.2e1 * t136) * t117 + 0.2e1 * (-t112 * MDP(10) + (MDP(6) - t134) * t124 - t153) * t127 + 0.2e1 * (t117 * t150 - t83 * t127) * MDP(17) + 0.2e1 * (t117 * t149 + t84 * t127) * MDP(18) + 0.2e1 * (-t74 * t127 + t99 * t95) * MDP(24) + 0.2e1 * (t75 * t127 + t99 * t97) * MDP(25) + (-t72 * t97 - t73 * t95) * t158; (-t81 * t127 - t72 * t95 + t73 * t97) * MDP(27); MDP(4) + (t95 ^ 2 + t119 + t94) * MDP(27); t97 * t161 + (-t82 - t151) * MDP(20) + (t99 * t103 + t115 * t95) * MDP(24) + (t99 * t104 + t115 * t97) * MDP(25) + (-t73 * t103 - t72 * t104 - t78 * t97 - t79 * t95) * MDP(26) + (t72 * t78 + t73 * t79 + t81 * t89) * MDP(27) + (-t111 * MDP(11) + MDP(8) - t164) * t127 + (MDP(7) - t111 * MDP(10) + MDP(12) * t146 + (-t116 + t118) * MDP(13) + (-pkin(3) * t123 - t149) * MDP(17) + (-pkin(3) * t126 + t150) * MDP(18)) * t124; -t140 + (-t82 + t151) * MDP(26) + (-t95 * t78 + t97 * t79) * MDP(27) + (MDP(10) - t130 + t132 - t142) * t127; MDP(9) + t116 * MDP(12) + 0.2e1 * t136 - t79 * t103 * t158 + (t78 ^ 2 + t79 ^ 2 + t89 ^ 2) * MDP(27) + (t103 * t163 - t78 * t158 + t161) * t104 + 0.2e1 * t132 * pkin(3) + 0.2e1 * t130 * t115; t83 * MDP(17) - t84 * MDP(18) + (-t114 * t97 - t95 * t156) * MDP(26) + (t72 * t114 + t73 * t156) * MDP(27) + (-t138 - t162) * t127 + t134 * t124 + t160; (-t95 * t114 + t97 * t156) * MDP(27) - t131 * t124 + t154; (-t103 * t156 - t114 * t104) * MDP(26) + (t78 * t114 + t79 * t156) * MDP(27) + t164; t114 ^ 2 * MDP(27) + (0.2e1 * t139 + (MDP(27) * t156 - 0.2e1 * MDP(25)) * t122) * pkin(4) + t138; -t127 * MDP(23) + (-t97 * MDP(26) + t72 * MDP(27)) * pkin(5) + t160; -t95 * t152 + t154; (-t104 * MDP(26) + t78 * MDP(27)) * pkin(5) + t135; t114 * t152 + MDP(23) + t162; MDP(27) * pkin(5) ^ 2 + MDP(23); t81 * MDP(27); -t127 * MDP(27); t142; 0; 0; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
