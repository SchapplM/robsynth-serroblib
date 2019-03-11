% Calculate joint inertia matrix for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPRP3_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:09:39
% EndTime: 2019-03-09 03:09:40
% DurationCPUTime: 0.57s
% Computational Cost: add. (804->163), mult. (1479->232), div. (0->0), fcn. (1458->8), ass. (0->70)
t164 = MDP(15) * qJ(4);
t118 = sin(pkin(10));
t120 = cos(pkin(10));
t146 = t118 ^ 2 + t120 ^ 2;
t163 = t146 * t164;
t122 = sin(qJ(5));
t155 = cos(qJ(5));
t139 = t155 * t120;
t103 = t118 * t122 - t139;
t104 = t118 * t155 + t122 * t120;
t143 = t120 * MDP(12);
t144 = t118 * MDP(13);
t152 = pkin(3) * MDP(15);
t128 = -t143 + t144 - t152;
t140 = MDP(22) - MDP(25);
t141 = MDP(21) + MDP(23);
t113 = -pkin(4) * t120 - pkin(3);
t89 = pkin(5) * t103 - qJ(6) * t104 + t113;
t162 = t89 * MDP(26) + t103 * t141 + t104 * t140 + t128;
t160 = -2 * MDP(17);
t159 = 2 * MDP(22);
t158 = 2 * MDP(23);
t157 = 2 * MDP(24);
t156 = 2 * MDP(25);
t124 = cos(qJ(3));
t154 = pkin(5) * t124;
t153 = pkin(8) + qJ(4);
t123 = sin(qJ(3));
t119 = sin(pkin(9));
t111 = pkin(1) * t119 + pkin(7);
t149 = t111 * t118;
t121 = cos(pkin(9));
t112 = -pkin(1) * t121 - pkin(2);
t101 = -pkin(3) * t124 - qJ(4) * t123 + t112;
t99 = t120 * t101;
t86 = -pkin(8) * t120 * t123 + t99 + (-pkin(4) - t149) * t124;
t147 = t118 * t123;
t148 = t111 * t124;
t92 = t118 * t101 + t120 * t148;
t88 = -pkin(8) * t147 + t92;
t82 = t122 * t86 + t155 * t88;
t96 = t104 * t123;
t151 = t104 * t96;
t150 = qJ(6) * t124;
t100 = pkin(4) * t147 + t123 * t111;
t145 = MDP(16) * t104;
t142 = MDP(15) + MDP(26);
t138 = t153 * t118;
t81 = -t122 * t88 + t155 * t86;
t137 = -MDP(26) * pkin(5) - MDP(23);
t136 = t146 * MDP(14);
t134 = -MDP(21) + t137;
t91 = -t118 * t148 + t99;
t133 = -t118 * t91 + t120 * t92;
t132 = MDP(26) * qJ(6) - t140;
t97 = -t122 * t147 + t123 * t139;
t131 = t97 * MDP(18) - t96 * MDP(19);
t130 = MDP(12) * t118 + MDP(13) * t120;
t129 = t104 * MDP(18) - t103 * MDP(19);
t117 = t124 ^ 2;
t116 = t123 ^ 2;
t106 = t153 * t120;
t95 = t97 ^ 2;
t94 = t106 * t155 - t122 * t138;
t93 = t106 * t122 + t138 * t155;
t90 = t97 * t103;
t83 = pkin(5) * t96 - qJ(6) * t97 + t100;
t80 = -t81 + t154;
t79 = -t150 + t82;
t1 = [MDP(1) + t116 * MDP(5) + (t111 ^ 2 * t116 + t91 ^ 2 + t92 ^ 2) * MDP(15) + t95 * MDP(16) + t97 * t96 * t160 + t117 * MDP(20) + (t79 ^ 2 + t80 ^ 2 + t83 ^ 2) * MDP(26) + (t119 ^ 2 + t121 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * (-t112 * MDP(10) + t123 * MDP(6) - t131) * t124 + 0.2e1 * (t116 * t149 - t124 * t91) * MDP(12) + 0.2e1 * (t111 * t116 * t120 + t124 * t92) * MDP(13) + 0.2e1 * (t100 * t96 - t124 * t81) * MDP(21) + (t100 * t97 + t124 * t82) * t159 + (t124 * t80 + t83 * t96) * t158 + (-t79 * t96 + t80 * t97) * t157 + (-t124 * t79 - t83 * t97) * t156 + 0.2e1 * (t112 * MDP(11) + (-t118 * t92 - t120 * t91) * MDP(14)) * t123; (-t124 * t83 + t79 * t97 + t80 * t96) * MDP(26) + (t133 - t148) * MDP(15) * t123; MDP(4) + (t116 * t146 + t117) * MDP(15) + (t96 ^ 2 + t117 + t95) * MDP(26); t97 * t145 + (-t90 - t151) * MDP(17) + (t100 * t103 + t113 * t96) * MDP(21) + (t100 * t104 + t113 * t97) * MDP(22) + (t103 * t83 + t89 * t96) * MDP(23) + (-t103 * t79 + t104 * t80 + t93 * t97 - t94 * t96) * MDP(24) + (-t104 * t83 - t89 * t97) * MDP(25) + (t79 * t94 + t80 * t93 + t83 * t89) * MDP(26) + (-t111 * MDP(11) + qJ(4) * t130 + t140 * t94 + t141 * t93 + MDP(8) - t129) * t124 + (MDP(7) - t130 * pkin(3) + (-MDP(10) + t128) * t111) * t123 + (MDP(14) + t164) * t133; (-t90 + t151) * MDP(24) + (t93 * t96 + t94 * t97) * MDP(26) + (-MDP(11) + t136 + t163) * t123 + (MDP(10) - t162) * t124; MDP(9) + (t89 ^ 2 + t93 ^ 2 + t94 ^ 2) * MDP(26) + (0.2e1 * t143 - 0.2e1 * t144 + t152) * pkin(3) + 0.2e1 * (t113 * MDP(21) + t89 * MDP(23) - t94 * MDP(24)) * t103 + (-0.2e1 * MDP(25) * t89 + t103 * t160 + t113 * t159 + t157 * t93 + t145) * t104 + (0.2e1 * t136 + t163) * qJ(4); t83 * MDP(26) + t140 * t97 + t141 * t96 + (t111 * MDP(15) + t130) * t123; -t142 * t124; t162; t142; -t124 * MDP(20) + t81 * MDP(21) - t82 * MDP(22) + (t81 - 0.2e1 * t154) * MDP(23) + (-pkin(5) * t97 - qJ(6) * t96) * MDP(24) + (-0.2e1 * t150 + t82) * MDP(25) + (-pkin(5) * t80 + qJ(6) * t79) * MDP(26) + t131; t132 * t97 + t134 * t96; (-pkin(5) * t104 - qJ(6) * t103) * MDP(24) + t132 * t94 + t134 * t93 + t129; 0; MDP(20) + pkin(5) * t158 + qJ(6) * t156 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(26); t124 * MDP(23) + t97 * MDP(24) + t80 * MDP(26); t96 * MDP(26); t104 * MDP(24) + MDP(26) * t93; 0; t137; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
