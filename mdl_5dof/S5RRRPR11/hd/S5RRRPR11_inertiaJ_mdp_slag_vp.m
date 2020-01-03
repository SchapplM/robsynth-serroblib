% Calculate joint inertia matrix for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR11_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:15
% EndTime: 2019-12-31 21:35:18
% DurationCPUTime: 0.59s
% Computational Cost: add. (391->137), mult. (741->191), div. (0->0), fcn. (644->6), ass. (0->61)
t106 = sin(qJ(5));
t109 = cos(qJ(5));
t111 = cos(qJ(2));
t108 = sin(qJ(2));
t110 = cos(qJ(3));
t140 = t108 * t110;
t101 = t111 * pkin(3);
t107 = sin(qJ(3));
t141 = t107 * t111;
t94 = -pkin(2) * t111 - pkin(7) * t108 - pkin(1);
t80 = -pkin(6) * t141 + t110 * t94;
t79 = t101 - t80;
t73 = pkin(4) * t111 - pkin(8) * t140 + t79;
t81 = t110 * t111 * pkin(6) + t107 * t94;
t78 = -qJ(4) * t111 + t81;
t74 = pkin(8) * t107 * t108 + t78;
t142 = t107 * t109;
t83 = t106 * t140 - t108 * t142;
t88 = t106 * t107 + t109 * t110;
t84 = t88 * t108;
t154 = t84 * MDP(24) - t83 * MDP(25) - (t106 * t74 - t109 * t73) * MDP(27) - (t106 * t73 + t109 * t74) * MDP(28);
t153 = 2 * MDP(19);
t152 = t80 * MDP(16) - t81 * MDP(17) - t154;
t151 = 0.2e1 * t111;
t144 = qJ(4) * t107;
t147 = pkin(3) + pkin(4);
t85 = t147 * t110 + pkin(2) + t144;
t149 = 0.2e1 * t85;
t148 = -2 * MDP(23);
t146 = pkin(7) - pkin(8);
t145 = pkin(3) * MDP(21);
t143 = qJ(4) * t110;
t137 = t84 * MDP(22);
t135 = (qJ(4) * t106 + t109 * t147) * MDP(27);
t134 = (t109 * qJ(4) - t106 * t147) * MDP(28);
t102 = t107 ^ 2;
t104 = t110 ^ 2;
t133 = t102 + t104;
t132 = MDP(12) * t110;
t131 = qJ(4) * MDP(20);
t130 = MDP(15) + MDP(26);
t129 = MDP(17) - MDP(20);
t128 = t146 * t107;
t127 = -pkin(3) * t110 - t144;
t126 = -pkin(3) * t107 + t143;
t124 = t107 * t79 + t110 * t78;
t93 = -pkin(2) + t127;
t123 = pkin(2) * MDP(16) - t93 * MDP(18);
t122 = -MDP(17) * pkin(2) - MDP(20) * t93;
t118 = t110 * MDP(13) - t107 * MDP(14);
t117 = t109 * MDP(27) - t106 * MDP(28);
t116 = -MDP(26) - t134 - t135;
t115 = MDP(18) + t117;
t89 = -t106 * t110 + t142;
t95 = t146 * t110;
t114 = t89 * MDP(24) - t88 * MDP(25) - (t106 * t95 - t109 * t128) * MDP(27) - (t106 * t128 + t109 * t95) * MDP(28);
t113 = -t107 * MDP(13) + t114;
t98 = pkin(7) * t141;
t82 = (pkin(6) - t126) * t108;
t75 = (-t147 * t107 - pkin(6) + t143) * t108;
t1 = [MDP(1) + (t78 ^ 2 + t79 ^ 2 + t82 ^ 2) * MDP(21) + (t83 * t148 + t137) * t84 + t130 * t111 ^ 2 + 0.2e1 * (t83 * MDP(27) + t84 * MDP(28)) * t75 + (t79 * MDP(18) - t78 * MDP(20) + pkin(1) * MDP(9) - t152) * t151 + (-0.2e1 * pkin(1) * MDP(10) + (MDP(5) - t118) * t151 + (-t107 * t78 + t110 * t79) * t153 + 0.2e1 * (t107 * MDP(18) - t110 * MDP(20)) * t82 + (t104 * MDP(11) - 0.2e1 * t107 * t132 + MDP(4) + 0.2e1 * (t107 * MDP(16) + t110 * MDP(17)) * pkin(6)) * t108) * t108; t98 * MDP(16) + (-t110 * t82 + t98) * MDP(18) + t124 * MDP(19) - t107 * t82 * MDP(20) + (t124 * pkin(7) + t82 * t93) * MDP(21) + t89 * t137 + (-t83 * t89 - t84 * t88) * MDP(23) + (t75 * t88 + t83 * t85) * MDP(27) + (t75 * t89 + t84 * t85) * MDP(28) + (-pkin(6) * MDP(10) + MDP(7) + (t129 * pkin(7) - MDP(14)) * t110 + t113) * t111 + (MDP(6) - pkin(6) * MDP(9) + (-t102 + t104) * MDP(12) + (-pkin(6) * MDP(16) + t122) * t110 + (t110 * MDP(11) + pkin(6) * MDP(17) - t123) * t107) * t108; MDP(8) + t102 * MDP(11) + (t133 * pkin(7) ^ 2 + t93 ^ 2) * MDP(21) + t88 * MDP(27) * t149 + t133 * pkin(7) * t153 + (MDP(22) * t89 + MDP(28) * t149 + t88 * t148) * t89 + 0.2e1 * t123 * t110 + 0.2e1 * (t122 + t132) * t107; (-0.2e1 * t101 + t80) * MDP(18) + t81 * MDP(20) + (-pkin(3) * t79 + qJ(4) * t78) * MDP(21) + (-MDP(15) + t116 - 0.2e1 * t131) * t111 + (t127 * MDP(19) + t118) * t108 + t152; t110 * MDP(14) + t126 * MDP(19) + ((MDP(21) * qJ(4) - t129) * t110 + (-MDP(16) - MDP(18) - t145) * t107) * pkin(7) - t113; 0.2e1 * pkin(3) * MDP(18) + 0.2e1 * t131 + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(21) + 0.2e1 * t135 + 0.2e1 * t134 + t130; MDP(19) * t140 + t79 * MDP(21) + t115 * t111; (MDP(21) * pkin(7) + MDP(19)) * t107; -t115 - t145; MDP(21); t111 * MDP(26) + t154; t114; t116; t117; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
