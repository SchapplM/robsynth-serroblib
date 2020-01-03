% Calculate joint inertia matrix for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR9_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:40
% EndTime: 2019-12-31 21:24:42
% DurationCPUTime: 0.50s
% Computational Cost: add. (642->134), mult. (1299->206), div. (0->0), fcn. (1347->8), ass. (0->67)
t137 = sin(qJ(3));
t140 = cos(qJ(3));
t136 = sin(qJ(5));
t139 = cos(qJ(5));
t161 = -qJ(4) - pkin(7);
t123 = t161 * t137;
t124 = t161 * t140;
t134 = sin(pkin(9));
t135 = cos(pkin(9));
t104 = t135 * t123 + t124 * t134;
t116 = t134 * t140 + t135 * t137;
t92 = -pkin(8) * t116 + t104;
t105 = t134 * t123 - t135 * t124;
t115 = -t134 * t137 + t135 * t140;
t93 = pkin(8) * t115 + t105;
t98 = -t139 * t115 + t116 * t136;
t99 = t115 * t136 + t116 * t139;
t147 = t99 * MDP(22) - t98 * MDP(23) + (-t136 * t93 + t139 * t92) * MDP(25) - (t136 * t92 + t139 * t93) * MDP(26);
t171 = t147 - (MDP(16) * t137 + MDP(17) * t140) * pkin(7) + t137 * MDP(13) + t140 * MDP(14);
t138 = sin(qJ(2));
t110 = t116 * t138;
t155 = t138 * t140;
t157 = t137 * t138;
t111 = -t134 * t157 + t135 * t155;
t90 = t139 * t110 + t111 * t136;
t91 = -t110 * t136 + t111 * t139;
t160 = t91 * MDP(22) - t90 * MDP(23);
t141 = cos(qJ(2));
t122 = -pkin(2) * t141 - pkin(7) * t138 - pkin(1);
t117 = t140 * t122;
t164 = pkin(6) * t137;
t100 = -qJ(4) * t155 + t117 + (-pkin(3) - t164) * t141;
t162 = pkin(6) * t141;
t149 = t140 * t162;
t103 = t149 + (-qJ(4) * t138 + t122) * t137;
t86 = t135 * t100 - t103 * t134;
t82 = -pkin(4) * t141 - pkin(8) * t111 + t86;
t87 = t134 * t100 + t135 * t103;
t85 = -pkin(8) * t110 + t87;
t77 = -t136 * t85 + t139 * t82;
t78 = t136 * t82 + t139 * t85;
t170 = t77 * MDP(25) - t78 * MDP(26) + t160;
t168 = 2 * MDP(18);
t167 = -2 * MDP(21);
t166 = 0.2e1 * MDP(26);
t165 = pkin(3) * t134;
t163 = pkin(6) * t140;
t159 = MDP(20) * t99;
t158 = MDP(25) * t98;
t156 = t137 * t140;
t121 = pkin(3) * t157 + t138 * pkin(6);
t127 = pkin(3) * t135 + pkin(4);
t152 = MDP(25) * (t127 * t139 - t136 * t165);
t151 = MDP(26) * (t127 * t136 + t139 * t165);
t150 = MDP(15) + MDP(24);
t128 = -pkin(3) * t140 - pkin(2);
t148 = MDP(12) * t156;
t146 = MDP(13) * t140 - MDP(14) * t137;
t143 = MDP(24) - t151 + t152;
t132 = t140 ^ 2;
t131 = t138 ^ 2;
t130 = t137 ^ 2;
t109 = -pkin(4) * t115 + t128;
t108 = t122 * t137 + t149;
t107 = -t137 * t162 + t117;
t101 = pkin(4) * t110 + t121;
t1 = [MDP(1) - 0.2e1 * pkin(1) * t138 * MDP(10) + (t121 ^ 2 + t86 ^ 2 + t87 ^ 2) * MDP(19) + (MDP(20) * t91 + t90 * t167) * t91 + t150 * t141 ^ 2 + (t132 * MDP(11) + MDP(4) - 0.2e1 * t148) * t131 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t146) * t138 - t160) * t141 + 0.2e1 * (-t107 * t141 + t131 * t164) * MDP(16) + 0.2e1 * (t108 * t141 + t131 * t163) * MDP(17) + (-t110 * t87 - t111 * t86) * t168 + 0.2e1 * (t101 * t90 - t141 * t77) * MDP(25) + (t101 * t91 + t141 * t78) * t166; (-t104 * t111 - t105 * t110 + t115 * t87 - t116 * t86) * MDP(18) + (t104 * t86 + t105 * t87 + t121 * t128) * MDP(19) + t91 * t159 + (-t90 * t99 - t91 * t98) * MDP(21) + (t101 * t98 + t109 * t90) * MDP(25) + (t101 * t99 + t109 * t91) * MDP(26) + (-pkin(6) * MDP(10) + MDP(7) - t171) * t141 + (MDP(6) - pkin(6) * MDP(9) + MDP(11) * t156 + (-t130 + t132) * MDP(12) + (-pkin(2) * t137 - t163) * MDP(16) + (-pkin(2) * t140 + t164) * MDP(17)) * t138; MDP(8) + t130 * MDP(11) + 0.2e1 * t148 + (-t104 * t116 + t105 * t115) * t168 + (t104 ^ 2 + t105 ^ 2 + t128 ^ 2) * MDP(19) + 0.2e1 * t109 * t158 + 0.2e1 * (MDP(16) * t140 - MDP(17) * t137) * pkin(2) + (t109 * t166 + t98 * t167 + t159) * t99; t107 * MDP(16) - t108 * MDP(17) + (-MDP(15) - t143) * t141 + t146 * t138 + ((-t110 * t134 - t111 * t135) * MDP(18) + (t134 * t87 + t135 * t86) * MDP(19)) * pkin(3) + t170; ((t115 * t134 - t116 * t135) * MDP(18) + (t104 * t135 + t105 * t134) * MDP(19)) * pkin(3) + t171; (t134 ^ 2 + t135 ^ 2) * MDP(19) * pkin(3) ^ 2 + 0.2e1 * t152 - 0.2e1 * t151 + t150; MDP(19) * t121 + t90 * MDP(25) + t91 * MDP(26); MDP(19) * t128 + MDP(26) * t99 + t158; 0; MDP(19); -t141 * MDP(24) + t170; t147; t143; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
