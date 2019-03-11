% Calculate joint inertia matrix for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR5_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:21:07
% EndTime: 2019-03-08 21:21:08
% DurationCPUTime: 0.47s
% Computational Cost: add. (450->140), mult. (893->198), div. (0->0), fcn. (904->10), ass. (0->68)
t169 = -2 * MDP(21);
t168 = pkin(4) + pkin(8);
t128 = sin(qJ(6));
t131 = cos(qJ(6));
t127 = -pkin(3) - qJ(5);
t132 = cos(qJ(3));
t129 = sin(qJ(3));
t149 = -t129 * qJ(4) - pkin(2);
t102 = t127 * t132 + t149;
t111 = t168 * t129;
t125 = cos(pkin(11));
t107 = t125 * t111;
t123 = sin(pkin(11));
t87 = t129 * pkin(5) + t107 + (pkin(9) * t132 - t102) * t123;
t157 = t125 * t132;
t90 = t125 * t102 + t123 * t111;
t88 = -pkin(9) * t157 + t90;
t105 = -t128 * t123 + t131 * t125;
t95 = t105 * t132;
t104 = t131 * t123 + t128 * t125;
t96 = t104 * t132;
t167 = (-t128 * t88 + t131 * t87) * MDP(25) - (t128 * t87 + t131 * t88) * MDP(26) - t96 * MDP(22) - t95 * MDP(23);
t164 = MDP(19) * t127 - MDP(18);
t126 = cos(pkin(6));
t124 = sin(pkin(6));
t130 = sin(qJ(2));
t159 = t124 * t130;
t100 = t126 * t129 + t132 * t159;
t163 = t100 ^ 2;
t115 = t123 * pkin(5) + qJ(4);
t162 = 0.2e1 * t115;
t161 = -pkin(9) + t127;
t160 = MDP(15) * pkin(8);
t133 = cos(qJ(2));
t158 = t124 * t133;
t112 = t168 * t132;
t114 = t123 ^ 2 + t125 ^ 2;
t153 = MDP(20) * t105;
t152 = qJ(4) * MDP(19);
t151 = t104 * MDP(25);
t150 = -MDP(11) + MDP(14);
t148 = -pkin(3) * MDP(15) + MDP(13);
t147 = -MDP(10) + t148;
t89 = -t123 * t102 + t107;
t83 = t90 * t123 + t89 * t125;
t99 = -t126 * t132 + t129 * t159;
t93 = t123 * t158 + t99 * t125;
t94 = t99 * t123 - t125 * t158;
t84 = t94 * t123 + t93 * t125;
t144 = (-t128 * t94 + t131 * t93) * MDP(25) - (t128 * t93 + t131 * t94) * MDP(26);
t143 = t95 * MDP(25) - t96 * MDP(26);
t142 = t125 * MDP(16) - t123 * MDP(17);
t141 = t123 * MDP(16) + t125 * MDP(17);
t140 = t105 * MDP(25) - t104 * MDP(26);
t139 = MDP(12) + t142;
t138 = t112 * MDP(19) + t143;
t108 = t161 * t123;
t109 = t161 * t125;
t137 = t105 * MDP(22) - t104 * MDP(23) + (-t128 * t108 + t131 * t109) * MDP(25) - (t131 * t108 + t128 * t109) * MDP(26);
t136 = t105 * MDP(26) + t141 + t151;
t135 = pkin(8) ^ 2;
t134 = qJ(4) ^ 2;
t122 = t132 ^ 2;
t121 = t129 ^ 2;
t110 = -t132 * pkin(3) + t149;
t103 = t114 * t127;
t98 = pkin(5) * t157 + t112;
t1 = [MDP(1) + (t124 ^ 2 * t133 ^ 2 + t99 ^ 2 + t163) * MDP(15) + (t93 ^ 2 + t94 ^ 2 + t163) * MDP(19); (t93 * t89 + t94 * t90) * MDP(19) + t138 * t100 + (t99 * MDP(12) + t93 * MDP(16) - t94 * MDP(17) + t144) * t129 + ((t123 * t93 - t125 * t94) * MDP(18) + t139 * t100) * t132 + (t100 * t132 + t99 * t129) * t160 + (-t130 * MDP(4) + (-MDP(15) * t110 + MDP(3) + (MDP(10) - MDP(13)) * t132 + t150 * t129) * t133) * t124; MDP(2) + (t110 ^ 2 + t122 * t135) * MDP(15) + (t112 ^ 2 + t89 ^ 2 + t90 ^ 2) * MDP(19) - (-t96 * MDP(20) + t95 * t169) * t96 + (t135 * MDP(15) + MDP(24) + MDP(5)) * t121 + 0.2e1 * t143 * t98 + 0.2e1 * (t121 + t122) * MDP(12) * pkin(8) + 0.2e1 * (-pkin(2) * MDP(11) - t110 * MDP(14) + t89 * MDP(16) - t90 * MDP(17) + t132 * MDP(6) + t167) * t129 + 0.2e1 * ((t123 * t89 - t125 * t90) * MDP(18) + t142 * t112 + pkin(2) * MDP(10) + t110 * MDP(13)) * t132; t147 * t99 + ((MDP(15) + MDP(19)) * qJ(4) + t136 + t150) * t100 + t164 * t84; -t96 * t153 + (t96 * t104 - t105 * t95) * MDP(21) + (t98 * t104 + t115 * t95) * MDP(25) + (t98 * t105 - t115 * t96) * MDP(26) + (t141 + t152) * t112 + (t139 * qJ(4) + MDP(8)) * t132 + (-pkin(3) * MDP(12) + t142 * t127 + MDP(7) + t137) * t129 + ((qJ(4) * MDP(15) + t150) * t132 + t147 * t129) * pkin(8) + t164 * t83; MDP(9) - 0.2e1 * pkin(3) * MDP(13) + (pkin(3) ^ 2 + t134) * MDP(15) - 0.2e1 * t103 * MDP(18) + (t114 * t127 ^ 2 + t134) * MDP(19) + t151 * t162 + (MDP(26) * t162 + t104 * t169 + t153) * t105 + 0.2e1 * (MDP(14) + t141) * qJ(4); t99 * MDP(15) + t84 * MDP(19); t83 * MDP(19) + (t139 + t140 + t160) * t129; -t114 * MDP(18) + t103 * MDP(19) + t148; t114 * MDP(19) + MDP(15); t100 * MDP(19); t142 * t132 + t138; t136 + t152; 0; MDP(19); t144; t129 * MDP(24) + t167; t137; t140; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
