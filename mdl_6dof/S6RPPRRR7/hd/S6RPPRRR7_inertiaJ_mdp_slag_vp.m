% Calculate joint inertia matrix for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPPRRR7_inertiaJ_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:49
% EndTime: 2019-03-09 02:33:51
% DurationCPUTime: 0.36s
% Computational Cost: add. (610->102), mult. (1048->142), div. (0->0), fcn. (1214->8), ass. (0->57)
t122 = sin(qJ(6));
t125 = cos(qJ(6));
t133 = t125 * MDP(30) - t122 * MDP(31);
t160 = MDP(23) + t133;
t119 = sin(pkin(10));
t120 = cos(pkin(10));
t124 = sin(qJ(4));
t126 = cos(qJ(4));
t101 = t126 * t119 + t124 * t120;
t102 = -t124 * t119 + t126 * t120;
t123 = sin(qJ(5));
t155 = cos(qJ(5));
t92 = t123 * t101 - t155 * t102;
t159 = t92 ^ 2;
t148 = t122 * MDP(27) + t125 * MDP(28);
t158 = MDP(30) * t122 + MDP(31) * t125;
t107 = t119 * pkin(3) + qJ(2);
t96 = t101 * pkin(4) + t107;
t157 = 0.2e1 * t96;
t156 = 0.2e1 * t107;
t121 = -pkin(1) - qJ(3);
t153 = -pkin(7) + t121;
t103 = t153 * t119;
t104 = t153 * t120;
t140 = -t124 * t103 + t126 * t104;
t81 = -t102 * pkin(8) + t140;
t135 = -t126 * t103 - t124 * t104;
t82 = -t101 * pkin(8) - t135;
t77 = t123 * t82 - t155 * t81;
t74 = t77 * t122;
t152 = t77 * t125;
t151 = t122 * t125;
t130 = t101 * t155 + t123 * t102;
t150 = t130 * MDP(29);
t149 = t92 * MDP(24);
t105 = t119 ^ 2 + t120 ^ 2;
t145 = t101 * MDP(16);
t142 = MDP(26) * t151;
t117 = t122 ^ 2;
t141 = t117 * MDP(25) + MDP(22) + 0.2e1 * t142;
t139 = pkin(5) * t92 - pkin(9) * t130;
t108 = t123 * pkin(4) + pkin(9);
t109 = -pkin(4) * t155 - pkin(5);
t137 = -t108 * t130 - t109 * t92;
t136 = t119 * MDP(7) + t120 * MDP(8);
t134 = MDP(27) * t125 - MDP(28) * t122;
t131 = -MDP(24) * t130 - t160 * t92;
t118 = t125 ^ 2;
t78 = t123 * t81 + t155 * t82;
t129 = -t77 * MDP(23) - t78 * MDP(24) + (-(-t117 + t118) * MDP(26) - MDP(25) * t151 - MDP(20)) * t92 + (-MDP(21) + t148) * t130;
t128 = (MDP(23) * t155 - t123 * MDP(24)) * pkin(4);
t127 = qJ(2) ^ 2;
t100 = t105 * t121;
t79 = pkin(5) * t130 + pkin(9) * t92 + t96;
t73 = t122 * t79 + t125 * t78;
t72 = -t122 * t78 + t125 * t79;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2) + t127) * MDP(6) + (t105 * t121 ^ 2 + t127) * MDP(10) + t145 * t156 - t149 * t157 + (t118 * MDP(25) + MDP(18) - 0.2e1 * t142) * t159 + (MDP(11) * t102 - 0.2e1 * t101 * MDP(12) + MDP(17) * t156) * t102 + (MDP(23) * t157 + t150) * t130 - 0.2e1 * t100 * MDP(9) + 0.2e1 * (t130 * t72 - t74 * t92) * MDP(30) + 0.2e1 * (-t130 * t73 - t152 * t92) * MDP(31) + 0.2e1 * (MDP(5) + t136) * qJ(2) - 0.2e1 * (-MDP(19) + t134) * t92 * t130; t100 * MDP(10) - pkin(1) * MDP(6) - t105 * MDP(9) + MDP(4) + t158 * (-t130 ^ 2 - t159); t105 * MDP(10) + MDP(6); qJ(2) * MDP(10) + t102 * MDP(17) + t160 * t130 + t136 + t145 - t149; 0; MDP(10); t102 * MDP(13) - t101 * MDP(14) + t140 * MDP(16) + t135 * MDP(17) + (t122 * t137 - t152) * MDP(30) + (t125 * t137 + t74) * MDP(31) + t129; t102 * MDP(16) - t101 * MDP(17) + t131; 0; -0.2e1 * t109 * t133 + MDP(15) + 0.2e1 * t128 + t141; (t122 * t139 - t152) * MDP(30) + (t125 * t139 + t74) * MDP(31) + t129; t131; 0; t128 + t141 + t133 * (pkin(5) - t109); 0.2e1 * pkin(5) * t133 + t141; t72 * MDP(30) - t73 * MDP(31) - t134 * t92 + t150; -t158 * t130; t133; -t108 * t158 + t148; -pkin(9) * t158 + t148; MDP(29);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
