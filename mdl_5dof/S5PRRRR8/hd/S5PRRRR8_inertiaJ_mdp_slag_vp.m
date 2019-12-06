% Calculate joint inertia matrix for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR8_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:43
% EndTime: 2019-12-05 17:16:44
% DurationCPUTime: 0.29s
% Computational Cost: add. (269->88), mult. (573->134), div. (0->0), fcn. (613->10), ass. (0->53)
t103 = cos(qJ(5));
t99 = sin(qJ(5));
t126 = t99 * MDP(21) + t103 * MDP(22);
t110 = t103 * MDP(24) - t99 * MDP(25);
t105 = cos(qJ(3));
t91 = -t105 * pkin(3) - pkin(2);
t130 = 0.2e1 * t91;
t129 = pkin(7) + pkin(8);
t100 = sin(qJ(4));
t101 = sin(qJ(3));
t104 = cos(qJ(4));
t84 = t100 * t105 + t104 * t101;
t127 = t84 * t99;
t102 = sin(qJ(2));
t97 = sin(pkin(5));
t125 = t102 * t97;
t124 = t103 * t99;
t106 = cos(qJ(2));
t123 = t106 * t97;
t85 = t129 * t101;
t86 = t129 * t105;
t73 = t100 * t86 + t104 * t85;
t122 = t73 * t103;
t121 = MDP(18) * t84;
t83 = t100 * t101 - t104 * t105;
t120 = t83 * MDP(23);
t118 = MDP(10) * t105;
t115 = MDP(20) * t124;
t95 = t99 ^ 2;
t116 = t95 * MDP(19) + MDP(16) + 0.2e1 * t115;
t114 = -pkin(4) * t84 - pkin(9) * t83;
t89 = t100 * pkin(3) + pkin(9);
t90 = -t104 * pkin(3) - pkin(4);
t113 = -t83 * t89 + t84 * t90;
t112 = MDP(21) * t103 - MDP(22) * t99;
t111 = -MDP(24) * t99 - MDP(25) * t103;
t98 = cos(pkin(5));
t78 = -t101 * t125 + t98 * t105;
t79 = t98 * t101 + t105 * t125;
t66 = t100 * t79 - t104 * t78;
t67 = t100 * t78 + t104 * t79;
t109 = -t67 * MDP(18) + (-MDP(17) - t110) * t66;
t108 = (t104 * MDP(17) - t100 * MDP(18)) * pkin(3);
t74 = -t100 * t85 + t104 * t86;
t96 = t103 ^ 2;
t107 = -t73 * MDP(17) - t74 * MDP(18) + ((-t95 + t96) * MDP(20) + MDP(19) * t124 + MDP(14)) * t84 + (-MDP(15) + t126) * t83;
t70 = t73 * t99;
t69 = t83 * pkin(4) - t84 * pkin(9) + t91;
t62 = t103 * t67 - t99 * t123;
t61 = -t103 * t123 - t99 * t67;
t60 = t103 * t74 + t99 * t69;
t59 = t103 * t69 - t99 * t74;
t1 = [MDP(1); (t66 * t127 + t61 * t83) * MDP(24) + (t66 * t103 * t84 - t62 * t83) * MDP(25) + (-t102 * MDP(4) + (-MDP(11) * t101 - MDP(17) * t83 + MDP(3) + t118 - t121) * t106) * t97; 0.2e1 * pkin(2) * t118 + t121 * t130 + MDP(2) + (MDP(17) * t130 + t120 + 0.2e1 * (-MDP(13) + t112) * t84) * t83 + 0.2e1 * (t73 * t127 + t59 * t83) * MDP(24) + 0.2e1 * (t84 * t122 - t60 * t83) * MDP(25) + (t96 * MDP(19) + MDP(12) - 0.2e1 * t115) * t84 ^ 2 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t101 + 0.2e1 * t105 * MDP(6)) * t101; t78 * MDP(10) - t79 * MDP(11) + t109; t101 * MDP(7) + t105 * MDP(8) + (t113 * t99 - t122) * MDP(24) + (t113 * t103 + t70) * MDP(25) + (-t101 * MDP(10) - t105 * MDP(11)) * pkin(7) + t107; -0.2e1 * t110 * t90 + MDP(9) + 0.2e1 * t108 + t116; t109; (t114 * t99 - t122) * MDP(24) + (t114 * t103 + t70) * MDP(25) + t107; t108 + t116 + t110 * (pkin(4) - t90); 0.2e1 * pkin(4) * t110 + t116; t61 * MDP(24) - t62 * MDP(25); t59 * MDP(24) - t60 * MDP(25) + t112 * t84 + t120; t111 * t89 + t126; t111 * pkin(9) + t126; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
