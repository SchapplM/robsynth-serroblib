% Calculate joint inertia matrix for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPPR3_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:25
% EndTime: 2019-03-09 02:45:26
% DurationCPUTime: 0.25s
% Computational Cost: add. (229->101), mult. (354->125), div. (0->0), fcn. (252->6), ass. (0->45)
t73 = sin(pkin(9));
t60 = pkin(1) * t73 + pkin(7);
t103 = -qJ(5) + t60;
t80 = -pkin(3) - pkin(4);
t74 = cos(pkin(9));
t61 = -t74 * pkin(1) - pkin(2);
t79 = cos(qJ(3));
t55 = t103 * t79;
t102 = t55 * t79;
t76 = sin(qJ(6));
t78 = cos(qJ(6));
t101 = t76 * t78;
t77 = sin(qJ(3));
t94 = t78 * MDP(25);
t100 = (MDP(16) + t94) * t77;
t63 = t77 * qJ(4);
t99 = t79 * pkin(3) + t63;
t70 = t77 ^ 2;
t72 = t79 ^ 2;
t98 = t70 + t72;
t97 = pkin(3) * MDP(15);
t96 = t60 ^ 2 * MDP(15);
t95 = t76 * MDP(26);
t93 = MDP(13) - MDP(18);
t92 = MDP(14) - MDP(11);
t91 = MDP(15) + MDP(19);
t90 = 0.2e1 * t79;
t89 = MDP(21) * t101;
t53 = -t99 + t61;
t88 = t80 * MDP(19) + MDP(17);
t52 = t79 * pkin(4) - t53;
t87 = MDP(12) - t88;
t86 = -MDP(22) * t78 + MDP(23) * t76;
t85 = t94 - t95;
t84 = MDP(25) * t76 + MDP(26) * t78;
t83 = -t76 * MDP(22) - t78 * MDP(23) - t84 * (-pkin(8) + t80);
t81 = qJ(4) ^ 2;
t75 = qJ(4) + pkin(5);
t71 = t78 ^ 2;
t69 = t76 ^ 2;
t54 = t103 * t77;
t51 = pkin(5) * t77 + pkin(8) * t79 + t52;
t50 = t51 * t76 + t54 * t78;
t49 = t51 * t78 - t54 * t76;
t1 = [MDP(1) + t53 ^ 2 * MDP(15) + (t52 ^ 2 + t54 ^ 2 + t55 ^ 2) * MDP(19) + (t73 ^ 2 + t74 ^ 2) * MDP(4) * pkin(1) ^ 2 + (t71 * MDP(20) - 0.2e1 * t89 + t96) * t72 + (MDP(24) + MDP(5) + t96) * t70 + (-t61 * MDP(10) - t53 * MDP(12) - t52 * MDP(17)) * t90 + (0.2e1 * t61 * MDP(11) - 0.2e1 * t53 * MDP(14) + 0.2e1 * t52 * MDP(16) + (MDP(6) + t86) * t90) * t77 + 0.2e1 * (-t54 * t77 - t102) * MDP(18) + 0.2e1 * (-t76 * t102 + t49 * t77) * MDP(25) + 0.2e1 * (-t78 * t102 - t50 * t77) * MDP(26) + 0.2e1 * t98 * MDP(13) * t60; (-t54 * t79 + t55 * t77) * MDP(19); t91 * t98 + MDP(4); t88 * t54 + (qJ(4) * MDP(19) + MDP(16) + t85) * t55 + (-pkin(3) * MDP(13) - t80 * MDP(18) + MDP(7) + t83) * t77 + (MDP(8) + MDP(20) * t101 + (-t69 + t71) * MDP(21) - t84 * t75 + t93 * qJ(4)) * t79 + ((qJ(4) * MDP(15) + t92) * t79 + (-MDP(10) - MDP(12) - t97) * t77) * t60; t99 * MDP(15) + t63 * MDP(19) + (MDP(10) + t87) * t79 + (t92 - t95) * t77 + t100; MDP(9) + 0.2e1 * pkin(3) * MDP(12) + (pkin(3) ^ 2 + t81) * MDP(15) + 0.2e1 * t80 * MDP(17) + (t80 ^ 2 + t81) * MDP(19) + t69 * MDP(20) + 0.2e1 * t89 + 0.2e1 * t85 * t75 + 0.2e1 * (MDP(14) + MDP(16)) * qJ(4); t54 * MDP(19) + (MDP(15) * t60 - t84 + t93) * t77; -t91 * t79; -t87 - t97; t91; -t79 * MDP(17) + MDP(19) * t52 - t77 * t95 + t100; 0; 0; 0; MDP(19); MDP(24) * t77 + t49 * MDP(25) - t50 * MDP(26) + t86 * t79; t84 * t79; t83; -t84; t85; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
