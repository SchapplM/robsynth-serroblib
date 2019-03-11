% Calculate joint inertia matrix for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPPR2_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:49
% EndTime: 2019-03-09 02:42:50
% DurationCPUTime: 0.30s
% Computational Cost: add. (416->97), mult. (717->144), div. (0->0), fcn. (728->8), ass. (0->51)
t118 = MDP(13) + MDP(17);
t85 = sin(pkin(10));
t87 = cos(pkin(10));
t90 = sin(qJ(3));
t92 = cos(qJ(3));
t71 = t85 * t90 - t87 * t92;
t69 = t71 ^ 2;
t73 = t85 * t92 + t87 * t90;
t70 = t73 ^ 2;
t88 = cos(pkin(9));
t82 = -pkin(1) * t88 - pkin(2);
t75 = -pkin(3) * t92 + t82;
t96 = -qJ(5) * t73 + t75;
t62 = pkin(4) * t71 + t96;
t117 = -0.2e1 * t62;
t86 = sin(pkin(9));
t80 = pkin(1) * t86 + pkin(7);
t113 = qJ(4) + t80;
t105 = t113 * t90;
t68 = t113 * t92;
t65 = -t85 * t105 + t87 * t68;
t61 = -t71 * pkin(5) + t65;
t116 = t61 * t71;
t89 = sin(qJ(6));
t91 = cos(qJ(6));
t115 = t89 * t91;
t114 = MDP(13) * pkin(3);
t77 = pkin(3) * t85 + qJ(5);
t112 = t77 * MDP(17);
t81 = -pkin(3) * t87 - pkin(4);
t111 = t81 * MDP(17);
t110 = t89 * MDP(23);
t109 = t91 * MDP(24);
t108 = t92 * MDP(10);
t107 = MDP(19) * t115;
t63 = t87 * t105 + t68 * t85;
t106 = t63 ^ 2 + t65 ^ 2;
t104 = MDP(15) + t111;
t101 = MDP(20) * t89 + MDP(21) * t91;
t100 = MDP(23) * t91 - t89 * MDP(24);
t99 = t109 + t110;
t98 = MDP(14) + t100;
t97 = -MDP(16) - t99;
t95 = t91 * MDP(20) - t89 * MDP(21) + t100 * (-pkin(8) + t81);
t84 = t91 ^ 2;
t83 = t89 ^ 2;
t60 = pkin(5) * t73 + t63;
t59 = (pkin(4) + pkin(8)) * t71 + t96;
t58 = t59 * t91 + t60 * t89;
t57 = -t59 * t89 + t60 * t91;
t1 = [MDP(1) - 0.2e1 * t82 * t108 + (t75 ^ 2 + t106) * MDP(13) + t73 * MDP(16) * t117 + (t62 ^ 2 + t106) * MDP(17) + t70 * MDP(22) + (t86 ^ 2 + t88 ^ 2) * MDP(4) * pkin(1) ^ 2 + (t83 * MDP(18) + 0.2e1 * t107) * t69 + (0.2e1 * t82 * MDP(11) + MDP(5) * t90 + 0.2e1 * t92 * MDP(6)) * t90 + (MDP(15) * t117 + 0.2e1 * t101 * t73) * t71 + 0.2e1 * (-t91 * t116 + t57 * t73) * MDP(23) + 0.2e1 * (t89 * t116 - t58 * t73) * MDP(24) + 0.2e1 * (MDP(12) + MDP(14)) * (t63 * t73 - t65 * t71); t118 * (t63 * t71 + t65 * t73); MDP(4) + t118 * (t69 + t70); t90 * MDP(7) + t92 * MDP(8) + t63 * MDP(15) + t65 * MDP(16) + (t63 * t81 + t65 * t77) * MDP(17) + (-t90 * MDP(10) - t92 * MDP(11)) * t80 + t99 * t61 + (t81 * MDP(14) + t95) * t73 + (MDP(18) * t115 + (-t83 + t84) * MDP(19) - t98 * t77) * t71 + ((-t71 * t85 - t73 * t87) * MDP(12) + (-t63 * t87 + t65 * t85) * MDP(13)) * pkin(3); t108 - t90 * MDP(11) + (-t87 * t114 + t104) * t71 + (t85 * t114 + t112 - t97) * t73; -0.2e1 * t107 + t84 * MDP(18) + MDP(9) + (t85 ^ 2 + t87 ^ 2) * MDP(13) * pkin(3) ^ 2 + (0.2e1 * MDP(15) + t111) * t81 + (0.2e1 * MDP(16) + 0.2e1 * t109 + 0.2e1 * t110 + t112) * t77; MDP(13) * t75 - t71 * MDP(15) + MDP(17) * t62 + t97 * t73; 0; 0; t118; t63 * MDP(17) + t98 * t73; t71 * MDP(17); t104; 0; MDP(17); MDP(22) * t73 + t57 * MDP(23) - t58 * MDP(24) + t101 * t71; t100 * t71; t95; -t99; t100; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
