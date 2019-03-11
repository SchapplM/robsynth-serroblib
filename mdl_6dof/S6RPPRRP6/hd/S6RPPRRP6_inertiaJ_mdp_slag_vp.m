% Calculate joint inertia matrix for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRRP6_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:19
% EndTime: 2019-03-09 02:11:20
% DurationCPUTime: 0.34s
% Computational Cost: add. (324->120), mult. (506->158), div. (0->0), fcn. (376->4), ass. (0->48)
t81 = cos(qJ(4));
t116 = -t81 * MDP(16) - MDP(8);
t115 = MDP(27) * pkin(8) + MDP(25);
t100 = MDP(23) - MDP(26);
t101 = MDP(22) + MDP(24);
t78 = sin(qJ(5));
t80 = cos(qJ(5));
t114 = -t100 * t80 - t101 * t78;
t113 = t100 * t78 - t101 * t80 - MDP(15);
t112 = 2 * pkin(5);
t75 = -pkin(7) + qJ(2);
t110 = t75 * t78;
t109 = t75 * t80;
t79 = sin(qJ(4));
t108 = t78 * t79;
t76 = pkin(1) + qJ(3);
t67 = pkin(4) * t79 - pkin(8) * t81 + t76;
t61 = t109 * t79 + t67 * t78;
t71 = t78 ^ 2;
t73 = t80 ^ 2;
t107 = t71 + t73;
t94 = -pkin(5) * t80 - qJ(6) * t78;
t68 = -pkin(4) + t94;
t104 = t68 * MDP(27);
t103 = t80 * MDP(18);
t99 = t107 * MDP(25);
t98 = t107 * MDP(27);
t97 = -MDP(27) * pkin(5) - MDP(24);
t95 = 0.2e1 * qJ(6) * MDP(26) + MDP(21);
t93 = -pkin(5) * t78 + qJ(6) * t80;
t58 = qJ(6) * t79 + t61;
t64 = t80 * t67;
t59 = -t64 + (-pkin(5) + t110) * t79;
t92 = t58 * t80 + t59 * t78;
t91 = -t58 * t78 + t59 * t80;
t90 = -MDP(22) + t97;
t89 = MDP(27) * qJ(6) - t100;
t88 = t80 * MDP(19) - t78 * MDP(20);
t87 = t78 * MDP(19) + t80 * MDP(20);
t86 = (-t108 * t75 + t64) * MDP(22) - t61 * MDP(23);
t85 = t78 * MDP(24) - t80 * MDP(26);
t84 = t78 * t90 + t80 * t89;
t82 = (qJ(2) ^ 2);
t74 = t81 ^ 2;
t72 = t79 ^ 2;
t70 = t73 * t81;
t62 = (-t75 - t93) * t81;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2 + t82) * MDP(6)) + (t76 ^ 2 + t82) * MDP(9) + t72 * MDP(21) + (t58 ^ 2 + t59 ^ 2 + t62 ^ 2) * MDP(27) + 0.2e1 * (-t59 * MDP(24) + t58 * MDP(26) + t86) * t79 + (t73 * MDP(17) - 0.2e1 * t78 * t103 + MDP(10) + 0.2e1 * (-t78 * MDP(22) - t80 * MDP(23)) * t75) * t74 + (2 * (MDP(5) + MDP(7)) * qJ(2)) + 0.2e1 * (t79 * MDP(15) - t116) * t76 + 0.2e1 * (t91 * MDP(25) + t85 * t62 + (-MDP(11) + t88) * t79) * t81; MDP(4) - (pkin(1) * MDP(6)) - t76 * MDP(9) + (t71 * t81 + t70) * MDP(25) + t91 * MDP(27) + t113 * t79 + t116; MDP(6) + MDP(9) + t98; MDP(7) + qJ(2) * MDP(9) + (-t62 * t81 + t79 * t92) * MDP(27) + t114 * (t72 + t74); 0; MDP(9) + (t107 * t72 + t74) * MDP(27); t70 * MDP(18) + (-t80 * MDP(24) - t78 * MDP(26) + t104) * t62 + (-t75 * MDP(16) + pkin(8) * t114 - MDP(13) + t87) * t79 + (MDP(12) + t75 * MDP(15) + t80 * t78 * MDP(17) - t71 * MDP(18) + (-pkin(4) * t78 + t109) * MDP(22) + (-pkin(4) * t80 - t110) * MDP(23) + t85 * t68) * t81 + t115 * t92; 0; (pkin(8) * t98 - MDP(16) + t99) * t79 + (-t104 - t113) * t81; MDP(14) + t71 * MDP(17) + (pkin(8) ^ 2 * t107 + t68 ^ 2) * MDP(27) + 0.2e1 * pkin(8) * t99 + 0.2e1 * (MDP(22) * pkin(4) - MDP(24) * t68) * t80 + 0.2e1 * (-MDP(23) * pkin(4) - MDP(26) * t68 + t103) * t78; t64 * MDP(24) + t61 * MDP(26) + (-pkin(5) * t59 + qJ(6) * t58) * MDP(27) + ((t112 - t110) * MDP(24) + t95) * t79 + (MDP(25) * t94 + t88) * t81 + t86; -t78 * t89 + t80 * t90; t84 * t79; MDP(25) * t93 + pkin(8) * t84 + t87; MDP(24) * t112 + ((pkin(5) ^ 2) + qJ(6) ^ 2) * MDP(27) + t95; MDP(25) * t80 * t81 - MDP(24) * t79 + MDP(27) * t59; t80 * MDP(27); MDP(27) * t108; t115 * t78; t97; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
