% Calculate joint inertia matrix for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPPRR1_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:15
% EndTime: 2019-03-08 19:16:16
% DurationCPUTime: 0.25s
% Computational Cost: add. (317->86), mult. (697->141), div. (0->0), fcn. (785->12), ass. (0->57)
t92 = sin(qJ(6));
t95 = cos(qJ(6));
t102 = MDP(22) * t95 - MDP(23) * t92;
t100 = MDP(15) + t102;
t117 = cos(qJ(5));
t86 = sin(pkin(12));
t89 = cos(pkin(12));
t93 = sin(qJ(5));
t74 = -t117 * t89 + t86 * t93;
t121 = t100 * t74;
t87 = sin(pkin(11));
t88 = sin(pkin(6));
t90 = cos(pkin(11));
t94 = sin(qJ(2));
t96 = cos(qJ(2));
t68 = (t87 * t94 - t90 * t96) * t88;
t120 = t68 ^ 2;
t81 = -pkin(2) * t90 - pkin(3);
t76 = -pkin(4) * t89 + t81;
t119 = 0.2e1 * t76;
t79 = pkin(2) * t87 + qJ(4);
t118 = pkin(8) + t79;
t116 = MDP(5) * pkin(2);
t75 = t117 * t86 + t93 * t89;
t115 = t75 * t92;
t114 = t75 * t95;
t113 = t92 * t95;
t112 = t86 ^ 2 + t89 ^ 2;
t111 = MDP(9) * t81;
t110 = t86 * MDP(7);
t109 = t89 * MDP(6);
t108 = MDP(21) * t74;
t107 = t75 * MDP(16);
t106 = MDP(18) * t113;
t105 = t112 * MDP(9);
t103 = MDP(19) * t95 - MDP(20) * t92;
t101 = -MDP(22) * t92 - MDP(23) * t95;
t99 = t107 - t109 + t110 + t111;
t98 = t92 * MDP(19) + t95 * MDP(20) + t101 * pkin(9);
t91 = cos(pkin(6));
t85 = t95 ^ 2;
t84 = t92 ^ 2;
t72 = t118 * t89;
t71 = t118 * t86;
t70 = (t87 * t96 + t90 * t94) * t88;
t67 = t70 * t89 + t86 * t91;
t66 = -t70 * t86 + t89 * t91;
t65 = t117 * t72 - t93 * t71;
t64 = t117 * t71 + t93 * t72;
t63 = pkin(5) * t74 - pkin(9) * t75 + t76;
t62 = t117 * t67 + t93 * t66;
t61 = -t117 * t66 + t93 * t67;
t60 = t63 * t92 + t65 * t95;
t59 = t63 * t95 - t65 * t92;
t58 = t62 * t95 + t68 * t92;
t57 = -t62 * t92 + t68 * t95;
t1 = [MDP(1) + (t70 ^ 2 + t91 ^ 2 + t120) * MDP(5) + (t66 ^ 2 + t67 ^ 2 + t120) * MDP(9); t70 * t87 * t116 + (t61 * t115 + t57 * t74) * MDP(22) + (t61 * t114 - t58 * t74) * MDP(23) + (MDP(3) * t96 - MDP(4) * t94) * t88 + (t74 * MDP(15) - t90 * t116 + t99) * t68 + (MDP(9) * t79 + MDP(8)) * (-t66 * t86 + t67 * t89); t107 * t119 + MDP(2) + (t87 ^ 2 + t90 ^ 2) * MDP(5) * pkin(2) ^ 2 + (-0.2e1 * t109 + 0.2e1 * t110 + t111) * t81 + (t85 * MDP(17) + MDP(10) - 0.2e1 * t106) * t75 ^ 2 + (MDP(15) * t119 + t108 + 0.2e1 * (-MDP(11) + t103) * t75) * t74 + 0.2e1 * (t64 * t115 + t59 * t74) * MDP(22) + 0.2e1 * (t64 * t114 - t60 * t74) * MDP(23) + (0.2e1 * t112 * MDP(8) + t105 * t79) * t79; t91 * MDP(5) + (t66 * t89 + t67 * t86) * MDP(9); 0; MDP(5) + t105; t68 * MDP(9); t99 + t121; 0; MDP(9); -MDP(16) * t62 - t100 * t61; -t65 * MDP(16) - t100 * t64 + (-MDP(13) + t98) * t74 + (MDP(12) + MDP(17) * t113 + (-t84 + t85) * MDP(18) + t101 * pkin(5)) * t75; -t107 - t121; 0; MDP(17) * t84 + 0.2e1 * pkin(5) * t102 + MDP(14) + 0.2e1 * t106; MDP(22) * t57 - MDP(23) * t58; t59 * MDP(22) - t60 * MDP(23) + t103 * t75 + t108; t101 * t75; t102; t98; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
