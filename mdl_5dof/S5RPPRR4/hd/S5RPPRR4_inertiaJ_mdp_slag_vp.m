% Calculate joint inertia matrix for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR4_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:59
% EndTime: 2022-01-23 09:17:00
% DurationCPUTime: 0.26s
% Computational Cost: add. (366->74), mult. (769->109), div. (0->0), fcn. (791->8), ass. (0->41)
t89 = sin(pkin(8));
t91 = cos(pkin(8));
t81 = -pkin(2) * t91 - t89 * qJ(3) - pkin(1);
t90 = cos(pkin(9));
t77 = t90 * t81;
t88 = sin(pkin(9));
t66 = -t90 * t89 * pkin(6) + t77 + (-qJ(2) * t88 - pkin(3)) * t91;
t113 = t88 * t89;
t109 = qJ(2) * t91;
t72 = t90 * t109 + t88 * t81;
t70 = -pkin(6) * t113 + t72;
t93 = sin(qJ(4));
t95 = cos(qJ(4));
t57 = t95 * t66 - t93 * t70;
t78 = -t93 * t88 + t95 * t90;
t75 = t78 * t89;
t55 = -t91 * pkin(4) - t75 * pkin(7) + t57;
t58 = t93 * t66 + t95 * t70;
t79 = t95 * t88 + t93 * t90;
t74 = t79 * t89;
t56 = -t74 * pkin(7) + t58;
t92 = sin(qJ(5));
t94 = cos(qJ(5));
t61 = t94 * t74 + t92 * t75;
t62 = -t92 * t74 + t94 * t75;
t118 = -(t92 * t55 + t94 * t56) * MDP(23) + (t94 * t55 - t92 * t56) * MDP(22) + t62 * MDP(19) - t61 * MDP(20);
t117 = t75 * MDP(12) - t74 * MDP(13) + t57 * MDP(15) - t58 * MDP(16) + t118;
t97 = (MDP(22) * t94 - MDP(23) * t92) * pkin(4);
t111 = (t94 * t78 - t92 * t79) * MDP(22) - (t92 * t78 + t94 * t79) * MDP(23);
t115 = t78 * MDP(15) - t79 * MDP(16) + t111;
t110 = pkin(3) * t113 + t89 * qJ(2);
t106 = MDP(14) + MDP(21);
t105 = t88 * MDP(7) + t90 * MDP(8);
t102 = t74 * MDP(15) + t75 * MDP(16);
t99 = t61 * MDP(22) + t62 * MDP(23);
t96 = qJ(2) ^ 2;
t87 = t91 ^ 2;
t86 = t89 ^ 2;
t84 = t86 * t96;
t71 = -t88 * t109 + t77;
t1 = [MDP(1) + (pkin(1) ^ 2 + t84) * MDP(6) + (t71 ^ 2 + t72 ^ 2 + t84) * MDP(9) + (t96 * MDP(6) + t106) * t87 + (MDP(10) * t75 - 0.2e1 * t74 * MDP(11)) * t75 + (MDP(17) * t62 - 0.2e1 * t61 * MDP(18)) * t62 + 0.2e1 * t102 * t110 + 0.2e1 * t99 * (t74 * pkin(4) + t110) + 0.2e1 * (t87 * MDP(5) + (MDP(5) + t105) * t86) * qJ(2) + 0.2e1 * (pkin(1) * MDP(4) - t71 * MDP(7) + t72 * MDP(8) - t117) * t91; -pkin(1) * MDP(6) + (t71 * t90 + t72 * t88) * MDP(9) + (-t90 * MDP(7) + t88 * MDP(8) - MDP(4) - t115) * t91; MDP(6) + (t88 ^ 2 + t90 ^ 2) * MDP(9); (MDP(9) * qJ(2) + t105) * t89 + t99 + t102; 0; MDP(9); (-t106 - t97) * t91 + t117; t115; 0; t106 + 0.2e1 * t97; -t91 * MDP(21) + t118; t111; 0; MDP(21) + t97; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
