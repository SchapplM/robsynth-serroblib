% Calculate joint inertia matrix for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRRP4_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:21
% EndTime: 2019-12-05 16:46:22
% DurationCPUTime: 0.23s
% Computational Cost: add. (187->67), mult. (349->91), div. (0->0), fcn. (278->6), ass. (0->38)
t71 = sin(qJ(3));
t60 = t71 * pkin(2) + pkin(7);
t70 = sin(qJ(4));
t68 = t70 ^ 2;
t73 = cos(qJ(4));
t93 = t73 ^ 2 + t68;
t95 = t93 * t60;
t99 = 2 * MDP(16);
t74 = cos(qJ(3));
t98 = t74 * pkin(2);
t61 = -pkin(3) - t98;
t97 = pkin(3) - t61;
t55 = -t73 * pkin(4) - t70 * qJ(5) - pkin(3);
t48 = t55 - t98;
t96 = -t48 - t55;
t94 = pkin(7) * t93;
t92 = MDP(18) * t70;
t91 = t48 * MDP(18);
t90 = t55 * MDP(18);
t89 = t70 * MDP(14);
t88 = t70 * MDP(17);
t87 = t70 * MDP(10) + t73 * MDP(11) + (-t70 * pkin(4) + t73 * qJ(5)) * MDP(16);
t86 = 0.2e1 * t70 * t73 * MDP(9) + t68 * MDP(8) + MDP(5);
t85 = t93 * MDP(18);
t72 = sin(qJ(2));
t75 = cos(qJ(2));
t52 = t71 * t72 - t74 * t75;
t53 = t71 * t75 + t74 * t72;
t84 = (MDP(16) * t93 - MDP(7)) * t53 + (-MDP(6) + t89) * t52;
t83 = -MDP(18) * pkin(4) - MDP(15);
t82 = t53 * t85;
t81 = t73 * MDP(13) - t89;
t80 = -0.2e1 * t73 * MDP(15) - 0.2e1 * t88;
t79 = (t74 * MDP(6) - t71 * MDP(7)) * pkin(2);
t78 = (-MDP(13) - MDP(15)) * t73 - t88;
t77 = (MDP(18) * qJ(5) - MDP(14) + MDP(17)) * t73 + (-MDP(13) + t83) * t70;
t62 = t70 * MDP(16);
t1 = [MDP(1) + (t93 * t53 ^ 2 + t52 ^ 2) * MDP(18); t75 * MDP(3) - t72 * MDP(4) + t60 * t82 + (t78 + t91) * t52 + t84; MDP(2) + t95 * t99 + t60 ^ 2 * t85 + (t80 + t91) * t48 + t86 - 0.2e1 * t81 * t61 + 0.2e1 * t79; pkin(7) * t82 + (t78 + t90) * t52 + t84; (t94 + t95) * MDP(16) + (pkin(7) * t95 + t48 * t55) * MDP(18) + t79 + (t97 * MDP(13) + t96 * MDP(15)) * t73 + (-t97 * MDP(14) + t96 * MDP(17)) * t70 + t86; t94 * t99 + pkin(7) ^ 2 * t85 + (t80 + t90) * t55 + 0.2e1 * t81 * pkin(3) + t86; t77 * t53; t77 * t60 + t87; t77 * pkin(7) + t87; MDP(12) + 0.2e1 * pkin(4) * MDP(15) + 0.2e1 * qJ(5) * MDP(17) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(18); t53 * t92; t60 * t92 + t62; pkin(7) * t92 + t62; t83; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
