% Calculate joint inertia matrix for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRPR1_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:43
% EndTime: 2019-12-05 17:47:44
% DurationCPUTime: 0.12s
% Computational Cost: add. (230->63), mult. (373->96), div. (0->0), fcn. (389->6), ass. (0->32)
t78 = sin(pkin(8));
t79 = cos(pkin(8));
t81 = sin(qJ(3));
t83 = cos(qJ(3));
t68 = -t78 * t83 - t79 * t81;
t76 = t81 * pkin(3) + qJ(2);
t97 = -0.2e1 * t68 * pkin(4) + 0.2e1 * t76;
t96 = pkin(3) * t78;
t95 = (pkin(1) * MDP(6));
t67 = -t78 * t81 + t79 * t83;
t80 = sin(qJ(5));
t82 = cos(qJ(5));
t55 = t80 * t67 - t82 * t68;
t56 = t82 * t67 + t80 * t68;
t94 = t56 * MDP(21) - t55 * MDP(22);
t84 = -pkin(1) - pkin(6);
t93 = -qJ(4) + t84;
t70 = t93 * t81;
t71 = t93 * t83;
t58 = t79 * t70 + t78 * t71;
t92 = t55 * MDP(21);
t75 = t79 * pkin(3) + pkin(4);
t91 = (t82 * t75 - t80 * t96) * MDP(21);
t90 = (-t80 * t75 - t82 * t96) * MDP(22);
t89 = t67 ^ 2 + t68 ^ 2;
t57 = -t78 * t70 + t79 * t71;
t49 = -t67 * pkin(7) + t57;
t50 = t68 * pkin(7) + t58;
t88 = t56 * MDP(18) - t55 * MDP(19) + (t82 * t49 - t80 * t50) * MDP(21) + (-t80 * t49 - t82 * t50) * MDP(22);
t87 = t57 * t67 - t58 * t68;
t86 = t67 * t79 - t68 * t78;
t1 = [MDP(1) - 0.2e1 * t87 * MDP(14) + (t57 ^ 2 + t58 ^ 2 + t76 ^ 2) * MDP(15) + t92 * t97 + (MDP(7) * t83 - 0.2e1 * t81 * MDP(8)) * t83 + ((-2 * MDP(4) + t95) * pkin(1)) + (MDP(16) * t56 - 0.2e1 * t55 * MDP(17) + MDP(22) * t97) * t56 + (0.2e1 * t81 * MDP(12) + 0.2e1 * t83 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2); -t89 * MDP(14) + t87 * MDP(15) + MDP(4) - t95; t89 * MDP(15) + MDP(6); (t84 * MDP(12) + MDP(9)) * t83 + (-t84 * MDP(13) - MDP(10)) * t81 + (-t86 * MDP(14) + (t57 * t79 + t58 * t78) * MDP(15)) * pkin(3) + t88; t86 * MDP(15) * pkin(3) + t83 * MDP(12) - t81 * MDP(13) + t94; MDP(11) + MDP(20) + (t78 ^ 2 + t79 ^ 2) * MDP(15) * pkin(3) ^ 2 + 0.2e1 * t91 + 0.2e1 * t90; t76 * MDP(15) + t56 * MDP(22) + t92; 0; 0; MDP(15); t88; t94; MDP(20) + t90 + t91; 0; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
