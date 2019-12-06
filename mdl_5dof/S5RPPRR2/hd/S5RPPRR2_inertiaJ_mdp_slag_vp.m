% Calculate joint inertia matrix for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPPRR2_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:57
% EndTime: 2019-12-05 17:39:57
% DurationCPUTime: 0.13s
% Computational Cost: add. (191->55), mult. (311->79), div. (0->0), fcn. (324->6), ass. (0->31)
t66 = sin(pkin(8));
t67 = cos(pkin(8));
t70 = sin(qJ(4));
t72 = cos(qJ(4));
t56 = t72 * t66 + t70 * t67;
t61 = t66 * pkin(3) + qJ(2);
t84 = 0.2e1 * t56 * pkin(4) + 0.2e1 * t61;
t83 = 0.2e1 * t61;
t68 = -pkin(1) - qJ(3);
t82 = -pkin(6) + t68;
t57 = -t70 * t66 + t72 * t67;
t69 = sin(qJ(5));
t71 = cos(qJ(5));
t49 = t71 * t56 + t69 * t57;
t50 = -t69 * t56 + t71 * t57;
t81 = t50 * MDP(23) - t49 * MDP(24);
t60 = t66 ^ 2 + t67 ^ 2;
t80 = t49 * MDP(23);
t79 = t56 * MDP(16);
t58 = t82 * t66;
t59 = t82 * t67;
t78 = -t70 * t58 + t72 * t59;
t43 = -t57 * pkin(7) + t78;
t76 = -t72 * t58 - t70 * t59;
t44 = -t56 * pkin(7) - t76;
t77 = t50 * MDP(20) - t49 * MDP(21) + (t71 * t43 - t69 * t44) * MDP(23) + (-t69 * t43 - t71 * t44) * MDP(24);
t75 = t66 * MDP(7) + t67 * MDP(8);
t74 = (MDP(23) * t71 - MDP(24) * t69) * pkin(4);
t73 = qJ(2) ^ 2;
t55 = t60 * t68;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2) + t73) * MDP(6) - 0.2e1 * t55 * MDP(9) + (t60 * t68 ^ 2 + t73) * MDP(10) + t79 * t83 + t80 * t84 + (MDP(11) * t57 - 0.2e1 * t56 * MDP(12) + MDP(17) * t83) * t57 + (MDP(18) * t50 - 0.2e1 * t49 * MDP(19) + MDP(24) * t84) * t50 + 0.2e1 * (MDP(5) + t75) * qJ(2); t55 * MDP(10) - pkin(1) * MDP(6) - t60 * MDP(9) + MDP(4); t60 * MDP(10) + MDP(6); qJ(2) * MDP(10) + t57 * MDP(17) + t50 * MDP(24) + t75 + t79 + t80; 0; MDP(10); t57 * MDP(13) - t56 * MDP(14) + MDP(16) * t78 + MDP(17) * t76 + t77; t57 * MDP(16) - t56 * MDP(17) + t81; 0; MDP(15) + MDP(22) + 0.2e1 * t74; t77; t81; 0; MDP(22) + t74; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
