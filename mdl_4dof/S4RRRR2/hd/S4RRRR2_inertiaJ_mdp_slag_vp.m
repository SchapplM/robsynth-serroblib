% Calculate joint inertia matrix for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S4RRRR2_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:19
% EndTime: 2019-12-31 17:23:19
% DurationCPUTime: 0.13s
% Computational Cost: add. (133->48), mult. (238->61), div. (0->0), fcn. (207->6), ass. (0->30)
t64 = sin(qJ(2));
t55 = t64 * pkin(1) + pkin(6);
t63 = sin(qJ(3));
t46 = (-pkin(7) - t55) * t63;
t66 = cos(qJ(3));
t61 = t66 * pkin(7);
t47 = t66 * t55 + t61;
t62 = sin(qJ(4));
t65 = cos(qJ(4));
t85 = (t65 * t46 - t62 * t47) * MDP(19) + (-t62 * t46 - t65 * t47) * MDP(20);
t51 = (-pkin(6) - pkin(7)) * t63;
t52 = t66 * pkin(6) + t61;
t84 = (t65 * t51 - t62 * t52) * MDP(19) + (-t62 * t51 - t65 * t52) * MDP(20);
t72 = t66 * MDP(12) - t63 * MDP(13);
t48 = t62 * t63 - t65 * t66;
t49 = t62 * t66 + t65 * t63;
t83 = t48 * MDP(19) + t49 * MDP(20);
t67 = cos(qJ(2));
t82 = t67 * pkin(1);
t80 = t49 * MDP(16) - t48 * MDP(17);
t57 = -t66 * pkin(3) - pkin(2);
t74 = t66 * MDP(10) + t63 * MDP(9) + t80;
t73 = MDP(4) + (MDP(7) * t63 + 0.2e1 * MDP(8) * t66) * t63 + (MDP(14) * t49 - 0.2e1 * MDP(15) * t48) * t49;
t71 = -MDP(12) * t63 - MDP(13) * t66;
t70 = (t67 * MDP(5) - t64 * MDP(6)) * pkin(1);
t69 = (MDP(19) * t65 - MDP(20) * t62) * pkin(3);
t68 = 0.2e1 * t83;
t56 = -pkin(2) - t82;
t50 = t57 - t82;
t1 = [t50 * t68 - 0.2e1 * t56 * t72 + MDP(1) + 0.2e1 * t70 + t73; t70 + t73 + t72 * (pkin(2) - t56) + t83 * (t50 + t57); 0.2e1 * pkin(2) * t72 + t57 * t68 + t73; t55 * t71 + t74 + t85; pkin(6) * t71 + t74 + t84; MDP(11) + MDP(18) + 0.2e1 * t69; t80 + t85; t80 + t84; MDP(18) + t69; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
