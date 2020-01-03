% Calculate joint inertia matrix for
% S5RPPRR12
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
%   see S5RPPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPPRR12_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:22
% EndTime: 2019-12-31 18:07:23
% DurationCPUTime: 0.24s
% Computational Cost: add. (200->64), mult. (357->94), div. (0->0), fcn. (350->6), ass. (0->36)
t68 = sin(qJ(5));
t70 = cos(qJ(5));
t74 = t68 * MDP(23) + t70 * MDP(24);
t65 = sin(pkin(8));
t66 = cos(pkin(8));
t69 = sin(qJ(4));
t83 = cos(qJ(4));
t51 = t69 * t65 - t83 * t66;
t86 = t51 * MDP(17);
t85 = t51 ^ 2;
t52 = t83 * t65 + t69 * t66;
t48 = t52 ^ 2;
t67 = -pkin(1) - qJ(3);
t84 = -pkin(6) + t67;
t54 = t84 * t65;
t55 = t84 * t66;
t46 = t69 * t54 - t83 * t55;
t82 = t46 * t51;
t81 = t68 * t70;
t56 = t65 ^ 2 + t66 ^ 2;
t57 = t65 * pkin(3) + qJ(2);
t78 = MDP(19) * t81;
t77 = t65 * MDP(7) + t66 * MDP(8);
t76 = -MDP(20) * t70 + MDP(21) * t68;
t75 = t70 * MDP(23) - t68 * MDP(24);
t73 = MDP(16) + t75;
t72 = t68 * MDP(20) + t70 * MDP(21) - t74 * pkin(7);
t71 = qJ(2) ^ 2;
t64 = t70 ^ 2;
t63 = t68 ^ 2;
t50 = t56 * t67;
t47 = t83 * t54 + t69 * t55;
t45 = t52 * pkin(4) + t51 * pkin(7) + t57;
t44 = t68 * t45 + t70 * t47;
t43 = t70 * t45 - t68 * t47;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2) + t71) * MDP(6) + (t56 * t67 ^ 2 + t71) * MDP(10) + t48 * MDP(22) + (t64 * MDP(18) + MDP(11) - 0.2e1 * t78) * t85 - 0.2e1 * t50 * MDP(9) + 0.2e1 * (t43 * t52 - t68 * t82) * MDP(23) + 0.2e1 * (-t44 * t52 - t70 * t82) * MDP(24) + 0.2e1 * (t52 * MDP(16) - t86) * t57 + 0.2e1 * (MDP(5) + t77) * qJ(2) + 0.2e1 * (MDP(12) + t76) * t52 * t51; t50 * MDP(10) - pkin(1) * MDP(6) - t56 * MDP(9) + MDP(4) + t74 * (-t48 - t85); t56 * MDP(10) + MDP(6); qJ(2) * MDP(10) + t73 * t52 + t77 - t86; 0; MDP(10); -t47 * MDP(17) - t73 * t46 + (-MDP(14) + t72) * t52 + (-MDP(13) - MDP(18) * t81 + (t63 - t64) * MDP(19) + t74 * pkin(4)) * t51; -t52 * MDP(17) - t73 * t51; 0; t63 * MDP(18) + 0.2e1 * pkin(4) * t75 + MDP(15) + 0.2e1 * t78; t52 * MDP(22) + t43 * MDP(23) - t44 * MDP(24) + t76 * t51; -t74 * t52; t75; t72; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
