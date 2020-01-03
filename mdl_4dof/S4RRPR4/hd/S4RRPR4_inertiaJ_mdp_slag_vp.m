% Calculate joint inertia matrix for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RRPR4_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:34
% DurationCPUTime: 0.14s
% Computational Cost: add. (126->49), mult. (231->64), div. (0->0), fcn. (185->6), ass. (0->34)
t64 = sin(pkin(7));
t65 = cos(pkin(7));
t80 = t64 ^ 2 + t65 ^ 2;
t81 = t80 * qJ(3);
t66 = sin(qJ(4));
t68 = cos(qJ(4));
t46 = t66 * t64 - t68 * t65;
t47 = t68 * t64 + t66 * t65;
t89 = t46 * MDP(16) + t47 * MDP(17);
t88 = t65 * MDP(7) - t64 * MDP(8);
t87 = 2 * MDP(9);
t69 = cos(qJ(2));
t86 = t69 * pkin(1);
t84 = t47 * MDP(13) - t46 * MDP(14);
t67 = sin(qJ(2));
t55 = t67 * pkin(1) + qJ(3);
t82 = t80 * t55;
t79 = pkin(2) * MDP(10);
t57 = -pkin(2) - t86;
t77 = t57 * MDP(10);
t76 = MDP(4) + (MDP(11) * t47 - 0.2e1 * MDP(12) * t46) * t47;
t56 = -t65 * pkin(3) - pkin(2);
t75 = t80 * MDP(10);
t74 = 0.2e1 * t88;
t73 = -t88 + t89;
t72 = (t69 * MDP(5) - t67 * MDP(6)) * pkin(1);
t71 = 0.2e1 * t89;
t61 = t65 * pkin(6);
t50 = t65 * qJ(3) + t61;
t49 = (-pkin(6) - qJ(3)) * t64;
t48 = t56 - t86;
t45 = t65 * t55 + t61;
t44 = (-pkin(6) - t55) * t64;
t1 = [MDP(1) + t82 * t87 + t55 ^ 2 * t75 + (-t74 + t77) * t57 + t48 * t71 + 0.2e1 * t72 + t76; (t81 + t82) * MDP(9) + (-t57 * pkin(2) + t55 * t81) * MDP(10) + t72 + t76 + t88 * (pkin(2) - t57) + t89 * (t48 + t56); t81 * t87 + qJ(3) ^ 2 * t75 + t56 * t71 + (t74 + t79) * pkin(2) + t76; t73 + t77; t73 - t79; MDP(10); (t68 * t44 - t66 * t45) * MDP(16) + (-t66 * t44 - t68 * t45) * MDP(17) + t84; (t68 * t49 - t66 * t50) * MDP(16) + (-t66 * t49 - t68 * t50) * MDP(17) + t84; 0; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
