% Calculate Gravitation load on the joints for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:44
% EndTime: 2019-03-09 01:53:45
% DurationCPUTime: 0.37s
% Computational Cost: add. (182->66), mult. (224->87), div. (0->0), fcn. (194->10), ass. (0->36)
t70 = sin(qJ(1));
t71 = cos(qJ(1));
t92 = -g(1) * t70 + g(2) * t71;
t91 = MDP(17) - MDP(20);
t64 = pkin(9) + qJ(4);
t56 = sin(t64);
t58 = cos(t64);
t47 = -g(3) * t56 - t58 * t92;
t87 = g(3) * t58;
t63 = pkin(10) + qJ(6);
t55 = sin(t63);
t86 = t70 * t55;
t57 = cos(t63);
t85 = t70 * t57;
t65 = sin(pkin(10));
t84 = t70 * t65;
t67 = cos(pkin(10));
t83 = t70 * t67;
t82 = t71 * t55;
t81 = t71 * t57;
t80 = t71 * t65;
t79 = t71 * t67;
t78 = t71 * pkin(1) + t70 * qJ(2);
t77 = -MDP(10) - MDP(21);
t76 = g(2) * t78;
t54 = g(1) * t71 + g(2) * t70;
t74 = -t56 * pkin(4) + t58 * qJ(5);
t66 = sin(pkin(9));
t73 = pkin(3) * t66 - t74;
t69 = -pkin(7) - qJ(3);
t60 = t71 * qJ(2);
t51 = t56 * t81 - t86;
t50 = t56 * t82 + t85;
t49 = t56 * t85 + t82;
t48 = -t56 * t86 + t81;
t1 = [(-g(1) * (-t70 * pkin(1) + t60) - t76) * MDP(6) + (-g(1) * (t60 + (-pkin(1) - qJ(3)) * t70) - g(2) * (t71 * qJ(3) + t78)) * MDP(10) + (-g(1) * (t56 * t79 - t84) - g(2) * (t56 * t83 + t80)) * MDP(18) + (-g(1) * (-t56 * t80 - t83) - g(2) * (-t56 * t84 + t79)) * MDP(19) + (-g(1) * t60 - t76 + (-g(1) * t73 + g(2) * t69) * t71 + (-g(1) * (-pkin(1) + t69) - g(2) * t73) * t70) * MDP(21) + (-g(1) * t51 - g(2) * t49) * MDP(27) + (g(1) * t50 - g(2) * t48) * MDP(28) - (MDP(2) - MDP(4) + MDP(9)) * t92 + (-t56 * MDP(16) - t66 * MDP(7) - MDP(8) * cos(pkin(9)) - t91 * t58 + MDP(3) - MDP(5)) * t54; -(-MDP(6) + t77) * t92; t77 * t54; (-g(3) * t74 + t92 * (pkin(4) * t58 + qJ(5) * t56)) * MDP(21) + t91 * (-t56 * t92 + t87) + (-MDP(18) * t67 + MDP(19) * t65 - MDP(27) * t57 + MDP(28) * t55 - MDP(16)) * t47; t47 * MDP(21); (-g(1) * t48 - g(2) * t50 + t55 * t87) * MDP(27) + (g(1) * t49 - g(2) * t51 + t57 * t87) * MDP(28);];
taug  = t1;
