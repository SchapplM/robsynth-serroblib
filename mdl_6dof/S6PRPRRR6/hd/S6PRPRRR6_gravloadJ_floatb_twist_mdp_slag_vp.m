% Calculate Gravitation load on the joints for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6PRPRRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:48:47
% EndTime: 2019-03-08 20:48:48
% DurationCPUTime: 0.34s
% Computational Cost: add. (249->83), mult. (535->147), div. (0->0), fcn. (642->12), ass. (0->36)
t65 = sin(pkin(6));
t89 = g(3) * t65;
t63 = qJ(5) + qJ(6);
t61 = sin(t63);
t67 = sin(qJ(4));
t88 = t61 * t67;
t62 = cos(t63);
t87 = t62 * t67;
t64 = sin(pkin(11));
t86 = t64 * t65;
t68 = sin(qJ(2));
t85 = t65 * t68;
t71 = cos(qJ(2));
t84 = t65 * t71;
t66 = sin(qJ(5));
t83 = t66 * t67;
t82 = t67 * t68;
t69 = cos(qJ(5));
t81 = t67 * t69;
t80 = t68 * t69;
t78 = cos(pkin(6));
t76 = t64 * t78;
t77 = cos(pkin(11));
t55 = t68 * t77 + t71 * t76;
t70 = cos(qJ(4));
t49 = t55 * t67 + t70 * t86;
t74 = t78 * t77;
t53 = t64 * t68 - t71 * t74;
t75 = t65 * t77;
t51 = -t53 * t67 + t70 * t75;
t54 = t64 * t71 + t68 * t74;
t56 = -t68 * t76 + t71 * t77;
t58 = -t67 * t84 + t70 * t78;
t79 = (-g(1) * (-t49 * t61 + t56 * t62) - g(2) * (t51 * t61 + t54 * t62) - g(3) * (-t58 * t61 + t62 * t85)) * MDP(27) + (-g(1) * (-t49 * t62 - t56 * t61) - g(2) * (t51 * t62 - t54 * t61) - g(3) * (-t58 * t62 - t61 * t85)) * MDP(28);
t47 = -g(1) * t55 - g(2) * t53 + g(3) * t84;
t1 = [(-MDP(1) - MDP(7)) * g(3); (-g(1) * (-t55 * pkin(2) + t56 * qJ(3)) - g(2) * (-t53 * pkin(2) + t54 * qJ(3)) - (pkin(2) * t71 + qJ(3) * t68) * t89) * MDP(7) + (-g(1) * (-t55 * t66 + t56 * t81) - g(2) * (-t53 * t66 + t54 * t81) - (t66 * t71 + t67 * t80) * t89) * MDP(20) + (-g(1) * (-t55 * t69 - t56 * t83) - g(2) * (-t53 * t69 - t54 * t83) - (-t66 * t82 + t69 * t71) * t89) * MDP(21) + (-g(1) * (-t55 * t61 + t56 * t87) - g(2) * (-t53 * t61 + t54 * t87) - (t61 * t71 + t62 * t82) * t89) * MDP(27) + (-g(1) * (-t55 * t62 - t56 * t88) - g(2) * (-t53 * t62 - t54 * t88) - (-t61 * t82 + t62 * t71) * t89) * MDP(28) + (-MDP(3) + MDP(5)) * t47 + (-t67 * MDP(13) - t70 * MDP(14) + MDP(4) - MDP(6)) * (g(1) * t56 + g(2) * t54 + g(3) * t85); t47 * MDP(7); (g(1) * t49 - g(2) * t51 + g(3) * t58) * MDP(14) + (-MDP(20) * t69 + MDP(21) * t66 - MDP(27) * t62 + MDP(28) * t61 - MDP(13)) * (g(1) * (t55 * t70 - t67 * t86) + g(2) * (t53 * t70 + t67 * t75) + g(3) * (-t67 * t78 - t70 * t84)); (-g(1) * (-t49 * t66 + t56 * t69) - g(2) * (t51 * t66 + t54 * t69) - g(3) * (-t58 * t66 + t65 * t80)) * MDP(20) + (-g(1) * (-t49 * t69 - t56 * t66) - g(2) * (t51 * t69 - t54 * t66) - g(3) * (-t58 * t69 - t66 * t85)) * MDP(21) + t79; t79;];
taug  = t1;
