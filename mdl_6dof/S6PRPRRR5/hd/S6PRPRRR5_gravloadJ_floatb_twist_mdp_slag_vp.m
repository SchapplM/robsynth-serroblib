% Calculate Gravitation load on the joints for
% S6PRPRRR5
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
%   see S6PRPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6PRPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:43:36
% EndTime: 2019-03-08 20:43:37
% DurationCPUTime: 0.31s
% Computational Cost: add. (249->71), mult. (457->126), div. (0->0), fcn. (534->12), ass. (0->35)
t68 = sin(pkin(6));
t91 = g(3) * t68;
t66 = qJ(4) + qJ(5);
t64 = sin(t66);
t71 = sin(qJ(6));
t90 = t64 * t71;
t74 = cos(qJ(6));
t89 = t64 * t74;
t67 = sin(pkin(11));
t88 = t67 * t68;
t69 = cos(pkin(11));
t87 = t68 * t69;
t72 = sin(qJ(4));
t86 = t68 * t72;
t75 = cos(qJ(4));
t85 = t68 * t75;
t76 = cos(qJ(2));
t84 = t68 * t76;
t70 = cos(pkin(6));
t73 = sin(qJ(2));
t83 = t70 * t73;
t82 = t70 * t76;
t81 = t71 * t73;
t80 = t73 * t74;
t60 = t67 * t82 + t69 * t73;
t65 = cos(t66);
t52 = t60 * t64 + t65 * t88;
t58 = t67 * t73 - t69 * t82;
t54 = -t58 * t64 + t65 * t87;
t57 = -t64 * t84 + t70 * t65;
t79 = (g(1) * t52 - g(2) * t54 + g(3) * t57) * MDP(21) + (-t74 * MDP(27) + t71 * MDP(28) - MDP(20)) * (g(1) * (t60 * t65 - t64 * t88) + g(2) * (t58 * t65 + t64 * t87) + g(3) * (-t70 * t64 - t65 * t84));
t50 = -g(1) * t60 - g(2) * t58 + g(3) * t84;
t61 = -t67 * t83 + t69 * t76;
t59 = t67 * t76 + t69 * t83;
t1 = [(-MDP(1) - MDP(7)) * g(3); (-g(1) * (-t60 * pkin(2) + t61 * qJ(3)) - g(2) * (-t58 * pkin(2) + t59 * qJ(3)) - (pkin(2) * t76 + qJ(3) * t73) * t91) * MDP(7) + (-g(1) * (-t60 * t71 + t61 * t89) - g(2) * (-t58 * t71 + t59 * t89) - (t64 * t80 + t71 * t76) * t91) * MDP(27) + (-g(1) * (-t60 * t74 - t61 * t90) - g(2) * (-t58 * t74 - t59 * t90) - (-t64 * t81 + t74 * t76) * t91) * MDP(28) + (-MDP(3) + MDP(5)) * t50 + (-t72 * MDP(13) - t75 * MDP(14) - MDP(20) * t64 - MDP(21) * t65 + MDP(4) - MDP(6)) * (g(1) * t61 + g(2) * t59 + t73 * t91); t50 * MDP(7); (-g(1) * (t60 * t75 - t67 * t86) - g(2) * (t58 * t75 + t69 * t86) - g(3) * (-t70 * t72 - t75 * t84)) * MDP(13) + (-g(1) * (-t60 * t72 - t67 * t85) - g(2) * (-t58 * t72 + t69 * t85) - g(3) * (-t70 * t75 + t72 * t84)) * MDP(14) + t79; t79; (-g(1) * (-t52 * t71 + t61 * t74) - g(2) * (t54 * t71 + t59 * t74) - g(3) * (-t57 * t71 + t68 * t80)) * MDP(27) + (-g(1) * (-t52 * t74 - t61 * t71) - g(2) * (t54 * t74 - t59 * t71) - g(3) * (-t57 * t74 - t68 * t81)) * MDP(28);];
taug  = t1;
