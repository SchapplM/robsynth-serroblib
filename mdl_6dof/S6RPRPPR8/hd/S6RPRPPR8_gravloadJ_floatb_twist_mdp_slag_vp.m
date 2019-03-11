% Calculate Gravitation load on the joints for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRPPR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:00:00
% EndTime: 2019-03-09 03:00:01
% DurationCPUTime: 0.35s
% Computational Cost: add. (117->68), mult. (246->90), div. (0->0), fcn. (207->6), ass. (0->34)
t92 = -MDP(12) - MDP(14) + MDP(19);
t91 = MDP(13) - MDP(16) - MDP(18);
t90 = -pkin(1) - pkin(7);
t89 = -pkin(3) - pkin(4);
t74 = cos(qJ(1));
t88 = g(2) * t74;
t70 = sin(qJ(3));
t87 = g(3) * t70;
t73 = cos(qJ(3));
t86 = g(3) * t73;
t71 = sin(qJ(1));
t85 = t70 * t71;
t84 = t70 * t74;
t83 = t71 * t73;
t69 = sin(qJ(6));
t82 = t74 * t69;
t72 = cos(qJ(6));
t81 = t74 * t72;
t78 = qJ(4) * t70;
t80 = pkin(3) * t83 + t71 * t78;
t79 = t74 * pkin(1) + t71 * qJ(2);
t63 = t73 * qJ(4);
t77 = MDP(17) + MDP(21);
t76 = pkin(3) * t85 + t74 * pkin(7) + t79;
t55 = g(1) * t74 + g(2) * t71;
t54 = g(1) * t71 - t88;
t64 = t74 * qJ(2);
t75 = pkin(3) * t84 - t74 * t63 + t64;
t51 = t71 * t69 - t73 * t81;
t50 = t71 * t72 + t73 * t82;
t49 = t72 * t83 + t82;
t48 = t69 * t83 - t81;
t47 = g(1) * t83 - t73 * t88 - t87;
t1 = [(-g(1) * (-t71 * pkin(1) + t64) - g(2) * t79) * MDP(6) + (-g(1) * (t90 * t71 + t75) - g(2) * (-t71 * t63 + t76)) * MDP(17) + (-g(1) * (pkin(4) * t84 + t75) - g(2) * (-t74 * qJ(5) + t76) + (-g(1) * (qJ(5) + t90) - g(2) * (t70 * pkin(4) - t63)) * t71) * MDP(21) + (-g(1) * t51 + g(2) * t49) * MDP(27) + (-g(1) * t50 - g(2) * t48) * MDP(28) + (MDP(2) - MDP(4) + MDP(15) - MDP(20)) * t54 + (t92 * t70 - t91 * t73 + MDP(3) - MDP(5)) * t55; (-MDP(6) - t77) * t54; (-g(1) * t80 - g(3) * (-t70 * pkin(3) + t63) - (-pkin(3) * t73 - t78) * t88) * MDP(17) + (-g(1) * (pkin(4) * t83 + t80) - g(3) * (t89 * t70 + t63) - (t89 * t73 - t78) * t88) * MDP(21) + t92 * t47 + t91 * (g(1) * t85 - g(2) * t84 + t86) + (-MDP(27) * t72 + MDP(28) * t69) * (t54 * t70 + t86); t77 * t47; t55 * MDP(21); (-g(1) * t48 + g(2) * t50 + t69 * t87) * MDP(27) + (-g(1) * t49 - g(2) * t51 + t72 * t87) * MDP(28);];
taug  = t1;
