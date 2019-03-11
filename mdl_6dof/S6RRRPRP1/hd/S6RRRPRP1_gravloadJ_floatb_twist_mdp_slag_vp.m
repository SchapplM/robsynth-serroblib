% Calculate Gravitation load on the joints for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:33:01
% EndTime: 2019-03-09 16:33:03
% DurationCPUTime: 0.31s
% Computational Cost: add. (288->69), mult. (290->92), div. (0->0), fcn. (245->10), ass. (0->41)
t82 = qJ(2) + qJ(3);
t76 = pkin(10) + t82;
t72 = sin(t76);
t73 = cos(t76);
t87 = cos(qJ(5));
t75 = t87 * pkin(5) + pkin(4);
t83 = -qJ(6) - pkin(9);
t95 = -t72 * t83 + t73 * t75;
t86 = sin(qJ(1));
t89 = cos(qJ(1));
t71 = g(1) * t89 + g(2) * t86;
t113 = -g(3) * t73 + t71 * t72;
t77 = sin(t82);
t112 = pkin(3) * t77;
t107 = g(3) * t72;
t84 = sin(qJ(5));
t104 = t86 * t84;
t103 = t86 * t87;
t102 = t89 * t84;
t101 = t89 * t87;
t78 = cos(t82);
t74 = pkin(3) * t78;
t88 = cos(qJ(2));
t79 = t88 * pkin(2);
t100 = t74 + t79;
t81 = -qJ(4) - pkin(8) - pkin(7);
t98 = pkin(5) * t84 - t81;
t97 = t74 + t95;
t90 = -g(3) * t78 + t71 * t77;
t96 = (-t71 * t73 - t107) * MDP(27) + t90 * MDP(16) + (g(3) * t77 + t71 * t78) * MDP(17) + (t87 * MDP(25) - t84 * MDP(26)) * t113;
t70 = g(1) * t86 - g(2) * t89;
t94 = -t72 * t75 - t73 * t83;
t64 = -t73 * t102 + t103;
t62 = t73 * t104 + t101;
t85 = sin(qJ(2));
t69 = -t85 * pkin(2) - t112;
t68 = pkin(1) + t100;
t66 = t89 * t68;
t65 = t73 * t101 + t104;
t63 = -t73 * t103 + t102;
t1 = [(-g(1) * (-t86 * t68 - t89 * t81) - g(2) * (-t86 * t81 + t66)) * MDP(19) + (-g(1) * t63 - g(2) * t65) * MDP(25) + (-g(1) * t62 - g(2) * t64) * MDP(26) + (-g(2) * t66 + (-g(1) * t98 - g(2) * t95) * t89 + (-g(1) * (-t68 - t95) - g(2) * t98) * t86) * MDP(28) + (MDP(3) - MDP(18)) * t71 + (-t85 * MDP(10) + MDP(16) * t78 - MDP(17) * t77 + t72 * MDP(27) + t88 * MDP(9) + MDP(2)) * t70; (-g(3) * t88 + t71 * t85) * MDP(9) + (g(3) * t85 + t71 * t88) * MDP(10) + (-g(3) * t100 - t71 * t69) * MDP(19) + (-g(3) * (t79 + t97) + t71 * (-t69 - t94)) * MDP(28) + t96; t90 * pkin(3) * MDP(19) + (-g(3) * t97 + t71 * (-t94 + t112)) * MDP(28) + t96; (-MDP(19) - MDP(28)) * t70; (g(1) * t65 - g(2) * t63 + t87 * t107) * MDP(26) + (pkin(5) * MDP(28) + MDP(25)) * (-g(1) * t64 + g(2) * t62 + t84 * t107); -t113 * MDP(28);];
taug  = t1;
