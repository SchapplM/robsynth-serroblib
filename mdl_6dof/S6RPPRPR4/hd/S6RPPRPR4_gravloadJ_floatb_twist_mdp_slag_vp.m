% Calculate Gravitation load on the joints for
% S6RPPRPR4
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
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:08
% EndTime: 2019-03-09 01:47:09
% DurationCPUTime: 0.27s
% Computational Cost: add. (137->54), mult. (239->78), div. (0->0), fcn. (253->10), ass. (0->31)
t61 = qJ(4) + pkin(10);
t54 = sin(t61);
t84 = g(3) * t54;
t82 = cos(qJ(1));
t81 = sin(qJ(1));
t55 = cos(t61);
t63 = sin(qJ(6));
t80 = t55 * t63;
t65 = cos(qJ(6));
t79 = t55 * t65;
t78 = t82 * pkin(1) + t81 * qJ(2);
t77 = cos(pkin(9));
t76 = sin(pkin(9));
t75 = MDP(18) + MDP(9);
t74 = t82 * pkin(2) + t78;
t73 = -t81 * pkin(1) + t82 * qJ(2);
t46 = -t81 * t76 - t82 * t77;
t47 = t82 * t76 - t81 * t77;
t72 = g(1) * t47 - g(2) * t46;
t71 = g(1) * t46 + g(2) * t47;
t70 = t46 * t63 + t47 * t79;
t69 = -t46 * t65 + t47 * t80;
t68 = -t81 * pkin(2) + t73;
t66 = cos(qJ(4));
t64 = sin(qJ(4));
t62 = -qJ(5) - pkin(7);
t53 = pkin(4) * t66 + pkin(3);
t48 = g(1) * t81 - g(2) * t82;
t44 = -t46 * t79 + t47 * t63;
t43 = t46 * t80 + t47 * t65;
t1 = [(-g(1) * t73 - g(2) * t78) * MDP(6) + (-g(1) * t68 - g(2) * t74) * MDP(9) + (-g(1) * (-t46 * t62 + t47 * t53 + t68) - g(2) * (-t46 * t53 - t47 * t62 + t74)) * MDP(18) + (-g(1) * t70 - g(2) * t44) * MDP(24) + (g(1) * t69 - g(2) * t43) * MDP(25) + (MDP(3) - MDP(5)) * (g(1) * t82 + g(2) * t81) + (MDP(2) + MDP(4)) * t48 + (MDP(8) - MDP(17)) * t71 + (-MDP(15) * t66 + MDP(16) * t64 - MDP(7)) * t72; (-MDP(6) - t75) * t48; t75 * g(3); (-g(3) * t64 - t71 * t66) * MDP(16) + (MDP(18) * pkin(4) + MDP(15)) * (g(3) * t66 - t71 * t64) + (-MDP(24) * t65 + MDP(25) * t63) * (-g(3) * t55 + t71 * t54); -t72 * MDP(18); (-g(1) * t43 - g(2) * t69 - t63 * t84) * MDP(24) + (g(1) * t44 - g(2) * t70 - t65 * t84) * MDP(25);];
taug  = t1;
