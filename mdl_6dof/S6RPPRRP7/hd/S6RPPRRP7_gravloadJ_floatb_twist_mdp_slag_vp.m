% Calculate Gravitation load on the joints for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:13:54
% EndTime: 2019-03-09 02:13:55
% DurationCPUTime: 0.34s
% Computational Cost: add. (154->59), mult. (213->74), div. (0->0), fcn. (179->8), ass. (0->31)
t67 = sin(qJ(1));
t69 = cos(qJ(1));
t90 = -g(1) * t67 + g(2) * t69;
t89 = MDP(17) - MDP(25);
t61 = pkin(9) + qJ(4);
t55 = sin(t61);
t56 = cos(t61);
t46 = -g(3) * t55 - t56 * t90;
t82 = g(3) * t56;
t66 = sin(qJ(5));
t81 = t67 * t66;
t68 = cos(qJ(5));
t80 = t67 * t68;
t79 = t69 * t66;
t78 = t69 * t68;
t77 = t69 * pkin(1) + t67 * qJ(2);
t76 = -MDP(10) - MDP(26);
t74 = g(2) * t77;
t73 = pkin(5) * t66 + pkin(7) + qJ(3);
t53 = g(1) * t69 + g(2) * t67;
t54 = t68 * pkin(5) + pkin(4);
t64 = -qJ(6) - pkin(8);
t71 = -t55 * t54 - t56 * t64;
t50 = t55 * t79 + t80;
t48 = -t55 * t81 + t78;
t62 = sin(pkin(9));
t70 = pkin(3) * t62 - t71;
t58 = t69 * qJ(2);
t51 = t55 * t78 - t81;
t49 = t55 * t80 + t79;
t1 = [(-g(1) * (-t67 * pkin(1) + t58) - t74) * MDP(6) + (-g(1) * (t58 + (-pkin(1) - qJ(3)) * t67) - g(2) * (t69 * qJ(3) + t77)) * MDP(10) + (-g(1) * t51 - g(2) * t49) * MDP(23) + (g(1) * t50 - g(2) * t48) * MDP(24) + (-g(1) * t58 - t74 + (-g(1) * t70 - g(2) * t73) * t69 + (-g(1) * (-pkin(1) - t73) - g(2) * t70) * t67) * MDP(26) - (MDP(2) - MDP(4) + MDP(9)) * t90 + (-t55 * MDP(16) - t62 * MDP(7) - MDP(8) * cos(pkin(9)) - t89 * t56 + MDP(3) - MDP(5)) * t53; -(-MDP(6) + t76) * t90; t76 * t53; (-g(3) * t71 + t90 * (t54 * t56 - t55 * t64)) * MDP(26) + t89 * (-t55 * t90 + t82) + (-MDP(23) * t68 + MDP(24) * t66 - MDP(16)) * t46; (g(1) * t49 - g(2) * t51 + t68 * t82) * MDP(24) + (pkin(5) * MDP(26) + MDP(23)) * (-g(1) * t48 - g(2) * t50 + t66 * t82); t46 * MDP(26);];
taug  = t1;
