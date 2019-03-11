% Calculate Gravitation load on the joints for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:49:26
% EndTime: 2019-03-09 01:49:27
% DurationCPUTime: 0.31s
% Computational Cost: add. (124->61), mult. (215->84), div. (0->0), fcn. (186->8), ass. (0->33)
t63 = sin(qJ(4));
t65 = cos(qJ(4));
t68 = -t63 * pkin(4) + t65 * qJ(5);
t64 = sin(qJ(1));
t66 = cos(qJ(1));
t52 = g(1) * t66 + g(2) * t64;
t88 = MDP(16) - MDP(19);
t49 = -g(3) * t63 + t52 * t65;
t84 = g(3) * t65;
t60 = pkin(9) + qJ(6);
t53 = sin(t60);
t82 = t64 * t53;
t54 = cos(t60);
t81 = t64 * t54;
t61 = sin(pkin(9));
t80 = t64 * t61;
t62 = cos(pkin(9));
t79 = t64 * t62;
t78 = t66 * t53;
t77 = t66 * t54;
t76 = t66 * t61;
t75 = t66 * t62;
t74 = -pkin(1) - qJ(3);
t73 = t66 * pkin(1) + t64 * qJ(2);
t71 = -MDP(20) - MDP(9);
t70 = t66 * qJ(3) + t73;
t51 = g(1) * t64 - g(2) * t66;
t57 = t66 * qJ(2);
t47 = t63 * t77 - t82;
t46 = -t63 * t78 - t81;
t45 = -t63 * t81 - t78;
t44 = t63 * t82 - t77;
t1 = [(-g(1) * (-t64 * pkin(1) + t57) - g(2) * t73) * MDP(6) + (-g(1) * (t64 * t74 + t57) - g(2) * t70) * MDP(9) + (-g(1) * (-t63 * t79 - t76) - g(2) * (t63 * t75 - t80)) * MDP(17) + (-g(1) * (t63 * t80 - t75) - g(2) * (-t63 * t76 - t79)) * MDP(18) + (-g(1) * (-t66 * pkin(7) + t57) - g(2) * (-t68 * t66 + t70) + (-g(1) * (t68 + t74) + g(2) * pkin(7)) * t64) * MDP(20) + (-g(1) * t45 - g(2) * t47) * MDP(26) + (-g(1) * t44 - g(2) * t46) * MDP(27) + (MDP(3) - MDP(5) - MDP(7)) * t52 + (MDP(15) * t63 + t88 * t65 + MDP(2) - MDP(4) + MDP(8)) * t51; (-MDP(6) + t71) * t51; t71 * t52; (-g(3) * t68 - t52 * (pkin(4) * t65 + qJ(5) * t63)) * MDP(20) + t88 * (t52 * t63 + t84) + (-MDP(17) * t62 + MDP(18) * t61 - MDP(26) * t54 + MDP(27) * t53 - MDP(15)) * t49; t49 * MDP(20); (-g(1) * t46 + g(2) * t44 + t53 * t84) * MDP(26) + (g(1) * t47 - g(2) * t45 + t54 * t84) * MDP(27);];
taug  = t1;
