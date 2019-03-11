% Calculate Gravitation load on the joints for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:33
% EndTime: 2019-03-09 01:51:34
% DurationCPUTime: 0.24s
% Computational Cost: add. (93->52), mult. (193->70), div. (0->0), fcn. (161->6), ass. (0->25)
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t59 = -pkin(4) * t53 + qJ(5) * t56;
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t46 = g(1) * t57 + g(2) * t54;
t75 = MDP(15) - MDP(18);
t74 = MDP(16) - MDP(19);
t69 = g(3) * t53;
t67 = t54 * t56;
t66 = t56 * t57;
t65 = -pkin(1) - qJ(3);
t64 = t57 * pkin(1) + t54 * qJ(2);
t62 = -MDP(20) - MDP(9);
t61 = t57 * qJ(3) + t64;
t45 = g(1) * t54 - g(2) * t57;
t55 = cos(qJ(6));
t52 = sin(qJ(6));
t49 = t57 * qJ(2);
t42 = t52 * t67 - t55 * t57;
t41 = t52 * t57 + t55 * t67;
t40 = t52 * t66 + t54 * t55;
t39 = t54 * t52 - t55 * t66;
t38 = t46 * t56 - t69;
t1 = [(-g(1) * (-t54 * pkin(1) + t49) - g(2) * t64) * MDP(6) + (-g(1) * (t65 * t54 + t49) - g(2) * t61) * MDP(9) + (-g(1) * (-pkin(7) * t57 + t49) - g(2) * (-t59 * t57 + t61) + (-g(1) * (t59 + t65) + g(2) * pkin(7)) * t54) * MDP(20) + (-g(1) * t42 + g(2) * t40) * MDP(26) + (-g(1) * t41 - g(2) * t39) * MDP(27) + (MDP(3) - MDP(5) - MDP(7) + MDP(17)) * t46 + (t75 * t53 + t74 * t56 + MDP(2) - MDP(4) + MDP(8)) * t45; (-MDP(6) + t62) * t45; t62 * t46; (-g(3) * t59 - t46 * (pkin(4) * t56 + qJ(5) * t53)) * MDP(20) - t75 * t38 + (-MDP(26) * t52 - MDP(27) * t55 + t74) * (g(3) * t56 + t46 * t53); t38 * MDP(20); (-g(1) * t39 + g(2) * t41 - t55 * t69) * MDP(26) + (-g(1) * t40 - g(2) * t42 + t52 * t69) * MDP(27);];
taug  = t1;
