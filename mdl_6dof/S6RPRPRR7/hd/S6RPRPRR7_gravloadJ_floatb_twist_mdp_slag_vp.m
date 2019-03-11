% Calculate Gravitation load on the joints for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:04
% EndTime: 2019-03-09 03:56:05
% DurationCPUTime: 0.17s
% Computational Cost: add. (159->48), mult. (185->64), div. (0->0), fcn. (157->8), ass. (0->26)
t54 = sin(qJ(3));
t69 = pkin(3) * t54;
t47 = qJ(3) + pkin(10) + qJ(5);
t46 = cos(t47);
t67 = g(3) * t46;
t53 = sin(qJ(6));
t55 = sin(qJ(1));
t66 = t55 * t53;
t56 = cos(qJ(6));
t65 = t55 * t56;
t58 = cos(qJ(1));
t64 = t58 * t53;
t63 = t58 * t56;
t62 = t58 * pkin(1) + t55 * qJ(2);
t43 = g(1) * t55 - g(2) * t58;
t45 = sin(t47);
t61 = (t43 * t45 + t67) * MDP(22) + (-t56 * MDP(28) + t53 * MDP(29) - MDP(21)) * (-g(3) * t45 + t43 * t46);
t44 = g(1) * t58 + g(2) * t55;
t57 = cos(qJ(3));
t52 = -qJ(4) - pkin(7);
t49 = t58 * qJ(2);
t42 = t45 * t63 - t66;
t41 = t45 * t64 + t65;
t40 = t45 * t65 + t64;
t39 = -t45 * t66 + t63;
t1 = [(-g(1) * (-t55 * pkin(1) + t49) - g(2) * t62) * MDP(6) + (-g(1) * (t58 * t69 + t49 + (-pkin(1) + t52) * t55) - g(2) * (-t58 * t52 + t55 * t69 + t62)) * MDP(15) + (-g(1) * t42 - g(2) * t40) * MDP(28) + (g(1) * t41 - g(2) * t39) * MDP(29) + (MDP(2) - MDP(4) + MDP(14)) * t43 + (-t54 * MDP(12) - MDP(13) * t57 - MDP(21) * t45 - MDP(22) * t46 + MDP(3) - MDP(5)) * t44; (-MDP(15) - MDP(6)) * t43; (g(3) * t57 + t43 * t54) * MDP(13) + t61 + (MDP(15) * pkin(3) + MDP(12)) * (g(3) * t54 - t43 * t57); -t44 * MDP(15); t61; (-g(1) * t39 - g(2) * t41 + t53 * t67) * MDP(28) + (g(1) * t40 - g(2) * t42 + t56 * t67) * MDP(29);];
taug  = t1;
