% Calculate Gravitation load on the joints for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRPRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:19
% EndTime: 2019-12-05 15:41:21
% DurationCPUTime: 0.24s
% Computational Cost: add. (87->44), mult. (209->64), div. (0->0), fcn. (192->6), ass. (0->26)
t73 = MDP(13) + MDP(15);
t72 = -MDP(14) + MDP(17);
t51 = sin(pkin(7));
t52 = cos(pkin(7));
t71 = -g(1) * t52 - g(2) * t51;
t54 = sin(qJ(2));
t70 = pkin(2) * t54;
t56 = cos(qJ(2));
t66 = g(3) * t56;
t53 = sin(qJ(4));
t65 = t53 * t54;
t55 = cos(qJ(4));
t64 = t54 * t55;
t63 = t56 * pkin(2) + t54 * qJ(3);
t62 = qJ(3) * t56;
t61 = -MDP(18) - MDP(7);
t59 = pkin(4) * t53 - qJ(5) * t55;
t41 = t51 * t53 - t52 * t64;
t43 = t51 * t64 + t52 * t53;
t36 = g(1) * t41 - g(2) * t43 + t55 * t66;
t47 = t52 * t62;
t46 = t51 * t62;
t44 = -t51 * t65 + t52 * t55;
t42 = t51 * t55 + t52 * t65;
t39 = -t71 * t54 - t66;
t1 = [(-MDP(1) + t61) * g(3); (-g(1) * (-t52 * t70 + t47) - g(2) * (-t51 * t70 + t46) - g(3) * t63) * MDP(7) + (-g(1) * t47 - g(2) * t46 - g(3) * (t56 * pkin(6) + t59 * t54 + t63) + t71 * (t59 * t56 + (-pkin(2) - pkin(6)) * t54)) * MDP(18) + (MDP(3) - MDP(5) + MDP(16)) * t39 + (-t73 * t53 + t72 * t55 + MDP(4) - MDP(6)) * (g(3) * t54 - t71 * t56); t61 * t39; (-g(1) * (-t41 * pkin(4) + t42 * qJ(5)) - g(2) * (t43 * pkin(4) - t44 * qJ(5)) - (-pkin(4) * t55 - qJ(5) * t53) * t66) * MDP(18) + t72 * (-g(1) * t42 + g(2) * t44 + t53 * t66) + t73 * t36; -t36 * MDP(18);];
taug = t1;
