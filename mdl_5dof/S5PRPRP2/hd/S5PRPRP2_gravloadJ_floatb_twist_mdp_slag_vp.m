% Calculate Gravitation load on the joints for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRPRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:52
% EndTime: 2021-01-15 15:04:53
% DurationCPUTime: 0.17s
% Computational Cost: add. (140->40), mult. (158->58), div. (0->0), fcn. (145->6), ass. (0->24)
t53 = sin(pkin(8));
t54 = cos(pkin(8));
t57 = cos(qJ(4));
t72 = -t53 * (-qJ(5) - pkin(6)) + (t57 * pkin(4) + pkin(3)) * t54;
t71 = MDP(13) + MDP(15);
t70 = MDP(14) + MDP(16);
t66 = g(3) * t53;
t52 = pkin(7) + qJ(2);
t51 = cos(t52);
t56 = sin(qJ(4));
t64 = t51 * t56;
t62 = t54 * t56;
t61 = t54 * t57;
t50 = sin(t52);
t60 = t51 * pkin(2) + t50 * qJ(3);
t59 = -MDP(18) - MDP(7);
t44 = g(1) * t51 + g(2) * t50;
t43 = g(1) * t50 - g(2) * t51;
t41 = t50 * t57 - t51 * t62;
t39 = t50 * t62 + t51 * t57;
t46 = t51 * qJ(3);
t42 = t50 * t56 + t51 * t61;
t40 = -t50 * t61 + t64;
t1 = [(-MDP(1) + t59) * g(3); (-g(1) * (-t50 * pkin(2) + t46) - g(2) * t60) * MDP(7) + (-g(1) * (pkin(4) * t64 + t46) - g(2) * (t72 * t51 + t60) + (-g(1) * (-pkin(2) - t72) - g(2) * pkin(4) * t56) * t50) * MDP(18) + (MDP(4) - MDP(6)) * t44 + t71 * (-g(1) * t40 - g(2) * t42) + t70 * (-g(1) * t39 - g(2) * t41) + (t53 * MDP(17) + t54 * MDP(5) + MDP(3)) * t43; t59 * t43; t70 * (g(1) * t42 - g(2) * t40 + t57 * t66) + (pkin(4) * MDP(18) + t71) * (-g(1) * t41 + g(2) * t39 + t56 * t66); (g(3) * t54 - t44 * t53) * MDP(18);];
taug = t1;
