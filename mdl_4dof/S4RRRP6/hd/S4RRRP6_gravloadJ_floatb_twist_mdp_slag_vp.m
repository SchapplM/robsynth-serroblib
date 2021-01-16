% Calculate Gravitation load on the joints for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:39
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:39:15
% EndTime: 2021-01-15 14:39:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (97->40), mult. (207->57), div. (0->0), fcn. (191->6), ass. (0->25)
t76 = MDP(10) - MDP(20);
t75 = MDP(16) + MDP(18);
t74 = MDP(17) + MDP(19);
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t58 = sin(qJ(1));
t61 = cos(qJ(1));
t63 = g(1) * t61 + g(2) * t58;
t43 = -g(3) * t60 + t57 * t63;
t70 = g(3) * t57;
t68 = t58 * t60;
t56 = sin(qJ(3));
t67 = t61 * t56;
t59 = cos(qJ(3));
t66 = t61 * t59;
t54 = t59 * pkin(3) + pkin(2);
t55 = qJ(4) + pkin(6);
t64 = t54 * t60 + t55 * t57;
t49 = t58 * t59 - t60 * t67;
t47 = t56 * t68 + t66;
t53 = t56 * pkin(3) + pkin(5);
t50 = t58 * t56 + t60 * t66;
t48 = -t59 * t68 + t67;
t45 = pkin(1) + t64;
t1 = [t63 * MDP(3) + (-g(1) * (-t45 * t58 + t53 * t61) - g(2) * (t45 * t61 + t53 * t58)) * MDP(21) + t75 * (-g(1) * t48 - g(2) * t50) + t74 * (-g(1) * t47 - g(2) * t49) + (MDP(9) * t60 - t76 * t57 + MDP(2)) * (g(1) * t58 - g(2) * t61); (-g(3) * t64 - t63 * (-t57 * t54 + t55 * t60)) * MDP(21) + t76 * (t60 * t63 + t70) + (-t74 * t56 + t75 * t59 + MDP(9)) * t43; t74 * (g(1) * t50 - g(2) * t48 + t59 * t70) + (pkin(3) * MDP(21) + t75) * (-g(1) * t49 + g(2) * t47 + t56 * t70); -t43 * MDP(21);];
taug = t1;
