% Calculate Gravitation load on the joints for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:34
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:33:43
% EndTime: 2021-01-15 16:33:45
% DurationCPUTime: 0.21s
% Computational Cost: add. (175->47), mult. (247->71), div. (0->0), fcn. (233->8), ass. (0->25)
t71 = MDP(17) + MDP(19);
t70 = MDP(18) + MDP(20);
t52 = sin(pkin(8));
t53 = cos(pkin(8));
t60 = g(1) * t53 + g(2) * t52;
t55 = sin(qJ(2));
t57 = cos(qJ(2));
t41 = -g(3) * t57 + t60 * t55;
t67 = g(3) * t55;
t65 = t52 * t57;
t64 = t53 * t57;
t54 = sin(qJ(3));
t63 = t54 * t57;
t56 = cos(qJ(3));
t62 = t56 * t57;
t51 = qJ(3) + qJ(4);
t48 = cos(t51);
t45 = t56 * pkin(3) + pkin(4) * t48;
t47 = sin(t51);
t37 = -g(1) * (-t47 * t64 + t52 * t48) - g(2) * (-t47 * t65 - t53 * t48) + t47 * t67;
t61 = t70 * (-g(1) * (-t52 * t47 - t48 * t64) - g(2) * (t53 * t47 - t48 * t65) + t48 * t67) + t71 * t37;
t50 = -qJ(5) - pkin(7) - pkin(6);
t44 = -pkin(3) * t54 - pkin(4) * t47;
t43 = pkin(2) + t45;
t1 = [(-MDP(1) - MDP(22)) * g(3); (-g(3) * (t43 * t57 - t55 * t50) + t60 * (t43 * t55 + t50 * t57)) * MDP(22) + (MDP(4) - MDP(21)) * (t60 * t57 + t67) + (MDP(10) * t56 - MDP(11) * t54 - t70 * t47 + t71 * t48 + MDP(3)) * t41; (-g(1) * (t52 * t56 - t53 * t63) - g(2) * (-t52 * t63 - t53 * t56) + t54 * t67) * MDP(10) + (-g(1) * (-t52 * t54 - t53 * t62) - g(2) * (-t52 * t62 + t53 * t54) + t56 * t67) * MDP(11) + (-g(1) * (t44 * t64 + t52 * t45) - g(2) * (t44 * t65 - t53 * t45) - t44 * t67) * MDP(22) + t61; t37 * pkin(4) * MDP(22) + t61; -t41 * MDP(22);];
taug = t1;
