% Calculate Gravitation load on the joints for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:59:36
% EndTime: 2021-01-15 22:59:37
% DurationCPUTime: 0.12s
% Computational Cost: add. (200->39), mult. (164->51), div. (0->0), fcn. (125->10), ass. (0->22)
t58 = qJ(1) + qJ(2);
t55 = sin(t58);
t56 = cos(t58);
t44 = g(2) * t55 - g(3) * t56;
t57 = qJ(3) + pkin(9);
t54 = qJ(5) + t57;
t49 = sin(t54);
t50 = cos(t54);
t68 = (-g(1) * t50 + t44 * t49) * MDP(23) + (g(1) * t49 + t44 * t50) * MDP(24);
t62 = cos(qJ(3));
t51 = t62 * pkin(3) + pkin(2);
t59 = -qJ(4) - pkin(7);
t67 = t55 * t51 + t56 * t59;
t66 = t56 * t51 - t55 * t59;
t45 = g(2) * t56 + g(3) * t55;
t52 = sin(t57);
t53 = cos(t57);
t60 = sin(qJ(3));
t65 = (-MDP(16) + MDP(6)) * t44 + (-MDP(12) * t62 + t60 * MDP(13) - t53 * MDP(14) + t52 * MDP(15) - t50 * MDP(23) + t49 * MDP(24) - MDP(5)) * t45;
t63 = cos(qJ(1));
t61 = sin(qJ(1));
t1 = [(-g(2) * t63 - g(3) * t61) * MDP(2) + (g(2) * t61 - g(3) * t63) * MDP(3) + (-g(2) * (t63 * pkin(1) + t66) - g(3) * (t61 * pkin(1) + t67)) * MDP(17) + t65; (-g(2) * t66 - g(3) * t67) * MDP(17) + t65; (g(1) * t60 + t44 * t62) * MDP(13) + (-g(1) * t53 + t44 * t52) * MDP(14) + (g(1) * t52 + t44 * t53) * MDP(15) + t68 + (MDP(17) * pkin(3) + MDP(12)) * (-g(1) * t62 + t44 * t60); t45 * MDP(17); t68;];
taug = t1;
