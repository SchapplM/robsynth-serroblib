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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:43:16
% EndTime: 2022-01-20 11:43:16
% DurationCPUTime: 0.11s
% Computational Cost: add. (200->40), mult. (164->51), div. (0->0), fcn. (125->10), ass. (0->22)
t53 = qJ(1) + qJ(2);
t50 = sin(t53);
t51 = cos(t53);
t42 = g(1) * t51 + g(2) * t50;
t52 = qJ(3) + pkin(9);
t49 = qJ(5) + t52;
t44 = sin(t49);
t45 = cos(t49);
t63 = (-g(3) * t45 + t42 * t44) * MDP(23) + (g(3) * t44 + t42 * t45) * MDP(24);
t57 = cos(qJ(3));
t46 = pkin(3) * t57 + pkin(2);
t54 = -qJ(4) - pkin(7);
t62 = t51 * t46 - t50 * t54;
t41 = g(1) * t50 - g(2) * t51;
t61 = -t46 * t50 - t51 * t54;
t47 = sin(t52);
t48 = cos(t52);
t55 = sin(qJ(3));
t60 = (-MDP(16) + MDP(6)) * t42 + (MDP(12) * t57 - t55 * MDP(13) + t48 * MDP(14) - t47 * MDP(15) + MDP(23) * t45 - MDP(24) * t44 + MDP(5)) * t41;
t58 = cos(qJ(1));
t56 = sin(qJ(1));
t1 = [(g(1) * t56 - g(2) * t58) * MDP(2) + (g(1) * t58 + g(2) * t56) * MDP(3) + (-g(1) * (-pkin(1) * t56 + t61) - g(2) * (pkin(1) * t58 + t62)) * MDP(17) + t60; (-g(1) * t61 - g(2) * t62) * MDP(17) + t60; (g(3) * t55 + t42 * t57) * MDP(13) + (-g(3) * t48 + t42 * t47) * MDP(14) + (g(3) * t47 + t42 * t48) * MDP(15) + t63 + (MDP(17) * pkin(3) + MDP(12)) * (-g(3) * t57 + t42 * t55); -t41 * MDP(17); t63;];
taug = t1;
