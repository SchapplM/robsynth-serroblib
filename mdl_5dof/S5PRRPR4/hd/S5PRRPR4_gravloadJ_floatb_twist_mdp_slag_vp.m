% Calculate Gravitation load on the joints for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:52:42
% EndTime: 2021-01-15 15:52:44
% DurationCPUTime: 0.22s
% Computational Cost: add. (152->51), mult. (201->77), div. (0->0), fcn. (189->10), ass. (0->23)
t43 = sin(pkin(8));
t44 = cos(pkin(8));
t53 = g(1) * t44 + g(2) * t43;
t47 = sin(qJ(2));
t49 = cos(qJ(2));
t34 = -g(3) * t49 + t53 * t47;
t60 = g(3) * t47;
t58 = t43 * t49;
t57 = t44 * t49;
t46 = sin(qJ(3));
t56 = t46 * t49;
t48 = cos(qJ(3));
t55 = t48 * t49;
t42 = qJ(3) + pkin(9);
t41 = qJ(5) + t42;
t36 = sin(t41);
t37 = cos(t41);
t54 = (-g(1) * (-t36 * t57 + t43 * t37) - g(2) * (-t36 * t58 - t44 * t37) + t36 * t60) * MDP(21) + (-g(1) * (-t43 * t36 - t37 * t57) - g(2) * (t44 * t36 - t37 * t58) + t37 * t60) * MDP(22);
t45 = qJ(4) + pkin(6);
t40 = cos(t42);
t39 = sin(t42);
t38 = pkin(3) * t48 + pkin(2);
t1 = [(-MDP(1) - MDP(15)) * g(3); (-g(3) * (t38 * t49 + t47 * t45) + t53 * (t38 * t47 - t45 * t49)) * MDP(15) + (MDP(4) - MDP(14)) * (t53 * t49 + t60) + (t48 * MDP(10) - t46 * MDP(11) + t40 * MDP(12) - t39 * MDP(13) + t37 * MDP(21) - t36 * MDP(22) + MDP(3)) * t34; (-g(1) * (-t43 * t46 - t44 * t55) - g(2) * (-t43 * t55 + t44 * t46) + t48 * t60) * MDP(11) + (-g(1) * (-t39 * t57 + t43 * t40) - g(2) * (-t39 * t58 - t44 * t40) + t39 * t60) * MDP(12) + (-g(1) * (-t43 * t39 - t40 * t57) - g(2) * (t44 * t39 - t40 * t58) + t40 * t60) * MDP(13) + t54 + (pkin(3) * MDP(15) + MDP(10)) * (-g(1) * (t43 * t48 - t44 * t56) - g(2) * (-t43 * t56 - t44 * t48) + t46 * t60); -t34 * MDP(15); t54;];
taug = t1;
