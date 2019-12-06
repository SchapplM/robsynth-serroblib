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
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:18
% EndTime: 2019-12-05 16:49:19
% DurationCPUTime: 0.19s
% Computational Cost: add. (123->45), mult. (185->71), div. (0->0), fcn. (169->8), ass. (0->23)
t44 = sin(pkin(8));
t45 = cos(pkin(8));
t53 = g(1) * t45 + g(2) * t44;
t47 = sin(qJ(2));
t49 = cos(qJ(2));
t33 = -g(3) * t49 + t53 * t47;
t60 = g(3) * t47;
t58 = t44 * t49;
t57 = t45 * t49;
t46 = sin(qJ(3));
t56 = t46 * t49;
t48 = cos(qJ(3));
t55 = t48 * t49;
t43 = qJ(3) + qJ(4);
t39 = sin(t43);
t40 = cos(t43);
t50 = -g(1) * (-t39 * t57 + t44 * t40) - g(2) * (-t39 * t58 - t45 * t40) + t39 * t60;
t54 = t50 * MDP(17) + (-g(1) * (-t44 * t39 - t40 * t57) - g(2) * (t45 * t39 - t40 * t58) + t40 * t60) * MDP(18);
t37 = t48 * pkin(3) + pkin(4) * t40;
t42 = -qJ(5) - pkin(7) - pkin(6);
t36 = -t46 * pkin(3) - pkin(4) * t39;
t35 = pkin(2) + t37;
t1 = [(-MDP(1) - MDP(20)) * g(3); (-g(3) * (t49 * t35 - t47 * t42) + t53 * (t35 * t47 + t42 * t49)) * MDP(20) + (MDP(4) - MDP(19)) * (t53 * t49 + t60) + (MDP(10) * t48 - MDP(11) * t46 + MDP(17) * t40 - MDP(18) * t39 + MDP(3)) * t33; (-g(1) * (t44 * t48 - t45 * t56) - g(2) * (-t44 * t56 - t45 * t48) + t46 * t60) * MDP(10) + (-g(1) * (-t44 * t46 - t45 * t55) - g(2) * (-t44 * t55 + t45 * t46) + t48 * t60) * MDP(11) + (-g(1) * (t36 * t57 + t44 * t37) - g(2) * (t36 * t58 - t45 * t37) - t36 * t60) * MDP(20) + t54; t50 * MDP(20) * pkin(4) + t54; -t33 * MDP(20);];
taug = t1;
