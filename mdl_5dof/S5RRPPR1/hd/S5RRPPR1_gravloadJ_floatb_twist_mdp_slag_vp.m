% Calculate Gravitation load on the joints for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:43
% EndTime: 2022-01-20 09:51:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (143->36), mult. (103->48), div. (0->0), fcn. (72->9), ass. (0->22)
t49 = qJ(1) + qJ(2);
t45 = sin(t49);
t60 = pkin(2) * t45;
t51 = sin(qJ(1));
t59 = t51 * pkin(1);
t44 = pkin(8) + t49;
t39 = sin(t44);
t40 = cos(t44);
t46 = cos(t49);
t41 = pkin(2) * t46;
t58 = t40 * pkin(3) + t39 * qJ(4) + t41;
t57 = g(1) * t40 + g(2) * t39;
t56 = g(1) * t39 - g(2) * t40;
t55 = g(1) * t45 - g(2) * t46;
t48 = pkin(9) + qJ(5);
t42 = sin(t48);
t43 = cos(t48);
t54 = -t57 * MDP(9) + t55 * MDP(5) + (g(1) * t46 + g(2) * t45) * MDP(6) + (t43 * MDP(16) - t42 * MDP(17) + cos(pkin(9)) * MDP(8)) * t56;
t53 = -t39 * pkin(3) + t40 * qJ(4) - t60;
t52 = cos(qJ(1));
t47 = t52 * pkin(1);
t1 = [(g(1) * t51 - g(2) * t52) * MDP(2) + (g(1) * t52 + g(2) * t51) * MDP(3) + (-g(1) * (-t59 - t60) - g(2) * (t41 + t47)) * MDP(7) + (-g(1) * (t53 - t59) - g(2) * (t47 + t58)) * MDP(10) + t54; t55 * pkin(2) * MDP(7) + (-g(1) * t53 - g(2) * t58) * MDP(10) + t54; (-MDP(10) - MDP(7)) * g(3); -t56 * MDP(10); (-g(3) * t43 + t57 * t42) * MDP(16) + (g(3) * t42 + t57 * t43) * MDP(17);];
taug = t1;
