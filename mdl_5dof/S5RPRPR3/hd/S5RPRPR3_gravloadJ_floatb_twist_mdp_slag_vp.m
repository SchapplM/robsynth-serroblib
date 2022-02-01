% Calculate Gravitation load on the joints for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:55
% EndTime: 2022-01-23 09:20:55
% DurationCPUTime: 0.11s
% Computational Cost: add. (161->37), mult. (120->55), div. (0->0), fcn. (104->10), ass. (0->21)
t60 = g(3) * sin(pkin(9));
t49 = cos(pkin(9));
t50 = sin(qJ(5));
t59 = t49 * t50;
t52 = cos(qJ(5));
t58 = t49 * t52;
t47 = qJ(1) + pkin(8);
t46 = qJ(3) + t47;
t44 = sin(t46);
t45 = cos(t46);
t57 = t45 * pkin(3) + t44 * qJ(4);
t56 = -t44 * pkin(3) + t45 * qJ(4);
t39 = g(1) * t44 - g(2) * t45;
t32 = t44 * t59 + t45 * t52;
t33 = -t44 * t58 + t45 * t50;
t34 = t44 * t52 - t45 * t59;
t35 = t44 * t50 + t45 * t58;
t54 = (-g(1) * t32 - g(2) * t34) * MDP(17) + (-g(1) * t33 - g(2) * t35) * MDP(16) + (MDP(7) - MDP(9)) * (g(1) * t45 + g(2) * t44) + (t49 * MDP(8) + MDP(6)) * t39;
t53 = cos(qJ(1));
t51 = sin(qJ(1));
t1 = [(g(1) * t53 + g(2) * t51) * MDP(3) + (-g(1) * (-pkin(2) * sin(t47) - t51 * pkin(1) + t56) - g(2) * (pkin(2) * cos(t47) + t53 * pkin(1) + t57)) * MDP(10) + t54 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t51 - g(2) * t53); (-MDP(10) - MDP(4)) * g(3); (-g(1) * t56 - g(2) * t57) * MDP(10) + t54; -t39 * MDP(10); (-g(1) * t34 + g(2) * t32 + t50 * t60) * MDP(16) + (g(1) * t35 - g(2) * t33 + t52 * t60) * MDP(17);];
taug = t1;
