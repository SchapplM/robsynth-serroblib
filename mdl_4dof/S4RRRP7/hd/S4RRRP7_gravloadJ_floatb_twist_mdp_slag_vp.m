% Calculate Gravitation load on the joints for
% S4RRRP7
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
%   see S4RRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:05
% EndTime: 2019-12-31 17:21:05
% DurationCPUTime: 0.23s
% Computational Cost: add. (106->49), mult. (251->74), div. (0->0), fcn. (245->6), ass. (0->22)
t46 = sin(qJ(1));
t49 = cos(qJ(1));
t55 = g(1) * t49 + g(2) * t46;
t65 = MDP(10) - MDP(19);
t64 = MDP(16) + MDP(18);
t63 = MDP(17) - MDP(20);
t45 = sin(qJ(2));
t60 = g(3) * t45;
t48 = cos(qJ(2));
t59 = t46 * t48;
t44 = sin(qJ(3));
t58 = t49 * t44;
t47 = cos(qJ(3));
t57 = t49 * t47;
t53 = t48 * pkin(2) + t45 * pkin(6) + pkin(1);
t52 = pkin(3) * t47 + qJ(4) * t44 + pkin(2);
t37 = t44 * t59 + t57;
t39 = -t46 * t47 + t48 * t58;
t32 = g(1) * t39 + g(2) * t37 + t44 * t60;
t40 = t46 * t44 + t48 * t57;
t38 = t47 * t59 - t58;
t1 = [t55 * MDP(3) + (-g(1) * (-t38 * pkin(3) - t37 * qJ(4)) - g(2) * (t40 * pkin(3) + t39 * qJ(4)) + (-g(1) * pkin(5) - g(2) * t53) * t49 + (-g(2) * pkin(5) + g(1) * t53) * t46) * MDP(21) - t63 * (g(1) * t37 - g(2) * t39) + t64 * (g(1) * t38 - g(2) * t40) + (t48 * MDP(9) - t65 * t45 + MDP(2)) * (g(1) * t46 - g(2) * t49); ((-pkin(6) * t55 - g(3) * t52) * t48 + (-g(3) * pkin(6) + t55 * t52) * t45) * MDP(21) + t65 * (t48 * t55 + t60) + (-t63 * t44 + t64 * t47 + MDP(9)) * (-g(3) * t48 + t45 * t55); (-g(1) * (-t39 * pkin(3) + t40 * qJ(4)) - g(2) * (-t37 * pkin(3) + t38 * qJ(4)) - (-pkin(3) * t44 + qJ(4) * t47) * t60) * MDP(21) + t63 * (g(1) * t40 + g(2) * t38 + t47 * t60) + t64 * t32; -t32 * MDP(21);];
taug = t1;
