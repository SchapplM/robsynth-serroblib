% Calculate Gravitation load on the joints for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:22
% EndTime: 2019-12-31 19:36:23
% DurationCPUTime: 0.20s
% Computational Cost: add. (125->52), mult. (180->73), div. (0->0), fcn. (153->8), ass. (0->27)
t50 = sin(qJ(1));
t53 = cos(qJ(1));
t39 = g(1) * t53 + g(2) * t50;
t46 = qJ(2) + pkin(8);
t43 = cos(t46);
t62 = g(3) * t43;
t48 = sin(qJ(5));
t61 = t50 * t48;
t51 = cos(qJ(5));
t60 = t50 * t51;
t59 = t53 * t48;
t58 = t53 * t51;
t38 = g(1) * t50 - g(2) * t53;
t42 = sin(t46);
t57 = t43 * pkin(3) + t42 * qJ(4);
t52 = cos(qJ(2));
t49 = sin(qJ(2));
t47 = -qJ(3) - pkin(6);
t44 = t52 * pkin(2);
t41 = t44 + pkin(1);
t40 = t53 * t41;
t37 = -t42 * t61 + t58;
t36 = t42 * t60 + t59;
t35 = t42 * t59 + t60;
t34 = t42 * t58 - t61;
t33 = -t39 * t42 + t62;
t1 = [(-g(1) * (-t50 * t41 - t53 * t47) - g(2) * (-t50 * t47 + t40)) * MDP(12) + (-g(2) * t40 + (g(1) * t47 - g(2) * t57) * t53 + (-g(1) * (-t41 - t57) + g(2) * t47) * t50) * MDP(16) + (-g(1) * t37 - g(2) * t35) * MDP(22) + (g(1) * t36 - g(2) * t34) * MDP(23) + (MDP(3) - MDP(11) - MDP(13)) * t39 + (-t49 * MDP(10) - t43 * MDP(14) + t42 * MDP(15) + t52 * MDP(9) + MDP(2)) * t38; (g(3) * t49 + t39 * t52) * MDP(10) + t33 * MDP(14) + (-g(3) * (t44 + t57) + t39 * (pkin(2) * t49 + pkin(3) * t42 - qJ(4) * t43)) * MDP(16) + (pkin(2) * MDP(12) + MDP(9)) * (-g(3) * t52 + t39 * t49) + (MDP(22) * t48 + MDP(23) * t51 + MDP(15)) * (-g(3) * t42 - t39 * t43); (-MDP(12) - MDP(16)) * t38; t33 * MDP(16); (-g(1) * t34 - g(2) * t36 + t51 * t62) * MDP(22) + (g(1) * t35 - g(2) * t37 - t48 * t62) * MDP(23);];
taug = t1;
