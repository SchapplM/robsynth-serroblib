% Calculate Gravitation load on the joints for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:58
% EndTime: 2020-01-03 12:03:58
% DurationCPUTime: 0.06s
% Computational Cost: add. (177->32), mult. (139->43), div. (0->0), fcn. (106->10), ass. (0->18)
t54 = qJ(1) + qJ(2);
t51 = sin(t54);
t52 = cos(t54);
t41 = g(2) * t51 - g(3) * t52;
t53 = pkin(9) + qJ(4);
t50 = qJ(5) + t53;
t44 = sin(t50);
t45 = cos(t50);
t62 = (-g(1) * t45 + t41 * t44) * MDP(23) + (g(1) * t44 + t41 * t45) * MDP(24);
t61 = t52 * pkin(2) + t51 * qJ(3);
t60 = t51 * pkin(2) - t52 * qJ(3);
t42 = g(2) * t52 + g(3) * t51;
t48 = sin(t53);
t49 = cos(t53);
t59 = (MDP(6) - MDP(9)) * t41 + (-t49 * MDP(16) + t48 * MDP(17) - t45 * MDP(23) + t44 * MDP(24) - MDP(7) * cos(pkin(9)) + MDP(8) * sin(pkin(9)) - MDP(5)) * t42;
t58 = cos(qJ(1));
t57 = sin(qJ(1));
t1 = [(-g(2) * t58 - g(3) * t57) * MDP(2) + (g(2) * t57 - g(3) * t58) * MDP(3) + (-g(2) * (t58 * pkin(1) + t61) - g(3) * (t57 * pkin(1) + t60)) * MDP(10) + t59; (-g(2) * t61 - g(3) * t60) * MDP(10) + t59; t42 * MDP(10); (-g(1) * t49 + t41 * t48) * MDP(16) + (g(1) * t48 + t41 * t49) * MDP(17) + t62; t62;];
taug = t1;
