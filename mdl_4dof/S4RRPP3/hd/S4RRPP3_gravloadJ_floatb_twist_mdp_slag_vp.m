% Calculate Gravitation load on the joints for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RRPP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:36:06
% EndTime: 2021-01-15 10:36:07
% DurationCPUTime: 0.15s
% Computational Cost: add. (106->43), mult. (152->54), div. (0->0), fcn. (121->8), ass. (0->23)
t59 = MDP(11) + MDP(15);
t58 = MDP(12) - MDP(17);
t46 = qJ(2) + pkin(6);
t43 = cos(t46);
t57 = g(3) * t43;
t49 = -qJ(3) - pkin(5);
t51 = sin(qJ(1));
t56 = t51 * t49;
t53 = cos(qJ(1));
t55 = g(1) * t53 + g(2) * t51;
t38 = g(1) * t51 - g(2) * t53;
t52 = cos(qJ(2));
t50 = sin(qJ(2));
t48 = cos(pkin(6));
t47 = sin(pkin(6));
t44 = t52 * pkin(2);
t42 = sin(t46);
t41 = t44 + pkin(1);
t40 = t49 * t53;
t37 = -t47 * pkin(3) + qJ(4) * t48;
t36 = pkin(3) * t48 + qJ(4) * t47 + pkin(2);
t29 = t36 * t52 + t37 * t50 + pkin(1);
t1 = [(-g(1) * (-t51 * t41 - t40) - g(2) * (t53 * t41 - t56)) * MDP(14) + (-g(1) * (-t29 * t51 - t40) - g(2) * (t29 * t53 - t56)) * MDP(18) + (MDP(3) - MDP(13) - MDP(16)) * t55 + (-t50 * MDP(10) + t52 * MDP(9) - t58 * t42 + t59 * t43 + MDP(2)) * t38; (g(3) * t50 + t55 * t52) * MDP(10) + (-g(3) * (t43 * pkin(3) + t42 * qJ(4) + t44) - t55 * (-t36 * t50 + t37 * t52)) * MDP(18) + (MDP(14) * pkin(2) + MDP(9)) * (-g(3) * t52 + t55 * t50) + t58 * (g(3) * t42 + t55 * t43) + t59 * (t55 * t42 - t57); (-MDP(14) - MDP(18)) * t38; (t57 - t55 * (t47 * t52 + t48 * t50)) * MDP(18);];
taug = t1;
