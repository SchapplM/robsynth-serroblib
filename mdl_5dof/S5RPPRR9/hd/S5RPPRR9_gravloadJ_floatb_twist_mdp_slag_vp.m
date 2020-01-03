% Calculate Gravitation load on the joints for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:48
% EndTime: 2019-12-31 18:02:49
% DurationCPUTime: 0.15s
% Computational Cost: add. (88->40), mult. (185->64), div. (0->0), fcn. (200->8), ass. (0->22)
t48 = sin(qJ(4));
t65 = g(3) * t48;
t63 = cos(qJ(1));
t62 = sin(qJ(1));
t47 = sin(qJ(5));
t50 = cos(qJ(4));
t61 = t47 * t50;
t49 = cos(qJ(5));
t60 = t49 * t50;
t59 = t63 * pkin(1) + t62 * qJ(2);
t58 = cos(pkin(8));
t57 = sin(pkin(8));
t56 = -pkin(1) * t62 + t63 * qJ(2);
t36 = -t57 * t62 - t58 * t63;
t37 = t57 * t63 - t58 * t62;
t54 = g(1) * t36 + g(2) * t37;
t53 = t36 * t47 + t37 * t60;
t52 = -t36 * t49 + t37 * t61;
t38 = g(1) * t62 - g(2) * t63;
t35 = -t36 * t60 + t37 * t47;
t34 = t36 * t61 + t37 * t49;
t1 = [(-g(1) * t56 - g(2) * t59) * MDP(6) + t54 * MDP(8) + (-g(1) * (-pkin(2) * t62 + t56) - g(2) * (pkin(2) * t63 + t59)) * MDP(9) + (-g(1) * t53 - g(2) * t35) * MDP(22) + (g(1) * t52 - g(2) * t34) * MDP(23) + (MDP(3) - MDP(5)) * (g(1) * t63 + g(2) * t62) + (MDP(2) + MDP(4)) * t38 + (-MDP(15) * t50 + MDP(16) * t48 - MDP(7)) * (g(1) * t37 - g(2) * t36); (-MDP(6) - MDP(9)) * t38; g(3) * MDP(9); (-t50 * t54 - t65) * MDP(16) + (-MDP(22) * t49 + MDP(23) * t47 - MDP(15)) * (-g(3) * t50 + t48 * t54); (-g(1) * t34 - g(2) * t52 - t47 * t65) * MDP(22) + (g(1) * t35 - g(2) * t53 - t49 * t65) * MDP(23);];
taug = t1;
