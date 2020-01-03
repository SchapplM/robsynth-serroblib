% Calculate Gravitation load on the joints for
% S5RRPPR3
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
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:36
% EndTime: 2019-12-31 19:26:36
% DurationCPUTime: 0.08s
% Computational Cost: add. (129->35), mult. (99->47), div. (0->0), fcn. (68->8), ass. (0->20)
t44 = qJ(1) + qJ(2);
t41 = sin(t44);
t56 = pkin(2) * t41;
t46 = sin(qJ(1));
t55 = t46 * pkin(1);
t40 = pkin(8) + t44;
t37 = sin(t40);
t38 = cos(t40);
t42 = cos(t44);
t39 = pkin(2) * t42;
t54 = t38 * pkin(3) + t37 * qJ(4) + t39;
t52 = g(1) * t37 - g(2) * t38;
t51 = g(1) * t41 - g(2) * t42;
t45 = sin(qJ(5));
t47 = cos(qJ(5));
t50 = t51 * MDP(5) - t52 * MDP(8) + (g(1) * t42 + g(2) * t41) * MDP(6) + (t45 * MDP(16) + t47 * MDP(17) + MDP(9)) * (-g(1) * t38 - g(2) * t37);
t49 = -t37 * pkin(3) + t38 * qJ(4) - t56;
t48 = cos(qJ(1));
t43 = t48 * pkin(1);
t1 = [(g(1) * t46 - g(2) * t48) * MDP(2) + (g(1) * t48 + g(2) * t46) * MDP(3) + (-g(1) * (-t55 - t56) - g(2) * (t39 + t43)) * MDP(7) + (-g(1) * (t49 - t55) - g(2) * (t43 + t54)) * MDP(10) + t50; t51 * pkin(2) * MDP(7) + (-g(1) * t49 - g(2) * t54) * MDP(10) + t50; (-MDP(10) - MDP(7)) * g(3); -t52 * MDP(10); (g(3) * t45 - t52 * t47) * MDP(16) + (g(3) * t47 + t52 * t45) * MDP(17);];
taug = t1;
