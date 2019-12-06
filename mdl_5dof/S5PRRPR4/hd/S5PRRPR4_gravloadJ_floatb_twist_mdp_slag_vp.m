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
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:23
% EndTime: 2019-12-05 16:23:24
% DurationCPUTime: 0.18s
% Computational Cost: add. (120->38), mult. (161->59), div. (0->0), fcn. (149->8), ass. (0->20)
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t46 = g(1) * t37 + g(2) * t36;
t40 = sin(qJ(2));
t42 = cos(qJ(2));
t30 = -g(3) * t42 + t46 * t40;
t53 = g(3) * t40;
t51 = t36 * t42;
t50 = t37 * t42;
t39 = sin(qJ(3));
t49 = t39 * t42;
t41 = cos(qJ(3));
t48 = t41 * t42;
t35 = qJ(3) + pkin(9) + qJ(5);
t32 = sin(t35);
t33 = cos(t35);
t47 = (-g(1) * (-t32 * t50 + t36 * t33) - g(2) * (-t32 * t51 - t37 * t33) + t32 * t53) * MDP(19) + (-g(1) * (-t36 * t32 - t33 * t50) - g(2) * (t37 * t32 - t33 * t51) + t33 * t53) * MDP(20);
t38 = -qJ(4) - pkin(6);
t34 = t41 * pkin(3) + pkin(2);
t1 = [(-MDP(1) - MDP(13)) * g(3); (-g(3) * (t42 * t34 - t40 * t38) + t46 * (t34 * t40 + t38 * t42)) * MDP(13) + (MDP(4) - MDP(12)) * (t46 * t42 + t53) + (MDP(10) * t41 - MDP(11) * t39 + MDP(19) * t33 - MDP(20) * t32 + MDP(3)) * t30; (-g(1) * (-t36 * t39 - t37 * t48) - g(2) * (-t36 * t48 + t37 * t39) + t41 * t53) * MDP(11) + t47 + (pkin(3) * MDP(13) + MDP(10)) * (-g(1) * (t36 * t41 - t37 * t49) - g(2) * (-t36 * t49 - t37 * t41) + t39 * t53); -t30 * MDP(13); t47;];
taug = t1;
