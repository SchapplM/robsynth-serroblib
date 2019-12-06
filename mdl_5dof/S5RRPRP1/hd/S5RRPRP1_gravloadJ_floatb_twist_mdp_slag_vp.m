% Calculate Gravitation load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:23
% EndTime: 2019-12-05 18:22:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (132->37), mult. (106->49), div. (0->0), fcn. (73->8), ass. (0->23)
t39 = qJ(1) + qJ(2);
t37 = sin(t39);
t55 = pkin(2) * t37;
t38 = cos(t39);
t54 = pkin(2) * t38;
t42 = sin(qJ(1));
t53 = t42 * pkin(1);
t44 = cos(qJ(1));
t52 = t44 * pkin(1);
t41 = sin(qJ(4));
t43 = cos(qJ(4));
t48 = g(2) * t38 + g(3) * t37;
t36 = pkin(8) + t39;
t33 = sin(t36);
t34 = cos(t36);
t49 = g(2) * t33 - g(3) * t34;
t50 = g(2) * t34 + g(3) * t33;
t51 = t49 * MDP(15) + (-g(2) * t37 + g(3) * t38) * MDP(6) + t48 * MDP(5) + (t43 * MDP(13) - t41 * MDP(14)) * t50;
t35 = t43 * pkin(4) + pkin(3);
t40 = -qJ(5) - pkin(7);
t47 = t33 * t40 - t34 * t35 - t54;
t46 = -t33 * t35 - t34 * t40 - t55;
t1 = [(g(2) * t44 + g(3) * t42) * MDP(2) + (-g(2) * t42 + g(3) * t44) * MDP(3) + (-g(2) * (-t52 - t54) - g(3) * (-t53 - t55)) * MDP(7) + (-g(2) * (t47 - t52) - g(3) * (t46 - t53)) * MDP(16) + t51; t48 * pkin(2) * MDP(7) + (-g(2) * t47 - g(3) * t46) * MDP(16) + t51; (-MDP(16) - MDP(7)) * g(1); (g(1) * t41 - t49 * t43) * MDP(14) + (MDP(16) * pkin(4) + MDP(13)) * (-g(1) * t43 - t49 * t41); -t50 * MDP(16);];
taug = t1;
