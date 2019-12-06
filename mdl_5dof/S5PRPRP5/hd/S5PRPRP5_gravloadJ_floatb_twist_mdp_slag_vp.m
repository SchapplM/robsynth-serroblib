% Calculate Gravitation load on the joints for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:42
% EndTime: 2019-12-05 15:38:43
% DurationCPUTime: 0.22s
% Computational Cost: add. (146->43), mult. (219->60), div. (0->0), fcn. (201->8), ass. (0->24)
t71 = MDP(14) + MDP(16);
t70 = MDP(15) - MDP(18);
t52 = sin(pkin(7));
t54 = cos(pkin(7));
t62 = g(1) * t54 + g(2) * t52;
t56 = sin(qJ(2));
t57 = cos(qJ(2));
t42 = -g(3) * t57 + t62 * t56;
t67 = g(3) * t56;
t65 = t52 * t57;
t64 = t54 * t57;
t63 = -MDP(19) - MDP(8);
t50 = pkin(8) + qJ(4);
t47 = sin(t50);
t48 = cos(t50);
t53 = cos(pkin(8));
t60 = t53 * pkin(3) + pkin(4) * t48 + qJ(5) * t47 + pkin(2);
t38 = t47 * t65 + t54 * t48;
t40 = t47 * t64 - t52 * t48;
t35 = g(1) * t40 + g(2) * t38 + t47 * t67;
t55 = -pkin(6) - qJ(3);
t41 = t52 * t47 + t48 * t64;
t39 = -t54 * t47 + t48 * t65;
t1 = [(-MDP(1) + t63) * g(3); (-g(3) * (t57 * pkin(2) + t56 * qJ(3)) + t62 * (pkin(2) * t56 - qJ(3) * t57)) * MDP(8) + (-g(3) * (-t56 * t55 + t60 * t57) + t62 * (t55 * t57 + t60 * t56)) * MDP(19) + (MDP(4) - MDP(7) - MDP(17)) * (t62 * t57 + t67) + (t53 * MDP(5) - sin(pkin(8)) * MDP(6) - t70 * t47 + t71 * t48 + MDP(3)) * t42; t63 * t42; (-g(1) * (-t40 * pkin(4) + t41 * qJ(5)) - g(2) * (-t38 * pkin(4) + t39 * qJ(5)) - (-pkin(4) * t47 + qJ(5) * t48) * t67) * MDP(19) + t70 * (g(1) * t41 + g(2) * t39 + t48 * t67) + t71 * t35; -t35 * MDP(19);];
taug = t1;
