% Calculate Gravitation load on the joints for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRRPP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:13
% EndTime: 2019-12-31 17:41:13
% DurationCPUTime: 0.19s
% Computational Cost: add. (147->41), mult. (165->52), div. (0->0), fcn. (129->4), ass. (0->25)
t69 = MDP(10) + MDP(12) + MDP(16);
t68 = MDP(11) - MDP(14) - MDP(17);
t52 = pkin(7) + qJ(2);
t47 = sin(t52);
t48 = cos(t52);
t39 = g(1) * t48 + g(2) * t47;
t53 = sin(qJ(3));
t67 = t39 * t53;
t49 = t53 * qJ(4);
t54 = cos(qJ(3));
t59 = t54 * pkin(3) + t49;
t65 = pkin(3) * t53;
t64 = g(1) * t47;
t61 = t54 * pkin(4);
t60 = t48 * t54;
t58 = qJ(4) * t54;
t57 = -MDP(15) - MDP(19);
t56 = pkin(3) * t60 + t47 * pkin(6) + (pkin(2) + t49) * t48;
t38 = -g(2) * t48 + t64;
t55 = -pkin(2) - t59;
t45 = t48 * pkin(6);
t42 = t48 * t58;
t40 = t47 * t58;
t34 = -g(3) * t54 + t67;
t1 = [(-MDP(1) + t57) * g(3); (-g(1) * t45 - g(2) * t56 - t55 * t64) * MDP(15) + (-g(1) * (-t48 * qJ(5) + t45) - g(2) * (pkin(4) * t60 + t56) + (-g(1) * (t55 - t61) + g(2) * qJ(5)) * t47) * MDP(19) + (MDP(4) - MDP(13) + MDP(18)) * t39 + (-t68 * t53 + t69 * t54 + MDP(3)) * t38; (-g(1) * (-t48 * t65 + t42) - g(2) * (-t47 * t65 + t40) - g(3) * t59) * MDP(15) + (-g(1) * t42 - g(2) * t40 - g(3) * (t59 + t61) + (pkin(3) + pkin(4)) * t67) * MDP(19) + t68 * (g(3) * t53 + t39 * t54) + t69 * t34; t57 * t34; t38 * MDP(19);];
taug = t1;
