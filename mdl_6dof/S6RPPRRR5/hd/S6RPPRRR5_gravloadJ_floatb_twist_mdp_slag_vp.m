% Calculate Gravitation load on the joints for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPPRRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:29:01
% EndTime: 2019-03-09 02:29:01
% DurationCPUTime: 0.16s
% Computational Cost: add. (116->45), mult. (175->60), div. (0->0), fcn. (150->8), ass. (0->24)
t48 = qJ(4) + qJ(5);
t43 = cos(t48);
t62 = g(3) * t43;
t49 = sin(qJ(6));
t51 = sin(qJ(1));
t61 = t51 * t49;
t52 = cos(qJ(6));
t60 = t51 * t52;
t54 = cos(qJ(1));
t59 = t54 * t49;
t58 = t54 * t52;
t57 = t54 * pkin(1) + t51 * qJ(2);
t41 = g(1) * t54 + g(2) * t51;
t42 = sin(t48);
t56 = (t41 * t42 + t62) * MDP(23) + (-t52 * MDP(29) + t49 * MDP(30) - MDP(22)) * (-g(3) * t42 + t41 * t43);
t40 = g(1) * t51 - g(2) * t54;
t53 = cos(qJ(4));
t50 = sin(qJ(4));
t45 = t54 * qJ(2);
t39 = t42 * t58 - t61;
t38 = -t42 * t59 - t60;
t37 = -t42 * t60 - t59;
t36 = t42 * t61 - t58;
t1 = [(-g(1) * (-t51 * pkin(1) + t45) - g(2) * t57) * MDP(6) + (-g(1) * (t45 + (-pkin(1) - qJ(3)) * t51) - g(2) * (t54 * qJ(3) + t57)) * MDP(9) + (-g(1) * t37 - g(2) * t39) * MDP(29) + (-g(1) * t36 - g(2) * t38) * MDP(30) + (MDP(3) - MDP(5) - MDP(7)) * t41 + (t50 * MDP(15) + t53 * MDP(16) + MDP(22) * t42 + MDP(23) * t43 + MDP(2) - MDP(4) + MDP(8)) * t40; (-MDP(6) - MDP(9)) * t40; -t41 * MDP(9); (g(3) * t50 - t41 * t53) * MDP(15) + (g(3) * t53 + t41 * t50) * MDP(16) + t56; t56; (-g(1) * t38 + g(2) * t36 + t49 * t62) * MDP(29) + (g(1) * t39 - g(2) * t37 + t52 * t62) * MDP(30);];
taug  = t1;
