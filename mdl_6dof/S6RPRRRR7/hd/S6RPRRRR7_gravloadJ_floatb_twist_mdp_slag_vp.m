% Calculate Gravitation load on the joints for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% MDP [34x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(34,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [34 1]), ...
  'S6RPRRRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [34x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:17:55
% EndTime: 2019-03-09 07:17:56
% DurationCPUTime: 0.16s
% Computational Cost: add. (208->44), mult. (218->62), div. (0->0), fcn. (190->10), ass. (0->25)
t51 = qJ(3) + qJ(4);
t49 = qJ(5) + t51;
t46 = cos(t49);
t66 = g(3) * t46;
t52 = sin(qJ(6));
t57 = cos(qJ(1));
t65 = t52 * t57;
t54 = sin(qJ(1));
t64 = t54 * t52;
t55 = cos(qJ(6));
t63 = t54 * t55;
t62 = t55 * t57;
t43 = g(1) * t54 - g(2) * t57;
t45 = sin(t49);
t61 = (t43 * t45 + t66) * MDP(27) + (-t55 * MDP(33) + t52 * MDP(34) - MDP(26)) * (-g(3) * t45 + t43 * t46);
t47 = sin(t51);
t48 = cos(t51);
t59 = (g(3) * t47 - t43 * t48) * MDP(19) + (g(3) * t48 + t43 * t47) * MDP(20) + t61;
t56 = cos(qJ(3));
t53 = sin(qJ(3));
t42 = t45 * t62 - t64;
t41 = t45 * t65 + t63;
t40 = t45 * t63 + t65;
t39 = -t45 * t64 + t62;
t1 = [(-g(1) * (-t54 * pkin(1) + qJ(2) * t57) - g(2) * (pkin(1) * t57 + t54 * qJ(2))) * MDP(6) + (-g(1) * t42 - g(2) * t40) * MDP(33) + (g(1) * t41 - g(2) * t39) * MDP(34) + (MDP(2) - MDP(4)) * t43 + (-t53 * MDP(12) - t56 * MDP(13) - MDP(19) * t47 - MDP(20) * t48 - MDP(26) * t45 - MDP(27) * t46 + MDP(3) - MDP(5)) * (g(1) * t57 + g(2) * t54); -t43 * MDP(6); (g(3) * t53 - t43 * t56) * MDP(12) + (g(3) * t56 + t43 * t53) * MDP(13) + t59; t59; t61; (-g(1) * t39 - g(2) * t41 + t52 * t66) * MDP(33) + (g(1) * t40 - g(2) * t42 + t55 * t66) * MDP(34);];
taug  = t1;
