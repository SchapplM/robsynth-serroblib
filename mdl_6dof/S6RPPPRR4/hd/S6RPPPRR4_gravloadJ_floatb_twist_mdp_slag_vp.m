% Calculate Gravitation load on the joints for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPPRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:52
% EndTime: 2019-03-09 01:35:53
% DurationCPUTime: 0.21s
% Computational Cost: add. (116->50), mult. (231->72), div. (0->0), fcn. (248->8), ass. (0->25)
t56 = cos(qJ(5));
t73 = g(3) * t56;
t72 = cos(qJ(1));
t71 = sin(qJ(1));
t53 = sin(qJ(6));
t54 = sin(qJ(5));
t70 = t53 * t54;
t55 = cos(qJ(6));
t69 = t54 * t55;
t68 = t72 * pkin(1) + t71 * qJ(2);
t67 = cos(pkin(9));
t66 = sin(pkin(9));
t65 = MDP(12) + MDP(9);
t64 = t72 * pkin(2) + t68;
t63 = -t71 * pkin(1) + t72 * qJ(2);
t41 = -t71 * t66 - t72 * t67;
t42 = t72 * t66 - t71 * t67;
t62 = g(1) * t42 - g(2) * t41;
t60 = t41 * t69 + t42 * t53;
t59 = t41 * t70 - t42 * t55;
t58 = -t71 * pkin(2) + t63;
t43 = g(1) * t71 - g(2) * t72;
t39 = -t41 * t53 + t42 * t69;
t38 = -t41 * t55 - t42 * t70;
t1 = [(-g(1) * t63 - g(2) * t68) * MDP(6) + (-g(1) * t58 - g(2) * t64) * MDP(9) + (-g(1) * (t42 * pkin(3) + t41 * qJ(4) + t58) - g(2) * (-t41 * pkin(3) + qJ(4) * t42 + t64)) * MDP(12) + (-g(1) * t60 - g(2) * t39) * MDP(25) + (g(1) * t59 - g(2) * t38) * MDP(26) + (MDP(3) - MDP(5)) * (g(1) * t72 + g(2) * t71) + (MDP(2) + MDP(4)) * t43 + (-t54 * MDP(18) - MDP(19) * t56 - MDP(11) + MDP(8)) * (g(1) * t41 + g(2) * t42) + (-MDP(7) + MDP(10)) * t62; (-MDP(6) - t65) * t43; t65 * g(3); -t62 * MDP(12); (t62 * t54 - t73) * MDP(19) + (-MDP(25) * t55 + MDP(26) * t53 - MDP(18)) * (g(3) * t54 + t62 * t56); (-g(1) * t38 - g(2) * t59 - t53 * t73) * MDP(25) + (g(1) * t39 - g(2) * t60 - t55 * t73) * MDP(26);];
taug  = t1;
