% Calculate Gravitation load on the joints for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPPRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPPRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:42
% EndTime: 2019-03-09 01:37:43
% DurationCPUTime: 0.17s
% Computational Cost: add. (106->49), mult. (211->70), div. (0->0), fcn. (218->8), ass. (0->26)
t47 = sin(qJ(5));
t65 = g(3) * t47;
t46 = sin(qJ(6));
t50 = cos(qJ(5));
t63 = t46 * t50;
t49 = cos(qJ(6));
t62 = t49 * t50;
t48 = sin(qJ(1));
t51 = cos(qJ(1));
t61 = t51 * pkin(1) + t48 * qJ(2);
t60 = sin(pkin(9));
t59 = -MDP(12) - MDP(9);
t58 = t51 * qJ(3) + t61;
t45 = cos(pkin(9));
t35 = t51 * t45 - t48 * t60;
t36 = t48 * t45 + t51 * t60;
t57 = g(1) * t36 - g(2) * t35;
t42 = t51 * qJ(2);
t55 = t42 + (-pkin(1) - qJ(3)) * t48;
t54 = t35 * t62 + t36 * t46;
t53 = t35 * t63 - t36 * t49;
t38 = g(1) * t51 + g(2) * t48;
t37 = g(1) * t48 - g(2) * t51;
t34 = -t35 * t46 + t36 * t62;
t33 = -t35 * t49 - t36 * t63;
t1 = [(-g(1) * (-t48 * pkin(1) + t42) - g(2) * t61) * MDP(6) + (-g(1) * t55 - g(2) * t58) * MDP(9) + t57 * MDP(11) + (-g(1) * (t51 * pkin(3) + t55) - g(2) * (t48 * pkin(3) + t58)) * MDP(12) + (-g(1) * t54 - g(2) * t34) * MDP(25) + (g(1) * t53 - g(2) * t33) * MDP(26) + (MDP(3) - MDP(5) - MDP(7)) * t38 + (MDP(2) - MDP(4) + MDP(8)) * t37 + (-t50 * MDP(18) + MDP(19) * t47 - MDP(10)) * (g(1) * t35 + g(2) * t36); (-MDP(6) + t59) * t37; t59 * t38; -g(3) * MDP(12); (t57 * t50 + t65) * MDP(19) + (MDP(25) * t49 - MDP(26) * t46 + MDP(18)) * (-g(3) * t50 + t57 * t47); (-g(1) * t33 - g(2) * t53 + t46 * t65) * MDP(25) + (g(1) * t34 - g(2) * t54 + t49 * t65) * MDP(26);];
taug  = t1;
