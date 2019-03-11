% Calculate Gravitation load on the joints for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6RPPRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:58
% EndTime: 2019-03-09 01:44:59
% DurationCPUTime: 0.17s
% Computational Cost: add. (140->49), mult. (145->70), div. (0->0), fcn. (119->10), ass. (0->29)
t56 = sin(qJ(4));
t72 = pkin(4) * t56;
t52 = qJ(4) + pkin(10);
t49 = cos(t52);
t70 = g(3) * t49;
t53 = qJ(1) + pkin(9);
t48 = sin(t53);
t55 = sin(qJ(6));
t69 = t48 * t55;
t58 = cos(qJ(6));
t68 = t48 * t58;
t50 = cos(t53);
t67 = t50 * t55;
t66 = t50 * t58;
t65 = -MDP(16) - MDP(7);
t60 = cos(qJ(1));
t64 = t60 * pkin(1) + t50 * pkin(2) + t48 * qJ(3);
t57 = sin(qJ(1));
t63 = -t57 * pkin(1) + t50 * qJ(3);
t42 = -g(1) * t50 - g(2) * t48;
t41 = g(1) * t48 - g(2) * t50;
t59 = cos(qJ(4));
t54 = -qJ(5) - pkin(7);
t47 = sin(t52);
t40 = t47 * t66 - t69;
t39 = t47 * t67 + t68;
t38 = t47 * t68 + t67;
t37 = -t47 * t69 + t66;
t1 = [(g(1) * t60 + g(2) * t57) * MDP(3) + (-g(1) * (-t48 * pkin(2) + t63) - g(2) * t64) * MDP(7) + (-g(1) * (t50 * t72 + (-pkin(2) + t54) * t48 + t63) - g(2) * (t48 * t72 - t50 * t54 + t64)) * MDP(16) + (-g(1) * t40 - g(2) * t38) * MDP(22) + (g(1) * t39 - g(2) * t37) * MDP(23) + (-MDP(5) + MDP(15)) * t41 + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t57 - g(2) * t60) + (MDP(13) * t56 + MDP(14) * t59 + MDP(6)) * t42; (-MDP(4) + t65) * g(3); t65 * t41; (g(3) * t59 + t41 * t56) * MDP(14) + (MDP(16) * pkin(4) + MDP(13)) * (g(3) * t56 - t41 * t59) + (-MDP(22) * t58 + MDP(23) * t55) * (-g(3) * t47 + t41 * t49); t42 * MDP(16); (-g(1) * t37 - g(2) * t39 + t55 * t70) * MDP(22) + (g(1) * t38 - g(2) * t40 + t58 * t70) * MDP(23);];
taug  = t1;
