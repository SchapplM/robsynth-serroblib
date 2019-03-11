% Calculate Gravitation load on the joints for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRPR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:11
% EndTime: 2019-03-09 01:56:12
% DurationCPUTime: 0.30s
% Computational Cost: add. (146->57), mult. (202->71), div. (0->0), fcn. (169->8), ass. (0->30)
t61 = sin(qJ(1));
t63 = cos(qJ(1));
t82 = -g(1) * t61 + g(2) * t63;
t81 = -MDP(16) + MDP(19);
t80 = MDP(17) - MDP(20);
t56 = pkin(9) + qJ(4);
t50 = sin(t56);
t76 = g(3) * t50;
t60 = sin(qJ(6));
t74 = t61 * t60;
t62 = cos(qJ(6));
t73 = t61 * t62;
t72 = t63 * t60;
t71 = t63 * t62;
t70 = t63 * pkin(1) + t61 * qJ(2);
t69 = -MDP(10) - MDP(21);
t68 = g(2) * t70;
t49 = g(1) * t63 + g(2) * t61;
t51 = cos(t56);
t66 = -t50 * pkin(4) + t51 * qJ(5);
t57 = sin(pkin(9));
t65 = pkin(3) * t57 - t66;
t59 = -pkin(7) - qJ(3);
t53 = t63 * qJ(2);
t47 = -t51 * t74 + t71;
t46 = -t51 * t73 - t72;
t45 = -t51 * t72 - t73;
t44 = -t51 * t71 + t74;
t41 = -t51 * t82 - t76;
t1 = [(-g(1) * (-t61 * pkin(1) + t53) - t68) * MDP(6) + (-g(1) * (t53 + (-pkin(1) - qJ(3)) * t61) - g(2) * (t63 * qJ(3) + t70)) * MDP(10) + (-g(1) * t53 - t68 + (-g(1) * t65 + g(2) * t59) * t63 + (-g(1) * (-pkin(1) + t59) - g(2) * t65) * t61) * MDP(21) + (-g(1) * t45 - g(2) * t47) * MDP(27) + (-g(1) * t44 - g(2) * t46) * MDP(28) - (MDP(2) - MDP(4) + MDP(9) + MDP(18)) * t82 + (-t57 * MDP(7) - MDP(8) * cos(pkin(9)) + t81 * t50 - t80 * t51 + MDP(3) - MDP(5)) * t49; -(-MDP(6) + t69) * t82; t69 * t49; (-g(3) * t66 + t82 * (pkin(4) * t51 + qJ(5) * t50)) * MDP(21) + t81 * t41 + (-MDP(27) * t60 - MDP(28) * t62 + t80) * (g(3) * t51 - t50 * t82); t41 * MDP(21); (-g(1) * t46 + g(2) * t44 - t62 * t76) * MDP(27) + (g(1) * t47 - g(2) * t45 + t60 * t76) * MDP(28);];
taug  = t1;
