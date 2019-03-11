% Calculate Gravitation load on the joints for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRRP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:24
% EndTime: 2019-03-09 02:16:25
% DurationCPUTime: 0.40s
% Computational Cost: add. (216->71), mult. (307->90), div. (0->0), fcn. (285->8), ass. (0->30)
t68 = sin(qJ(1));
t70 = cos(qJ(1));
t90 = -g(1) * t68 + g(2) * t70;
t89 = MDP(17) - MDP(26);
t88 = MDP(23) + MDP(25);
t87 = MDP(24) - MDP(27);
t63 = pkin(9) + qJ(4);
t58 = cos(t63);
t82 = g(3) * t58;
t67 = sin(qJ(5));
t81 = t68 * t67;
t69 = cos(qJ(5));
t80 = t68 * t69;
t79 = t70 * t67;
t78 = t70 * t69;
t77 = t70 * pkin(1) + t68 * qJ(2);
t76 = -MDP(10) - MDP(28);
t55 = g(1) * t70 + g(2) * t68;
t74 = pkin(5) * t69 + qJ(6) * t67 + pkin(4);
t57 = sin(t63);
t64 = sin(pkin(9));
t73 = pkin(3) * t64 + t57 * pkin(4) - t58 * pkin(8);
t50 = t57 * t81 - t78;
t52 = t57 * t79 + t80;
t44 = g(1) * t50 - g(2) * t52 + t67 * t82;
t66 = -pkin(7) - qJ(3);
t60 = t70 * qJ(2);
t53 = t57 * t78 - t81;
t51 = t57 * t80 + t79;
t1 = [(-g(1) * (-t68 * pkin(1) + t60) - g(2) * t77) * MDP(6) + (-g(1) * (t60 + (-pkin(1) - qJ(3)) * t68) - g(2) * (t70 * qJ(3) + t77)) * MDP(10) + (-g(1) * (t53 * pkin(5) + t52 * qJ(6) + t60) - g(2) * (t51 * pkin(5) + t50 * qJ(6) + t77) + (-g(1) * t73 + g(2) * t66) * t70 + (-g(1) * (-pkin(1) + t66) - g(2) * t73) * t68) * MDP(28) + t87 * (g(1) * t52 + g(2) * t50) + t88 * (-g(1) * t53 - g(2) * t51) - (MDP(2) - MDP(4) + MDP(9)) * t90 + (-t57 * MDP(16) - t64 * MDP(7) - MDP(8) * cos(pkin(9)) - t89 * t58 + MDP(3) - MDP(5)) * t55; -(-MDP(6) + t76) * t90; t76 * t55; ((pkin(8) * t90 + g(3) * t74) * t57 + (-g(3) * pkin(8) + t90 * t74) * t58) * MDP(28) + t89 * (-t57 * t90 + t82) + (t87 * t67 - t88 * t69 - MDP(16)) * (-g(3) * t57 - t90 * t58); (-g(1) * (-t50 * pkin(5) + t51 * qJ(6)) - g(2) * (t52 * pkin(5) - t53 * qJ(6)) - (-pkin(5) * t67 + qJ(6) * t69) * t82) * MDP(28) + t87 * (g(1) * t51 - g(2) * t53 + t69 * t82) + t88 * t44; -t44 * MDP(28);];
taug  = t1;
