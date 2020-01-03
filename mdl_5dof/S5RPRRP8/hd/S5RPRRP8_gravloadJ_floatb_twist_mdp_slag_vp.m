% Calculate Gravitation load on the joints for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRRP8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:30
% EndTime: 2019-12-31 18:47:31
% DurationCPUTime: 0.20s
% Computational Cost: add. (140->39), mult. (279->52), div. (0->0), fcn. (299->6), ass. (0->20)
t68 = sin(qJ(3));
t69 = sin(qJ(1));
t70 = cos(qJ(3));
t71 = cos(qJ(1));
t49 = -t69 * t68 - t71 * t70;
t50 = t71 * t68 - t69 * t70;
t60 = sin(qJ(4));
t61 = cos(qJ(4));
t64 = -t61 * pkin(4) - t60 * qJ(5);
t62 = pkin(3) - t64;
t82 = -(-g(1) * pkin(7) + g(2) * t62) * t49 + (g(2) * pkin(7) + g(1) * t62) * t50;
t48 = g(1) * t49 + g(2) * t50;
t79 = MDP(16) - MDP(19);
t80 = MDP(15) + MDP(17);
t81 = (-t79 * t60 + t61 * t80 + MDP(8)) * (g(1) * t50 - g(2) * t49) + (MDP(18) - MDP(9)) * t48;
t67 = t71 * pkin(1) + t69 * qJ(2);
t66 = -t69 * pkin(1) + t71 * qJ(2);
t51 = g(1) * t69 - g(2) * t71;
t41 = g(3) * t61 - t48 * t60;
t1 = [(-g(1) * t66 - g(2) * t67) * MDP(6) + (-g(1) * (-t69 * pkin(2) + t66) - g(2) * (t71 * pkin(2) + t67) - t82) * MDP(20) + (MDP(3) - MDP(5)) * (g(1) * t71 + g(2) * t69) + (MDP(2) + MDP(4)) * t51 - t81; (-MDP(20) - MDP(6)) * t51; t82 * MDP(20) + t81; (-g(3) * t64 - t48 * (pkin(4) * t60 - qJ(5) * t61)) * MDP(20) + t79 * (-g(3) * t60 - t48 * t61) + t80 * t41; -t41 * MDP(20);];
taug = t1;
