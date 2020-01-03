% Calculate Gravitation load on the joints for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RRRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:26
% EndTime: 2019-12-31 21:49:26
% DurationCPUTime: 0.14s
% Computational Cost: add. (272->41), mult. (200->51), div. (0->0), fcn. (155->8), ass. (0->25)
t65 = sin(qJ(4));
t67 = cos(qJ(4));
t73 = t67 * pkin(4) + t65 * qJ(5);
t85 = -pkin(3) - t73;
t84 = MDP(15) + MDP(17);
t83 = MDP(16) - MDP(19);
t64 = qJ(1) + qJ(2);
t63 = qJ(3) + t64;
t59 = sin(t63);
t60 = cos(t63);
t50 = g(1) * t60 + g(2) * t59;
t82 = g(1) * t59;
t56 = t60 * pkin(8);
t61 = sin(t64);
t77 = -pkin(2) * t61 + t56;
t76 = t59 * pkin(8) - t85 * t60;
t62 = cos(t64);
t75 = pkin(2) * t62 + t76;
t71 = (-MDP(18) + MDP(9)) * t50 + (-t83 * t65 + t84 * t67 + MDP(8)) * (-g(2) * t60 + t82);
t70 = t85 * t82;
t69 = (g(1) * t61 - g(2) * t62) * MDP(5) + (g(1) * t62 + g(2) * t61) * MDP(6) + t71;
t68 = cos(qJ(1));
t66 = sin(qJ(1));
t39 = -g(3) * t67 + t50 * t65;
t1 = [(g(1) * t66 - g(2) * t68) * MDP(2) + (g(1) * t68 + g(2) * t66) * MDP(3) + (-g(1) * (-t66 * pkin(1) + t77) - g(2) * (t68 * pkin(1) + t75) - t70) * MDP(20) + t69; (-g(1) * t77 - g(2) * t75 - t70) * MDP(20) + t69; (-g(1) * t56 - g(2) * t76 - t70) * MDP(20) + t71; (-g(3) * t73 + t50 * (pkin(4) * t65 - qJ(5) * t67)) * MDP(20) + t83 * (g(3) * t65 + t50 * t67) + t84 * t39; -t39 * MDP(20);];
taug = t1;
