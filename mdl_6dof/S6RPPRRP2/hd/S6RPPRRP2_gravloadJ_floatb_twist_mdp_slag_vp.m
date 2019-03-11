% Calculate Gravitation load on the joints for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRRP2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:17
% EndTime: 2019-03-09 02:01:18
% DurationCPUTime: 0.33s
% Computational Cost: add. (309->70), mult. (291->94), div. (0->0), fcn. (271->10), ass. (0->33)
t65 = qJ(1) + pkin(9);
t60 = sin(t65);
t62 = cos(t65);
t79 = g(1) * t62 + g(2) * t60;
t92 = MDP(15) - MDP(24);
t91 = MDP(21) + MDP(23);
t90 = MDP(22) - MDP(25);
t64 = pkin(10) + qJ(4);
t59 = sin(t64);
t87 = g(3) * t59;
t70 = sin(qJ(1));
t86 = t70 * pkin(1);
t69 = sin(qJ(5));
t85 = t60 * t69;
t71 = cos(qJ(5));
t84 = t60 * t71;
t83 = t62 * t69;
t82 = t62 * t71;
t81 = MDP(26) + MDP(8);
t78 = g(1) * t60 - g(2) * t62;
t61 = cos(t64);
t67 = cos(pkin(10));
t76 = t67 * pkin(3) + t61 * pkin(4) + t59 * pkin(8) + pkin(2);
t75 = pkin(5) * t71 + qJ(6) * t69 + pkin(4);
t51 = t61 * t85 + t82;
t53 = t61 * t83 - t84;
t45 = g(1) * t53 + g(2) * t51 + t69 * t87;
t72 = cos(qJ(1));
t68 = -pkin(7) - qJ(3);
t63 = t72 * pkin(1);
t54 = t61 * t82 + t85;
t52 = t61 * t84 - t83;
t1 = [(g(1) * t72 + g(2) * t70) * MDP(3) - t79 * MDP(7) + (-g(1) * (-t60 * pkin(2) + t62 * qJ(3) - t86) - g(2) * (t62 * pkin(2) + t60 * qJ(3) + t63)) * MDP(8) + (-g(1) * (-t52 * pkin(5) - t51 * qJ(6) - t86) - g(2) * (t54 * pkin(5) + t53 * qJ(6) + t63) + (g(1) * t68 - g(2) * t76) * t62 + (g(1) * t76 + g(2) * t68) * t60) * MDP(26) - t90 * (g(1) * t51 - g(2) * t53) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t70 - g(2) * t72) + t91 * (g(1) * t52 - g(2) * t54) + (t61 * MDP(14) + MDP(5) * t67 - MDP(6) * sin(pkin(10)) - t92 * t59) * t78; (-MDP(4) - t81) * g(3); -t81 * t78; ((-pkin(8) * t79 - g(3) * t75) * t61 + (-g(3) * pkin(8) + t79 * t75) * t59) * MDP(26) + t92 * (t61 * t79 + t87) + (-t90 * t69 + t91 * t71 + MDP(14)) * (-g(3) * t61 + t59 * t79); (-g(1) * (-t53 * pkin(5) + t54 * qJ(6)) - g(2) * (-t51 * pkin(5) + t52 * qJ(6)) - (-pkin(5) * t69 + qJ(6) * t71) * t87) * MDP(26) + t90 * (g(1) * t54 + g(2) * t52 + t71 * t87) + t91 * t45; -t45 * MDP(26);];
taug  = t1;
