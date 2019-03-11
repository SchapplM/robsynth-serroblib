% Calculate Gravitation load on the joints for
% S6RPPRRP1
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
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPPRRP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:44
% EndTime: 2019-03-09 01:58:45
% DurationCPUTime: 0.28s
% Computational Cost: add. (211->59), mult. (197->79), div. (0->0), fcn. (165->10), ass. (0->33)
t62 = qJ(1) + pkin(9);
t57 = sin(t62);
t59 = cos(t62);
t76 = g(1) * t59 + g(2) * t57;
t92 = MDP(15) - MDP(23);
t61 = pkin(10) + qJ(4);
t56 = sin(t61);
t58 = cos(t61);
t46 = -g(3) * t58 + t56 * t76;
t86 = g(3) * t56;
t68 = sin(qJ(1));
t84 = t68 * pkin(1);
t67 = sin(qJ(5));
t83 = t57 * t67;
t69 = cos(qJ(5));
t82 = t57 * t69;
t81 = t59 * t67;
t80 = t59 * t69;
t79 = MDP(24) + MDP(8);
t77 = pkin(5) * t67 + pkin(7) + qJ(3);
t75 = g(1) * t57 - g(2) * t59;
t55 = t69 * pkin(5) + pkin(4);
t65 = -qJ(6) - pkin(8);
t73 = t58 * t55 - t56 * t65;
t51 = -t58 * t81 + t82;
t49 = t58 * t83 + t80;
t64 = cos(pkin(10));
t71 = t64 * pkin(3) + pkin(2) + t73;
t70 = cos(qJ(1));
t60 = t70 * pkin(1);
t52 = t58 * t80 + t83;
t50 = -t58 * t82 + t81;
t1 = [(g(1) * t70 + g(2) * t68) * MDP(3) - t76 * MDP(7) + (-g(1) * (-t57 * pkin(2) + t59 * qJ(3) - t84) - g(2) * (t59 * pkin(2) + t57 * qJ(3) + t60)) * MDP(8) + (-g(1) * t50 - g(2) * t52) * MDP(21) + (-g(1) * t49 - g(2) * t51) * MDP(22) + (g(1) * t84 - g(2) * t60 + (-g(1) * t77 - g(2) * t71) * t59 + (g(1) * t71 - g(2) * t77) * t57) * MDP(24) + (MDP(4) * pkin(1) + MDP(2)) * (g(1) * t68 - g(2) * t70) + (t58 * MDP(14) + MDP(5) * t64 - MDP(6) * sin(pkin(10)) - t92 * t56) * t75; (-MDP(4) - t79) * g(3); -t79 * t75; (-g(3) * t73 + t76 * (t55 * t56 + t58 * t65)) * MDP(24) + t92 * (t58 * t76 + t86) + (MDP(21) * t69 - MDP(22) * t67 + MDP(14)) * t46; (g(1) * t52 - g(2) * t50 + t69 * t86) * MDP(22) + (pkin(5) * MDP(24) + MDP(21)) * (-g(1) * t51 + g(2) * t49 + t67 * t86); -t46 * MDP(24);];
taug  = t1;
