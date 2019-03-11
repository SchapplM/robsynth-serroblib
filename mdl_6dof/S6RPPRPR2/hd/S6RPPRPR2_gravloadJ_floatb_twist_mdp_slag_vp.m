% Calculate Gravitation load on the joints for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRPR2_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRPR2_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:37
% EndTime: 2019-03-09 01:42:38
% DurationCPUTime: 0.25s
% Computational Cost: add. (203->57), mult. (186->76), div. (0->0), fcn. (155->10), ass. (0->32)
t61 = qJ(1) + pkin(9);
t56 = sin(t61);
t58 = cos(t61);
t75 = g(1) * t58 + g(2) * t56;
t88 = MDP(14) - MDP(17);
t87 = MDP(15) - MDP(18);
t60 = pkin(10) + qJ(4);
t57 = cos(t60);
t82 = g(3) * t57;
t66 = sin(qJ(1));
t81 = t66 * pkin(1);
t65 = sin(qJ(6));
t80 = t56 * t65;
t67 = cos(qJ(6));
t79 = t56 * t67;
t78 = t58 * t65;
t77 = t58 * t67;
t76 = MDP(19) + MDP(8);
t74 = g(1) * t56 - g(2) * t58;
t55 = sin(t60);
t72 = t57 * pkin(4) + t55 * qJ(5);
t63 = cos(pkin(10));
t70 = t63 * pkin(3) + pkin(2) + t72;
t68 = cos(qJ(1));
t64 = -pkin(7) - qJ(3);
t59 = t68 * pkin(1);
t50 = -t55 * t80 + t77;
t49 = t55 * t79 + t78;
t48 = t55 * t78 + t79;
t47 = t55 * t77 - t80;
t43 = t75 * t55 - t82;
t1 = [(g(1) * t68 + g(2) * t66) * MDP(3) + (-g(1) * (-t56 * pkin(2) + t58 * qJ(3) - t81) - g(2) * (t58 * pkin(2) + t56 * qJ(3) + t59)) * MDP(8) + (g(1) * t81 - g(2) * t59 + (g(1) * t64 - g(2) * t70) * t58 + (g(1) * t70 + g(2) * t64) * t56) * MDP(19) + (-g(1) * t50 - g(2) * t48) * MDP(25) + (g(1) * t49 - g(2) * t47) * MDP(26) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t66 - g(2) * t68) - (MDP(7) + MDP(16)) * t75 + (MDP(5) * t63 - MDP(6) * sin(pkin(10)) - t87 * t55 + t88 * t57) * t74; (-MDP(4) - t76) * g(3); -t76 * t74; (-g(3) * t72 + t75 * (pkin(4) * t55 - qJ(5) * t57)) * MDP(19) + t88 * t43 + (-MDP(25) * t65 - MDP(26) * t67 + t87) * (g(3) * t55 + t75 * t57); -t43 * MDP(19); (-g(1) * t47 - g(2) * t49 + t67 * t82) * MDP(25) + (g(1) * t48 - g(2) * t50 - t65 * t82) * MDP(26);];
taug  = t1;
