% Calculate Gravitation load on the joints for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:04
% EndTime: 2019-03-08 19:45:05
% DurationCPUTime: 0.55s
% Computational Cost: add. (301->76), mult. (527->122), div. (0->0), fcn. (604->12), ass. (0->39)
t108 = MDP(14) - MDP(17);
t107 = MDP(15) - MDP(18);
t79 = sin(pkin(6));
t104 = g(3) * t79;
t78 = sin(pkin(10));
t83 = sin(qJ(2));
t85 = cos(qJ(2));
t96 = cos(pkin(10));
t97 = cos(pkin(6));
t91 = t97 * t96;
t64 = t78 * t83 - t85 * t91;
t94 = t78 * t97;
t66 = t96 * t83 + t85 * t94;
t87 = -g(1) * t66 - g(2) * t64 + t85 * t104;
t76 = pkin(11) + qJ(4);
t74 = sin(t76);
t82 = sin(qJ(6));
t103 = t74 * t82;
t84 = cos(qJ(6));
t102 = t74 * t84;
t101 = t78 * t79;
t100 = t79 * t83;
t99 = t82 * t85;
t98 = t84 * t85;
t95 = MDP(19) + MDP(8);
t93 = t79 * t96;
t65 = t78 * t85 + t83 * t91;
t67 = -t83 * t94 + t96 * t85;
t92 = g(1) * t67 + g(2) * t65;
t75 = cos(t76);
t58 = t65 * t74 + t75 * t93;
t60 = -t75 * t101 + t67 * t74;
t62 = t74 * t100 - t97 * t75;
t89 = g(1) * t60 + g(2) * t58 + g(3) * t62;
t80 = cos(pkin(11));
t63 = t75 * t100 + t97 * t74;
t61 = t74 * t101 + t67 * t75;
t59 = t65 * t75 - t74 * t93;
t1 = [(-MDP(1) - t95) * g(3); (-g(1) * (-t66 * pkin(2) + t67 * qJ(3)) - g(2) * (-t64 * pkin(2) + t65 * qJ(3)) - (pkin(2) * t85 + qJ(3) * t83) * t104) * MDP(8) + (t83 * t104 + t92) * (-pkin(8) - qJ(3)) * MDP(19) + (-g(1) * (-t66 * t103 + t67 * t84) - g(2) * (-t64 * t103 + t65 * t84) - (t74 * t99 + t83 * t84) * t104) * MDP(25) + (-g(1) * (-t66 * t102 - t67 * t82) - g(2) * (-t64 * t102 - t65 * t82) - (t74 * t98 - t82 * t83) * t104) * MDP(26) + (MDP(4) - MDP(7) - MDP(16)) * (g(3) * t100 + t92) + (-t108 * t75 + t107 * t74 - MDP(3) - t80 * MDP(5) + MDP(6) * sin(pkin(11)) + (-t80 * pkin(3) - pkin(4) * t75 - qJ(5) * t74 - pkin(2)) * MDP(19)) * t87; t95 * t87; (-g(1) * (-t60 * pkin(4) + t61 * qJ(5)) - g(2) * (-t58 * pkin(4) + t59 * qJ(5)) - g(3) * (-t62 * pkin(4) + t63 * qJ(5))) * MDP(19) + (-MDP(25) * t82 - MDP(26) * t84 + t107) * (g(1) * t61 + g(2) * t59 + g(3) * t63) + t108 * t89; -t89 * MDP(19); (-g(1) * (t60 * t84 - t66 * t82) - g(2) * (t58 * t84 - t64 * t82) - g(3) * (t62 * t84 + t79 * t99)) * MDP(25) + (-g(1) * (-t60 * t82 - t66 * t84) - g(2) * (-t58 * t82 - t64 * t84) - g(3) * (-t62 * t82 + t79 * t98)) * MDP(26);];
taug  = t1;
