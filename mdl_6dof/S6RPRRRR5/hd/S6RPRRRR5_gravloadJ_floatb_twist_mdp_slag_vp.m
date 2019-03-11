% Calculate Gravitation load on the joints for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:09:50
% EndTime: 2019-03-09 07:09:51
% DurationCPUTime: 0.22s
% Computational Cost: add. (288->56), mult. (259->80), div. (0->0), fcn. (242->12), ass. (0->35)
t70 = pkin(11) + qJ(3);
t67 = qJ(4) + t70;
t63 = sin(t67);
t92 = g(3) * t63;
t71 = qJ(5) + qJ(6);
t68 = sin(t71);
t77 = cos(qJ(1));
t90 = t68 * t77;
t69 = cos(t71);
t89 = t69 * t77;
t74 = sin(qJ(5));
t88 = t74 * t77;
t75 = sin(qJ(1));
t87 = t75 * t68;
t86 = t75 * t69;
t85 = t75 * t74;
t76 = cos(qJ(5));
t84 = t75 * t76;
t83 = t76 * t77;
t64 = cos(t67);
t53 = t64 * t87 + t89;
t54 = -t64 * t86 + t90;
t55 = -t64 * t90 + t86;
t56 = t64 * t89 + t87;
t82 = (-g(1) * t55 + g(2) * t53 + t68 * t92) * MDP(34) + (g(1) * t56 - g(2) * t54 + t69 * t92) * MDP(35);
t81 = g(1) * t77 + g(2) * t75;
t61 = g(1) * t75 - g(2) * t77;
t80 = (t81 * t64 + t92) * MDP(21) + (t76 * MDP(27) - t74 * MDP(28) + t69 * MDP(34) - t68 * MDP(35) + MDP(20)) * (-g(3) * t64 + t81 * t63);
t66 = cos(t70);
t65 = sin(t70);
t60 = t64 * t83 + t85;
t59 = -t64 * t88 + t84;
t58 = -t64 * t84 + t88;
t57 = t64 * t85 + t83;
t1 = [(-g(1) * (-t75 * pkin(1) + qJ(2) * t77) - g(2) * (pkin(1) * t77 + t75 * qJ(2))) * MDP(7) + (-g(1) * t58 - g(2) * t60) * MDP(27) + (-g(1) * t57 - g(2) * t59) * MDP(28) + (-g(1) * t54 - g(2) * t56) * MDP(34) + (-g(1) * t53 - g(2) * t55) * MDP(35) + (MDP(3) - MDP(6)) * t81 + (t66 * MDP(13) - t65 * MDP(14) + MDP(20) * t64 - MDP(21) * t63 + MDP(4) * cos(pkin(11)) - MDP(5) * sin(pkin(11)) + MDP(2)) * t61; -t61 * MDP(7); (-g(3) * t66 + t81 * t65) * MDP(13) + (g(3) * t65 + t81 * t66) * MDP(14) + t80; t80; (-g(1) * t59 + g(2) * t57 + t74 * t92) * MDP(27) + (g(1) * t60 - g(2) * t58 + t76 * t92) * MDP(28) + t82; t82;];
taug  = t1;
