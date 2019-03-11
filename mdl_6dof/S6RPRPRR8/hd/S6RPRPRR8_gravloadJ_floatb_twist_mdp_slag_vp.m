% Calculate Gravitation load on the joints for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR8_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR8_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:31
% EndTime: 2019-03-09 03:59:32
% DurationCPUTime: 0.19s
% Computational Cost: add. (167->60), mult. (213->82), div. (0->0), fcn. (197->10), ass. (0->40)
t93 = MDP(15) * pkin(3) + MDP(12);
t66 = qJ(5) + qJ(6);
t59 = sin(t66);
t60 = cos(t66);
t68 = sin(qJ(5));
t71 = cos(qJ(5));
t92 = -MDP(21) * t71 + MDP(22) * t68 - MDP(28) * t60 + MDP(29) * t59;
t69 = sin(qJ(3));
t91 = pkin(3) * t69;
t65 = qJ(3) + pkin(10);
t58 = cos(t65);
t90 = g(3) * t58;
t73 = cos(qJ(1));
t89 = t59 * t73;
t88 = t60 * t73;
t87 = t68 * t73;
t70 = sin(qJ(1));
t86 = t70 * t59;
t85 = t70 * t60;
t84 = t70 * t68;
t83 = t70 * t71;
t82 = t71 * t73;
t57 = sin(t65);
t47 = -t57 * t86 + t88;
t48 = t57 * t85 + t89;
t49 = t57 * t89 + t85;
t50 = t57 * t88 - t86;
t81 = (-g(1) * t47 - g(2) * t49 + t59 * t90) * MDP(28) + (g(1) * t48 - g(2) * t50 + t60 * t90) * MDP(29);
t80 = t73 * pkin(1) + t70 * qJ(2);
t72 = cos(qJ(3));
t74 = t72 * MDP(13);
t56 = g(1) * t73 + g(2) * t70;
t55 = g(1) * t70 - g(2) * t73;
t67 = -qJ(4) - pkin(7);
t62 = t73 * qJ(2);
t54 = t57 * t82 - t84;
t53 = t57 * t87 + t83;
t52 = t57 * t83 + t87;
t51 = -t57 * t84 + t82;
t1 = [(-g(1) * (-t70 * pkin(1) + t62) - g(2) * t80) * MDP(6) + (-g(1) * (t73 * t91 + t62 + (-pkin(1) + t67) * t70) - g(2) * (-t73 * t67 + t70 * t91 + t80)) * MDP(15) + (-g(1) * t54 - g(2) * t52) * MDP(21) + (g(1) * t53 - g(2) * t51) * MDP(22) + (-g(1) * t50 - g(2) * t48) * MDP(28) + (g(1) * t49 - g(2) * t47) * MDP(29) + (MDP(2) - MDP(4) + MDP(14)) * t55 + (-t69 * MDP(12) + MDP(3) - MDP(5) - t74) * t56; (-MDP(15) - MDP(6)) * t55; (-t92 * t57 + t93 * t69 + t74) * g(3) + (MDP(13) * t69 + t92 * t58 - t93 * t72) * t55; -t56 * MDP(15); (-g(1) * t51 - g(2) * t53 + t68 * t90) * MDP(21) + (g(1) * t52 - g(2) * t54 + t71 * t90) * MDP(22) + t81; t81;];
taug  = t1;
