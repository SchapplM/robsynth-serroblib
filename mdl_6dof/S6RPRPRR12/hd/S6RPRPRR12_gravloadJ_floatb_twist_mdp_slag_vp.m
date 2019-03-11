% Calculate Gravitation load on the joints for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR12_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRPRR12_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:19:58
% EndTime: 2019-03-09 04:19:59
% DurationCPUTime: 0.31s
% Computational Cost: add. (143->62), mult. (247->87), div. (0->0), fcn. (227->8), ass. (0->32)
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t76 = t70 * pkin(3) - t73 * qJ(4);
t71 = sin(qJ(1));
t74 = cos(qJ(1));
t94 = -g(1) * t71 + g(2) * t74;
t93 = -MDP(12) + MDP(15);
t92 = MDP(13) - MDP(16);
t88 = g(3) * t70;
t85 = t71 * t73;
t68 = qJ(5) + qJ(6);
t62 = sin(t68);
t84 = t74 * t62;
t63 = cos(t68);
t83 = t74 * t63;
t69 = sin(qJ(5));
t82 = t74 * t69;
t72 = cos(qJ(5));
t81 = t74 * t72;
t48 = t71 * t62 - t73 * t83;
t49 = -t71 * t63 - t73 * t84;
t50 = -t63 * t85 - t84;
t51 = -t62 * t85 + t83;
t80 = (-g(1) * t50 + g(2) * t48 - t63 * t88) * MDP(30) + (g(1) * t51 - g(2) * t49 + t62 * t88) * MDP(31);
t79 = t74 * pkin(1) + t71 * qJ(2);
t65 = t74 * qJ(2);
t57 = -t69 * t85 + t81;
t56 = -t72 * t85 - t82;
t55 = -t71 * t72 - t73 * t82;
t54 = t71 * t69 - t73 * t81;
t53 = -t73 * t94 - t88;
t1 = [(-g(1) * (-t71 * pkin(1) + t65) - g(2) * t79) * MDP(6) + (-g(1) * (t76 * t74 + t65) - g(2) * (t74 * pkin(7) + t79) + (-g(1) * (-pkin(1) - pkin(7)) - g(2) * t76) * t71) * MDP(17) + (-g(1) * t55 - g(2) * t57) * MDP(23) + (-g(1) * t54 - g(2) * t56) * MDP(24) + (-g(1) * t49 - g(2) * t51) * MDP(30) + (-g(1) * t48 - g(2) * t50) * MDP(31) - (MDP(2) - MDP(4) + MDP(14)) * t94 + (t93 * t70 - t92 * t73 + MDP(3) - MDP(5)) * (g(1) * t74 + g(2) * t71); -(-MDP(17) - MDP(6)) * t94; (g(3) * t76 + t94 * (pkin(3) * t73 + qJ(4) * t70)) * MDP(17) + t93 * t53 + (-MDP(23) * t69 - MDP(24) * t72 - MDP(30) * t62 - MDP(31) * t63 + t92) * (g(3) * t73 - t70 * t94); t53 * MDP(17); (-g(1) * t56 + g(2) * t54 - t72 * t88) * MDP(23) + (g(1) * t57 - g(2) * t55 + t69 * t88) * MDP(24) + t80; t80;];
taug  = t1;
