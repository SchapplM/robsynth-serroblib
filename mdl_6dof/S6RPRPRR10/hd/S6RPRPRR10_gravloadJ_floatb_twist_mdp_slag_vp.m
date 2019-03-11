% Calculate Gravitation load on the joints for
% S6RPRPRR10
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
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:09:19
% EndTime: 2019-03-09 04:09:20
% DurationCPUTime: 0.35s
% Computational Cost: add. (208->71), mult. (269->102), div. (0->0), fcn. (252->10), ass. (0->41)
t82 = sin(qJ(3));
t84 = cos(qJ(3));
t87 = t82 * pkin(3) - t84 * qJ(4);
t83 = sin(qJ(1));
t85 = cos(qJ(1));
t110 = -g(1) * t83 + g(2) * t85;
t109 = MDP(13) - MDP(16);
t66 = -g(3) * t82 - t110 * t84;
t105 = g(3) * t84;
t79 = pkin(10) + qJ(5);
t74 = qJ(6) + t79;
t70 = sin(t74);
t103 = t83 * t70;
t71 = cos(t74);
t102 = t83 * t71;
t72 = sin(t79);
t101 = t83 * t72;
t73 = cos(t79);
t100 = t83 * t73;
t80 = sin(pkin(10));
t99 = t83 * t80;
t81 = cos(pkin(10));
t98 = t83 * t81;
t97 = t85 * t70;
t96 = t85 * t71;
t95 = t85 * t72;
t94 = t85 * t73;
t93 = t85 * t80;
t92 = t85 * t81;
t57 = -t82 * t103 + t96;
t58 = t82 * t102 + t97;
t59 = t82 * t97 + t102;
t60 = t82 * t96 - t103;
t91 = (-g(1) * t57 - g(2) * t59 + t70 * t105) * MDP(30) + (g(1) * t58 - g(2) * t60 + t71 * t105) * MDP(31);
t90 = t85 * pkin(1) + t83 * qJ(2);
t76 = t85 * qJ(2);
t64 = t82 * t94 - t101;
t63 = t82 * t95 + t100;
t62 = t82 * t100 + t95;
t61 = -t82 * t101 + t94;
t1 = [(-g(1) * (-t83 * pkin(1) + t76) - g(2) * t90) * MDP(6) + (-g(1) * (t82 * t92 - t99) - g(2) * (t82 * t98 + t93)) * MDP(14) + (-g(1) * (-t82 * t93 - t98) - g(2) * (-t82 * t99 + t92)) * MDP(15) + (-g(1) * (t87 * t85 + t76) - g(2) * (t85 * pkin(7) + t90) + (-g(1) * (-pkin(1) - pkin(7)) - g(2) * t87) * t83) * MDP(17) + (-g(1) * t64 - g(2) * t62) * MDP(23) + (g(1) * t63 - g(2) * t61) * MDP(24) + (-g(1) * t60 - g(2) * t58) * MDP(30) + (g(1) * t59 - g(2) * t57) * MDP(31) - (MDP(2) - MDP(4)) * t110 + (-t82 * MDP(12) - t109 * t84 + MDP(3) - MDP(5)) * (g(1) * t85 + g(2) * t83); -(-MDP(17) - MDP(6)) * t110; (g(3) * t87 + t110 * (pkin(3) * t84 + qJ(4) * t82)) * MDP(17) + t109 * (-t110 * t82 + t105) + (-t81 * MDP(14) + t80 * MDP(15) - t73 * MDP(23) + t72 * MDP(24) - t71 * MDP(30) + t70 * MDP(31) - MDP(12)) * t66; t66 * MDP(17); (-g(1) * t61 - g(2) * t63 + t72 * t105) * MDP(23) + (g(1) * t62 - g(2) * t64 + t73 * t105) * MDP(24) + t91; t91;];
taug  = t1;
