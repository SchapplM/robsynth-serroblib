% Calculate Gravitation load on the joints for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:49:36
% EndTime: 2019-03-09 03:49:37
% DurationCPUTime: 0.32s
% Computational Cost: add. (248->59), mult. (330->86), div. (0->0), fcn. (335->10), ass. (0->31)
t71 = pkin(10) + qJ(3);
t68 = sin(t71);
t69 = cos(t71);
t76 = sin(qJ(5));
t92 = cos(qJ(5));
t100 = t68 * t92 - t69 * t76;
t77 = sin(qJ(1));
t53 = t100 * t77;
t82 = t68 * t76 + t69 * t92;
t54 = t82 * t77;
t79 = cos(qJ(1));
t55 = t100 * t79;
t56 = t82 * t79;
t75 = sin(qJ(6));
t78 = cos(qJ(6));
t93 = g(3) * t100;
t101 = (MDP(31) * t78 - MDP(32) * t75 + MDP(24)) * (-g(1) * t55 - g(2) * t53 + g(3) * t82) + (g(1) * t56 + g(2) * t54 + t93) * MDP(25);
t64 = g(1) * t79 + g(2) * t77;
t98 = MDP(13) + MDP(15);
t97 = MDP(14) - MDP(17);
t63 = g(1) * t77 - g(2) * t79;
t87 = t69 * pkin(3) + t68 * qJ(4);
t85 = t54 * t78 + t79 * t75;
t84 = t54 * t75 - t79 * t78;
t73 = cos(pkin(10));
t83 = t73 * pkin(2) + pkin(1) + t87;
t74 = -pkin(7) - qJ(2);
t51 = -g(3) * t69 + t64 * t68;
t50 = t56 * t78 - t77 * t75;
t49 = -t56 * t75 - t77 * t78;
t1 = [(-g(1) * (-t77 * pkin(1) + t79 * qJ(2)) - g(2) * (t79 * pkin(1) + t77 * qJ(2))) * MDP(7) + ((g(1) * t74 - g(2) * t83) * t79 + (g(1) * t83 + g(2) * t74) * t77) * MDP(18) + (g(1) * t54 - g(2) * t56) * MDP(24) + (g(1) * t53 - g(2) * t55) * MDP(25) + (g(1) * t85 - g(2) * t50) * MDP(31) + (-g(1) * t84 - g(2) * t49) * MDP(32) + (MDP(3) - MDP(6) - MDP(16)) * t64 + (MDP(4) * t73 - MDP(5) * sin(pkin(10)) - t97 * t68 + t98 * t69 + MDP(2)) * t63; (-MDP(18) - MDP(7)) * t63; (-g(3) * t87 + t64 * (pkin(3) * t68 - qJ(4) * t69)) * MDP(18) + t97 * (g(3) * t68 + t64 * t69) + t98 * t51 - t101; -t51 * MDP(18); t101; (-g(1) * t49 + g(2) * t84 + t75 * t93) * MDP(31) + (g(1) * t50 + g(2) * t85 + t78 * t93) * MDP(32);];
taug  = t1;
