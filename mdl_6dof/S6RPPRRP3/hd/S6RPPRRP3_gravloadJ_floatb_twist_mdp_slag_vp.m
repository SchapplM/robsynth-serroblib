% Calculate Gravitation load on the joints for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPRRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:48
% EndTime: 2019-03-09 02:03:50
% DurationCPUTime: 0.37s
% Computational Cost: add. (229->65), mult. (284->91), div. (0->0), fcn. (265->8), ass. (0->28)
t66 = qJ(1) + pkin(9);
t63 = sin(t66);
t64 = cos(t66);
t98 = -g(1) * t63 + g(2) * t64;
t68 = sin(qJ(4));
t71 = cos(qJ(4));
t85 = t71 * pkin(8);
t97 = t68 * pkin(4) - t85;
t96 = MDP(14) - MDP(23);
t95 = MDP(20) + MDP(22);
t94 = MDP(21) - MDP(24);
t87 = g(3) * t71;
t67 = sin(qJ(5));
t84 = t67 * t68;
t70 = cos(qJ(5));
t83 = t68 * t70;
t82 = MDP(25) + MDP(7);
t72 = cos(qJ(1));
t81 = t72 * pkin(1) + t64 * pkin(2) + t63 * qJ(3);
t69 = sin(qJ(1));
t80 = -t69 * pkin(1) + t64 * qJ(3);
t75 = pkin(5) * t70 + qJ(6) * t67 + pkin(4);
t53 = t63 * t84 - t64 * t70;
t55 = t63 * t70 + t64 * t84;
t48 = g(1) * t53 - g(2) * t55 + t67 * t87;
t56 = -t63 * t67 + t64 * t83;
t54 = t63 * t83 + t64 * t67;
t1 = [(g(1) * t72 + g(2) * t69) * MDP(3) + t98 * MDP(5) + (-g(1) * (-t63 * pkin(2) + t80) - g(2) * t81) * MDP(7) + (-g(1) * (t56 * pkin(5) + t55 * qJ(6) + t97 * t64 + t80) - g(2) * (t54 * pkin(5) + t64 * pkin(7) + t53 * qJ(6) + t81) + (-g(1) * (-pkin(2) - pkin(7)) - g(2) * t97) * t63) * MDP(25) + t94 * (g(1) * t55 + g(2) * t53) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t69 - g(2) * t72) + t95 * (-g(1) * t56 - g(2) * t54) + (-t68 * MDP(13) - t96 * t71 - MDP(6)) * (g(1) * t64 + g(2) * t63); (-MDP(4) - t82) * g(3); t82 * t98; (-g(3) * (-t68 * t75 + t85) + t98 * (pkin(8) * t68 + t71 * t75)) * MDP(25) + t96 * (-t68 * t98 + t87) + (t94 * t67 - t95 * t70 - MDP(13)) * (-g(3) * t68 - t98 * t71); (-g(1) * (-t53 * pkin(5) + t54 * qJ(6)) - g(2) * (t55 * pkin(5) - t56 * qJ(6)) - (-pkin(5) * t67 + qJ(6) * t70) * t87) * MDP(25) + t94 * (g(1) * t54 - g(2) * t56 + t70 * t87) + t95 * t48; -t48 * MDP(25);];
taug  = t1;
