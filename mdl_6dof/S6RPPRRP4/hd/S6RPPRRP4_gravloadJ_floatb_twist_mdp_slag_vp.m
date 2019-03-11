% Calculate Gravitation load on the joints for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:20
% EndTime: 2019-03-09 02:06:21
% DurationCPUTime: 0.40s
% Computational Cost: add. (228->68), mult. (480->98), div. (0->0), fcn. (551->8), ass. (0->32)
t112 = MDP(16) - MDP(25);
t111 = MDP(22) + MDP(24);
t110 = MDP(23) - MDP(26);
t100 = cos(pkin(9));
t104 = sin(qJ(1));
t105 = cos(qJ(1));
t99 = sin(pkin(9));
t70 = -t105 * t100 - t104 * t99;
t71 = -t104 * t100 + t105 * t99;
t93 = g(1) * t70 + g(2) * t71;
t83 = sin(qJ(4));
t107 = g(3) * t83;
t106 = t83 * pkin(8);
t82 = sin(qJ(5));
t85 = cos(qJ(4));
t103 = t82 * t85;
t84 = cos(qJ(5));
t102 = t84 * t85;
t101 = t105 * pkin(1) + t104 * qJ(2);
t98 = MDP(27) + MDP(9);
t97 = t105 * pkin(2) + t101;
t96 = -t104 * pkin(1) + t105 * qJ(2);
t92 = t85 * pkin(4) + pkin(3) + t106;
t91 = pkin(5) * t84 + qJ(6) * t82 + pkin(4);
t61 = t71 * t102 + t70 * t82;
t60 = t71 * t103 - t70 * t84;
t64 = -t70 * t103 - t71 * t84;
t90 = -g(1) * t64 + g(2) * t60 + t82 * t107;
t88 = -t104 * pkin(2) + t96;
t72 = g(1) * t104 - g(2) * t105;
t65 = -t70 * t102 + t71 * t82;
t1 = [(-g(1) * t96 - g(2) * t101) * MDP(6) + t93 * MDP(8) + (-g(1) * t88 - g(2) * t97) * MDP(9) + (-g(1) * (t61 * pkin(5) + t60 * qJ(6) + t88) - g(2) * (t65 * pkin(5) + t64 * qJ(6) + t97) + (-g(2) * pkin(7) - g(1) * t92) * t71 + (-g(1) * pkin(7) + g(2) * t92) * t70) * MDP(27) + t110 * (g(1) * t60 + g(2) * t64) + (MDP(3) - MDP(5)) * (g(1) * t105 + g(2) * t104) + (MDP(2) + MDP(4)) * t72 + t111 * (-g(1) * t61 - g(2) * t65) + (-t85 * MDP(15) + t112 * t83 - MDP(7)) * (g(1) * t71 - g(2) * t70); (-MDP(6) - t98) * t72; t98 * g(3); (-g(3) * (-t91 * t85 - t106) + t93 * (pkin(8) * t85 - t91 * t83)) * MDP(27) + t112 * (-t93 * t85 - t107) + (-t110 * t82 + t111 * t84 + MDP(15)) * (g(3) * t85 - t93 * t83); (-g(1) * (-t64 * pkin(5) + t65 * qJ(6)) - g(2) * (pkin(5) * t60 - qJ(6) * t61) - (pkin(5) * t82 - qJ(6) * t84) * t107) * MDP(27) - t110 * (-g(1) * t65 + g(2) * t61 + t84 * t107) - t111 * t90; t90 * MDP(27);];
taug  = t1;
