% Calculate Gravitation load on the joints for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR14V3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR14V3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(1,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:10:07
% EndTime: 2019-04-12 15:10:10
% DurationCPUTime: 0.59s
% Computational Cost: add. (202->93), mult. (531->157), div. (0->0), fcn. (596->10), ass. (0->48)
t127 = MDP(9) + MDP(11);
t126 = -MDP(14) * qJ(3) + MDP(10) - MDP(13);
t86 = sin(qJ(5));
t88 = sin(qJ(2));
t115 = t86 * t88;
t87 = sin(qJ(4));
t94 = cos(qJ(1));
t105 = t94 * t87;
t92 = cos(qJ(4));
t93 = cos(qJ(2));
t107 = t92 * t93;
t89 = sin(qJ(1));
t76 = t107 * t89 - t105;
t91 = cos(qJ(5));
t65 = t115 * t89 + t76 * t91;
t106 = t92 * t94;
t109 = t89 * t87;
t75 = t109 * t93 + t106;
t85 = sin(qJ(6));
t90 = cos(qJ(6));
t125 = t65 * t85 - t75 * t90;
t124 = t65 * t90 + t75 * t85;
t123 = -g(1) * t94 - g(2) * t89;
t71 = -g(3) * t93 - t123 * t88;
t120 = g(3) * t88;
t116 = t85 * t91;
t114 = t87 * t88;
t113 = t87 * t93;
t112 = t88 * t91;
t111 = t88 * t92;
t110 = t88 * t94;
t108 = t90 * t91;
t103 = t85 * t114;
t102 = t90 * t114;
t101 = t88 * t105;
t64 = t112 * t89 - t76 * t86;
t74 = t111 * t91 - t86 * t93;
t98 = t111 * t86 + t91 * t93;
t79 = t106 * t93 + t109;
t78 = t105 * t93 - t89 * t92;
t77 = t107 * t91 + t115;
t70 = t74 * t94;
t69 = t74 * t89;
t68 = t110 * t86 + t79 * t91;
t67 = t110 * t91 - t79 * t86;
t63 = t68 * t90 + t78 * t85;
t62 = -t68 * t85 + t78 * t90;
t1 = [(g(1) * t76 - g(2) * t79) * MDP(20) + (-g(1) * t75 + g(2) * t78) * MDP(21) + (g(1) * t65 - g(2) * t68) * MDP(27) + (g(1) * t64 - g(2) * t67) * MDP(28) + (g(1) * t124 - g(2) * t63) * MDP(34) + (-g(1) * t125 - g(2) * t62) * MDP(35) - (MDP(3) - MDP(12)) * t123 + (-t126 * t88 + t127 * t93 + MDP(2)) * (g(1) * t89 - g(2) * t94); (g(1) * t70 + g(2) * t69 - g(3) * t77) * MDP(27) + (-g(3) * (-t107 * t86 + t112) + t123 * t98) * MDP(28) + (-g(1) * (-t101 * t85 - t70 * t90) - g(2) * (-t103 * t89 - t69 * t90) - g(3) * (t113 * t85 + t77 * t90)) * MDP(34) + (-g(1) * (-t101 * t90 + t70 * t85) - g(2) * (-t102 * t89 + t69 * t85) - g(3) * (t113 * t90 - t77 * t85)) * MDP(35) + t126 * (-t123 * t93 + t120) + (t92 * MDP(20) - t87 * MDP(21) + t127) * t71; -t71 * MDP(14); (g(1) * t79 + g(2) * t76 + g(3) * t111) * MDP(21) + (-g(1) * (-t108 * t78 + t79 * t85) - g(2) * (-t108 * t75 + t76 * t85) - (-t108 * t87 + t85 * t92) * t120) * MDP(34) + (-g(1) * (t116 * t78 + t79 * t90) - g(2) * (t116 * t75 + t76 * t90) - (t116 * t87 + t90 * t92) * t120) * MDP(35) + (MDP(27) * t91 - MDP(28) * t86 + MDP(20)) * (g(1) * t78 + g(2) * t75 + g(3) * t114); (g(1) * t68 + g(2) * t65 + g(3) * t74) * MDP(28) + (-MDP(34) * t90 + MDP(35) * t85 - MDP(27)) * (g(1) * t67 + g(2) * t64 - g(3) * t98); (-g(1) * t62 + g(2) * t125 - g(3) * (-t74 * t85 + t102)) * MDP(34) + (g(1) * t63 + g(2) * t124 - g(3) * (-t74 * t90 - t103)) * MDP(35);];
taug  = t1;
