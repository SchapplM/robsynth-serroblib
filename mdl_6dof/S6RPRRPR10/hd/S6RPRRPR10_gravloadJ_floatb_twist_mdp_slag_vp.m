% Calculate Gravitation load on the joints for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR10_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:37:11
% EndTime: 2019-03-09 05:37:12
% DurationCPUTime: 0.55s
% Computational Cost: add. (186->78), mult. (446->109), div. (0->0), fcn. (471->8), ass. (0->34)
t134 = MDP(13) - MDP(22);
t125 = MDP(19) + MDP(21);
t124 = MDP(20) - MDP(23);
t89 = sin(qJ(6));
t90 = sin(qJ(4));
t93 = cos(qJ(6));
t94 = cos(qJ(4));
t100 = t89 * t90 + t93 * t94;
t101 = t89 * t94 - t90 * t93;
t92 = sin(qJ(1));
t107 = t92 * t94;
t91 = sin(qJ(3));
t96 = cos(qJ(1));
t109 = t91 * t96;
t78 = t90 * t109 + t107;
t106 = t96 * t94;
t108 = t92 * t90;
t79 = t91 * t106 - t108;
t102 = t78 * t89 + t79 * t93;
t104 = -t78 * t93 + t79 * t89;
t95 = cos(qJ(3));
t112 = g(3) * t95;
t76 = t91 * t108 - t106;
t77 = t91 * t107 + t90 * t96;
t67 = t76 * t93 - t77 * t89;
t68 = t76 * t89 + t77 * t93;
t132 = (-g(1) * t67 - g(2) * t104 + t101 * t112) * MDP(30) - (-g(1) * t68 + g(2) * t102 - t100 * t112) * MDP(31);
t127 = -g(1) * t92 + g(2) * t96;
t111 = t95 * pkin(8);
t105 = t96 * pkin(1) + t92 * qJ(2);
t99 = pkin(4) * t94 + qJ(5) * t90 + pkin(3);
t66 = g(1) * t76 - g(2) * t78 + t90 * t112;
t86 = t96 * qJ(2);
t1 = [(-g(1) * (-t92 * pkin(1) + t86) - g(2) * t105) * MDP(6) + (-g(1) * (pkin(3) * t109 + t79 * pkin(4) + t78 * qJ(5) - t96 * t111 + t86) - g(2) * (t77 * pkin(4) + pkin(7) * t96 + t76 * qJ(5) + t105) + (-g(1) * (-pkin(1) - pkin(7)) - g(2) * (pkin(3) * t91 - t111)) * t92) * MDP(24) + (-g(1) * t102 - g(2) * t68) * MDP(30) + (g(1) * t104 - g(2) * t67) * MDP(31) - (MDP(2) - MDP(4)) * t127 + t125 * (-g(1) * t79 - g(2) * t77) + t124 * (g(1) * t78 + g(2) * t76) + (-t91 * MDP(12) - t134 * t95 + MDP(3) - MDP(5)) * (g(1) * t96 + g(2) * t92); -(-MDP(24) - MDP(6)) * t127; (-g(3) * (-t99 * t91 + t111) + t127 * (pkin(8) * t91 + t99 * t95)) * MDP(24) + t134 * (-t127 * t91 + t112) + (t100 * MDP(30) - t101 * MDP(31) - t124 * t90 + t125 * t94 + MDP(12)) * (g(3) * t91 + t127 * t95); (-g(1) * (-pkin(4) * t76 + qJ(5) * t77) - g(2) * (pkin(4) * t78 - qJ(5) * t79) - (-pkin(4) * t90 + qJ(5) * t94) * t112) * MDP(24) + t124 * (g(1) * t77 - g(2) * t79 + t94 * t112) + t125 * t66 - t132; -t66 * MDP(24); t132;];
taug  = t1;
