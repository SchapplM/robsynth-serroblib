% Calculate Gravitation load on the joints for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:16:25
% EndTime: 2019-03-09 03:16:27
% DurationCPUTime: 0.65s
% Computational Cost: add. (363->86), mult. (381->119), div. (0->0), fcn. (355->10), ass. (0->42)
t86 = sin(qJ(1));
t87 = cos(qJ(1));
t67 = g(1) * t87 + g(2) * t86;
t110 = MDP(24) + MDP(26);
t109 = MDP(25) - MDP(28);
t79 = pkin(9) + qJ(3);
t74 = sin(t79);
t76 = cos(t79);
t58 = -g(3) * t76 + t67 * t74;
t108 = MDP(14) - MDP(17) - MDP(27);
t105 = g(3) * t74;
t78 = pkin(10) + qJ(5);
t75 = cos(t78);
t103 = t75 * t87;
t80 = sin(pkin(10));
t102 = t80 * t87;
t82 = cos(pkin(10));
t101 = t82 * t87;
t73 = sin(t78);
t100 = t86 * t73;
t99 = t86 * t75;
t98 = t86 * t80;
t97 = t86 * t82;
t96 = t87 * t73;
t95 = -MDP(18) - MDP(29);
t85 = -pkin(7) - qJ(2);
t94 = pkin(4) * t80 - t85;
t66 = g(1) * t86 - g(2) * t87;
t92 = pkin(3) * t76 + qJ(4) * t74;
t71 = pkin(4) * t82 + pkin(3);
t84 = -pkin(8) - qJ(4);
t90 = t71 * t76 - t74 * t84;
t89 = pkin(5) * t75 + qJ(6) * t73 + t71;
t60 = t76 * t100 + t103;
t62 = t76 * t96 - t99;
t54 = g(1) * t62 + g(2) * t60 + t73 * t105;
t83 = cos(pkin(9));
t72 = pkin(2) * t83 + pkin(1);
t65 = t87 * t72;
t63 = t76 * t103 + t100;
t61 = t76 * t99 - t96;
t1 = [(-g(1) * (-t86 * pkin(1) + qJ(2) * t87) - g(2) * (pkin(1) * t87 + t86 * qJ(2))) * MDP(7) + (-g(1) * (-t76 * t97 + t102) - g(2) * (t76 * t101 + t98)) * MDP(15) + (-g(1) * (t76 * t98 + t101) - g(2) * (-t76 * t102 + t97)) * MDP(16) + (-g(2) * t65 + (g(1) * t85 - g(2) * t92) * t87 + (-g(1) * (-t72 - t92) + g(2) * t85) * t86) * MDP(18) + (-g(1) * (-t61 * pkin(5) - t60 * qJ(6)) - g(2) * (t63 * pkin(5) + t62 * qJ(6) + t65) + (-g(1) * t94 - g(2) * t90) * t87 + (-g(1) * (-t72 - t90) - g(2) * t94) * t86) * MDP(29) - t109 * (g(1) * t60 - g(2) * t62) + (MDP(3) - MDP(6)) * t67 + t110 * (g(1) * t61 - g(2) * t63) + (t76 * MDP(13) + MDP(4) * t83 - MDP(5) * sin(pkin(9)) - t108 * t74 + MDP(2)) * t66; (-MDP(7) + t95) * t66; (-g(3) * t92 + t67 * (pkin(3) * t74 - qJ(4) * t76)) * MDP(18) + ((-g(3) * t89 + t67 * t84) * t76 + (g(3) * t84 + t67 * t89) * t74) * MDP(29) + t108 * (t67 * t76 + t105) + (t82 * MDP(15) - MDP(16) * t80 - t109 * t73 + t110 * t75 + MDP(13)) * t58; t95 * t58; (-g(1) * (-pkin(5) * t62 + qJ(6) * t63) - g(2) * (-pkin(5) * t60 + qJ(6) * t61) - (-pkin(5) * t73 + qJ(6) * t75) * t105) * MDP(29) + t109 * (g(1) * t63 + g(2) * t61 + t75 * t105) + t110 * t54; -t54 * MDP(29);];
taug  = t1;
