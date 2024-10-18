% Calculate Gravitation load on the joints for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:57
% EndTime: 2024-09-27 21:45:57
% DurationCPUTime: 0.16s
% Computational Cost: add. (439->60), mult. (231->88), div. (0->0), fcn. (208->18), ass. (0->46)
t106 = g(3) * sin(pkin(5));
t89 = cos(pkin(5));
t90 = sin(qJ(3));
t105 = t89 * t90;
t91 = cos(qJ(3));
t104 = t89 * t91;
t87 = qJ(3) + qJ(4);
t82 = pkin(5) + t87;
t77 = qJ(5) + t82;
t102 = sin(t77) / 0.2e1;
t71 = cos(t77) / 0.2e1;
t100 = pkin(5) - t87;
t98 = -qJ(5) + t100;
t75 = cos(t98);
t66 = t75 / 0.2e1 + t71;
t85 = qJ(5) + t87;
t78 = sin(t85);
t86 = pkin(10) + qJ(2);
t80 = sin(t86);
t81 = cos(t86);
t57 = -t81 * t66 + t80 * t78;
t58 = -t80 * t66 - t81 * t78;
t92 = sin(t98);
t65 = t102 - t92 / 0.2e1;
t79 = cos(t85);
t95 = t80 * t65 - t81 * t79;
t96 = t81 * t65 + t80 * t79;
t103 = (-g(1) * t58 + g(2) * t57 - g(3) * (t102 + t92 / 0.2e1)) * MDP(24) + (-g(1) * t95 + g(2) * t96 - g(3) * (t71 - t75 / 0.2e1)) * MDP(25);
t101 = sin(t82) / 0.2e1;
t74 = cos(t82) / 0.2e1;
t76 = cos(t100);
t68 = t76 / 0.2e1 + t74;
t83 = sin(t87);
t59 = -t81 * t68 + t80 * t83;
t60 = -t80 * t68 - t81 * t83;
t97 = sin(t100);
t67 = t101 - t97 / 0.2e1;
t84 = cos(t87);
t93 = t80 * t67 - t81 * t84;
t94 = t81 * t67 + t80 * t84;
t99 = (-g(1) * t60 + g(2) * t59 - g(3) * (t101 + t97 / 0.2e1)) * MDP(17) + (-g(1) * t93 + g(2) * t94 - g(3) * (t74 - t76 / 0.2e1)) * MDP(18) + t103;
t64 = -t80 * t105 + t81 * t91;
t63 = -t80 * t104 - t81 * t90;
t62 = -t81 * t105 - t80 * t91;
t61 = -t81 * t104 + t80 * t90;
t1 = [-g(3) * MDP(1); (g(1) * t80 - g(2) * t81) * MDP(3) + (g(1) * t81 + g(2) * t80) * MDP(4) + (-g(1) * t62 - g(2) * t64) * MDP(10) + (-g(1) * t61 - g(2) * t63) * MDP(11) + (g(1) * t94 + g(2) * t93) * MDP(17) + (-g(1) * t59 - g(2) * t60) * MDP(18) + (g(1) * t96 + g(2) * t95) * MDP(24) + (-g(1) * t57 - g(2) * t58) * MDP(25); (-g(1) * t63 + g(2) * t61 - t91 * t106) * MDP(10) + (g(1) * t64 - g(2) * t62 + t90 * t106) * MDP(11) + t99; t99; t103;];
taug = t1;
