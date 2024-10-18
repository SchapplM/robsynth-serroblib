% Calculate Gravitation load on the joints for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR13_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR13_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRRR13_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:16
% EndTime: 2024-09-27 17:32:16
% DurationCPUTime: 0.12s
% Computational Cost: add. (413->47), mult. (234->71), div. (0->0), fcn. (216->16), ass. (0->37)
t95 = g(3) * sin(pkin(5));
t80 = cos(pkin(5));
t81 = sin(qJ(4));
t94 = t80 * t81;
t83 = cos(qJ(4));
t93 = t80 * t83;
t77 = qJ(4) + qJ(5);
t71 = pkin(5) + t77;
t67 = cos(t71) / 0.2e1;
t90 = pkin(5) - t77;
t68 = cos(t90);
t64 = t68 / 0.2e1 + t67;
t78 = qJ(1) + qJ(2);
t76 = qJ(3) + t78;
t69 = sin(t76);
t70 = cos(t76);
t72 = sin(t77);
t53 = -t70 * t64 + t69 * t72;
t54 = -t69 * t64 - t70 * t72;
t89 = sin(t90);
t91 = sin(t71) / 0.2e1;
t63 = t91 - t89 / 0.2e1;
t74 = cos(t77);
t87 = t69 * t63 - t70 * t74;
t88 = t70 * t63 + t69 * t74;
t92 = (-g(1) * t54 + g(2) * t53 - g(3) * (t91 + t89 / 0.2e1)) * MDP(22) + (-g(1) * t87 + g(2) * t88 - g(3) * (t67 - t68 / 0.2e1)) * MDP(23);
t55 = t69 * t81 - t70 * t93;
t56 = -t69 * t83 - t70 * t94;
t57 = -t69 * t93 - t70 * t81;
t58 = -t69 * t94 + t70 * t83;
t86 = (g(1) * t69 - g(2) * t70) * MDP(8) + (g(1) * t70 + g(2) * t69) * MDP(9) + (-g(1) * t56 - g(2) * t58) * MDP(15) + (-g(1) * t55 - g(2) * t57) * MDP(16) + (g(1) * t88 + g(2) * t87) * MDP(22) + (-g(1) * t53 - g(2) * t54) * MDP(23);
t73 = sin(t78);
t75 = cos(t78);
t85 = (g(1) * t73 - g(2) * t75) * MDP(5) + (g(1) * t75 + g(2) * t73) * MDP(6) + t86;
t84 = cos(qJ(1));
t82 = sin(qJ(1));
t1 = [(g(1) * t82 - g(2) * t84) * MDP(2) + (g(1) * t84 + g(2) * t82) * MDP(3) + t85; t85; t86; (-g(1) * t57 + g(2) * t55 - t83 * t95) * MDP(15) + (g(1) * t58 - g(2) * t56 + t81 * t95) * MDP(16) + t92; t92;];
taug = t1;
