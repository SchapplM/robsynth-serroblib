% Calculate Gravitation load on the joints for
% S6RPRPRR4
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
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:27
% EndTime: 2019-03-09 03:46:28
% DurationCPUTime: 0.23s
% Computational Cost: add. (215->59), mult. (235->90), div. (0->0), fcn. (217->10), ass. (0->35)
t63 = qJ(1) + pkin(10);
t58 = sin(t63);
t59 = cos(t63);
t78 = g(1) * t59 + g(2) * t58;
t90 = MDP(10) - MDP(13);
t89 = MDP(11) - MDP(14);
t69 = cos(qJ(3));
t84 = g(3) * t69;
t64 = qJ(5) + qJ(6);
t60 = sin(t64);
t66 = sin(qJ(3));
t83 = t60 * t66;
t61 = cos(t64);
t82 = t61 * t66;
t65 = sin(qJ(5));
t81 = t65 * t66;
t68 = cos(qJ(5));
t80 = t66 * t68;
t46 = -t58 * t60 + t59 * t82;
t47 = t58 * t61 + t59 * t83;
t48 = t58 * t82 + t59 * t60;
t49 = -t58 * t83 + t59 * t61;
t79 = (-g(1) * t46 - g(2) * t48 + t61 * t84) * MDP(28) + (g(1) * t47 - g(2) * t49 - t60 * t84) * MDP(29);
t67 = sin(qJ(1));
t70 = cos(qJ(1));
t76 = g(1) * t67 - g(2) * t70;
t75 = t69 * pkin(3) + t66 * qJ(4);
t73 = pkin(2) + t75;
t72 = t76 * pkin(1);
t55 = -t58 * t81 + t59 * t68;
t54 = t58 * t80 + t59 * t65;
t53 = t58 * t68 + t59 * t81;
t52 = -t58 * t65 + t59 * t80;
t50 = t78 * t66 - t84;
t1 = [t76 * MDP(2) + (g(1) * t70 + g(2) * t67) * MDP(3) + MDP(4) * t72 - t78 * MDP(12) + (t72 + (-g(1) * pkin(7) - g(2) * t73) * t59 + (-g(2) * pkin(7) + g(1) * t73) * t58) * MDP(15) + (-g(1) * t55 - g(2) * t53) * MDP(21) + (g(1) * t54 - g(2) * t52) * MDP(22) + (-g(1) * t49 - g(2) * t47) * MDP(28) + (g(1) * t48 - g(2) * t46) * MDP(29) + (-t89 * t66 + t90 * t69) * (g(1) * t58 - g(2) * t59); (-MDP(15) - MDP(4)) * g(3); (-g(3) * t75 + t78 * (pkin(3) * t66 - qJ(4) * t69)) * MDP(15) + t90 * t50 + (-MDP(21) * t65 - MDP(22) * t68 - MDP(28) * t60 - MDP(29) * t61 + t89) * (g(3) * t66 + t78 * t69); -t50 * MDP(15); (-g(1) * t52 - g(2) * t54 + t68 * t84) * MDP(21) + (g(1) * t53 - g(2) * t55 - t65 * t84) * MDP(22) + t79; t79;];
taug  = t1;
