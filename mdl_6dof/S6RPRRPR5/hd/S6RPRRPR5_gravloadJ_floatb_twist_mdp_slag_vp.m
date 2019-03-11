% Calculate Gravitation load on the joints for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:34
% EndTime: 2019-03-09 05:13:35
% DurationCPUTime: 0.31s
% Computational Cost: add. (281->63), mult. (259->85), div. (0->0), fcn. (220->10), ass. (0->36)
t77 = pkin(10) + qJ(3);
t73 = cos(t77);
t74 = qJ(4) + t77;
t70 = sin(t74);
t71 = cos(t74);
t89 = t71 * pkin(4) + t70 * qJ(5);
t99 = pkin(3) * t73 + t89;
t98 = MDP(20) - MDP(23);
t97 = MDP(21) - MDP(24);
t96 = pkin(4) * t70;
t94 = g(3) * t71;
t80 = sin(qJ(6));
t81 = sin(qJ(1));
t93 = t81 * t80;
t82 = cos(qJ(6));
t92 = t81 * t82;
t83 = cos(qJ(1));
t91 = t83 * t80;
t90 = t83 * t82;
t88 = qJ(5) * t71;
t72 = sin(t77);
t87 = -pkin(3) * t72 - t96;
t65 = g(1) * t83 + g(2) * t81;
t64 = g(1) * t81 - g(2) * t83;
t53 = t65 * t70 - t94;
t86 = t98 * t53 + (-MDP(31) * t80 - MDP(32) * t82 + t97) * (g(3) * t70 + t65 * t71);
t79 = cos(pkin(10));
t85 = t79 * pkin(2) + pkin(1) + t99;
t76 = -pkin(8) - pkin(7) - qJ(2);
t63 = t83 * t88;
t62 = t81 * t88;
t60 = -t70 * t93 + t90;
t59 = t70 * t92 + t91;
t58 = t70 * t91 + t92;
t57 = t70 * t90 - t93;
t1 = [(-g(1) * (-t81 * pkin(1) + t83 * qJ(2)) - g(2) * (t83 * pkin(1) + t81 * qJ(2))) * MDP(7) + ((g(1) * t76 - g(2) * t85) * t83 + (g(1) * t85 + g(2) * t76) * t81) * MDP(25) + (-g(1) * t60 - g(2) * t58) * MDP(31) + (g(1) * t59 - g(2) * t57) * MDP(32) + (MDP(3) - MDP(6) - MDP(22)) * t65 + (t98 * t71 - t97 * t70 + t73 * MDP(13) - t72 * MDP(14) + MDP(4) * t79 - MDP(5) * sin(pkin(10)) + MDP(2)) * t64; (-MDP(25) - MDP(7)) * t64; (-g(3) * t73 + t65 * t72) * MDP(13) + (g(3) * t72 + t65 * t73) * MDP(14) + (-g(1) * (t83 * t87 + t63) - g(2) * (t81 * t87 + t62) - g(3) * t99) * MDP(25) + t86; (-g(1) * (-t83 * t96 + t63) - g(2) * (-t81 * t96 + t62) - g(3) * t89) * MDP(25) + t86; -t53 * MDP(25); (-g(1) * t57 - g(2) * t59 + t82 * t94) * MDP(31) + (g(1) * t58 - g(2) * t60 - t80 * t94) * MDP(32);];
taug  = t1;
