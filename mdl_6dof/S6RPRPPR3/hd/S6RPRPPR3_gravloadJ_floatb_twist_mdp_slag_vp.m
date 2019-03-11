% Calculate Gravitation load on the joints for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:24
% EndTime: 2019-03-09 02:45:25
% DurationCPUTime: 0.33s
% Computational Cost: add. (192->61), mult. (232->87), div. (0->0), fcn. (195->8), ass. (0->36)
t99 = MDP(10) + MDP(12) - MDP(17);
t98 = -MDP(11) + MDP(14) + MDP(16);
t71 = qJ(1) + pkin(9);
t65 = sin(t71);
t66 = cos(t71);
t57 = g(1) * t66 + g(2) * t65;
t73 = sin(qJ(3));
t97 = t57 * t73;
t67 = t73 * qJ(4);
t76 = cos(qJ(3));
t85 = t76 * pkin(3) + t67;
t95 = pkin(3) * t73;
t94 = g(1) * t65;
t90 = g(3) * t76;
t89 = t76 * pkin(4);
t88 = t66 * t76;
t72 = sin(qJ(6));
t87 = t72 * t73;
t75 = cos(qJ(6));
t86 = t73 * t75;
t84 = qJ(4) * t76;
t83 = -MDP(15) - MDP(19);
t74 = sin(qJ(1));
t82 = -t74 * pkin(1) + t66 * pkin(7);
t77 = cos(qJ(1));
t81 = t77 * pkin(1) + pkin(3) * t88 + t65 * pkin(7) + (pkin(2) + t67) * t66;
t80 = -g(2) * t66 + t94;
t78 = -pkin(2) - t85;
t60 = t66 * t84;
t58 = t65 * t84;
t54 = -t65 * t72 + t66 * t86;
t53 = -t65 * t75 - t66 * t87;
t52 = -t65 * t86 - t66 * t72;
t51 = t65 * t87 - t66 * t75;
t49 = -t90 + t97;
t1 = [(g(1) * t77 + g(2) * t74) * MDP(3) + (-g(1) * t82 - g(2) * t81 - t78 * t94) * MDP(15) + (-g(1) * (-t66 * qJ(5) + t82) - g(2) * (pkin(4) * t88 + t81) + (-g(1) * (t78 - t89) + g(2) * qJ(5)) * t65) * MDP(19) + (-g(1) * t52 - g(2) * t54) * MDP(25) + (-g(1) * t51 - g(2) * t53) * MDP(26) + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t74 - g(2) * t77) + (-MDP(13) + MDP(18)) * t57 + (t98 * t73 + t99 * t76) * t80; (-MDP(4) + t83) * g(3); (-g(1) * (-t66 * t95 + t60) - g(2) * (-t65 * t95 + t58) - g(3) * t85) * MDP(15) + (-g(1) * t60 - g(2) * t58 - g(3) * (t85 + t89) + (pkin(3) + pkin(4)) * t97) * MDP(19) + t99 * t49 + (-MDP(25) * t75 + MDP(26) * t72 - t98) * (g(3) * t73 + t57 * t76); t83 * t49; t80 * MDP(19); (-g(1) * t53 + g(2) * t51 - t72 * t90) * MDP(25) + (g(1) * t54 - g(2) * t52 - t75 * t90) * MDP(26);];
taug  = t1;
