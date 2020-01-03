% Calculate Gravitation load on the joints for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:29
% EndTime: 2019-12-31 20:58:30
% DurationCPUTime: 0.32s
% Computational Cost: add. (227->59), mult. (267->69), div. (0->0), fcn. (216->6), ass. (0->31)
t70 = cos(qJ(2));
t67 = qJ(2) + qJ(3);
t64 = sin(t67);
t65 = cos(t67);
t81 = t65 * pkin(3) + t64 * qJ(4);
t78 = t70 * pkin(2) + t81;
t91 = pkin(1) + t78;
t71 = cos(qJ(1));
t90 = g(2) * t71;
t89 = MDP(16) + MDP(18) + MDP(22);
t88 = MDP(17) - MDP(20) - MDP(23);
t69 = sin(qJ(1));
t83 = g(1) * t71;
t53 = g(2) * t69 + t83;
t87 = t53 * t64;
t68 = sin(qJ(2));
t85 = pkin(2) * t68;
t84 = pkin(3) * t64;
t60 = t65 * pkin(4);
t80 = qJ(4) * t65;
t72 = -pkin(7) - pkin(6);
t79 = qJ(5) + t72;
t77 = t91 * t90;
t76 = -t84 - t85;
t52 = g(1) * t69 - t90;
t48 = -g(3) * t65 + t87;
t75 = t88 * (g(3) * t64 + t53 * t65) + t89 * t48;
t73 = (pkin(3) + pkin(4)) * t87;
t56 = t71 * t80;
t54 = t69 * t80;
t1 = [(t72 * t83 - t77 + (g(1) * t91 + g(2) * t72) * t69) * MDP(21) + (-t77 + (g(1) * t79 - g(2) * t60) * t71 + (-g(1) * (-t91 - t60) + g(2) * t79) * t69) * MDP(25) + (MDP(3) - MDP(19) + MDP(24)) * t53 + (-t68 * MDP(10) + t70 * MDP(9) - t88 * t64 + t89 * t65 + MDP(2)) * t52; (-g(3) * t70 + t53 * t68) * MDP(9) + (g(3) * t68 + t53 * t70) * MDP(10) + (-g(1) * (t71 * t76 + t56) - g(2) * (t69 * t76 + t54) - g(3) * t78) * MDP(21) + (-g(1) * (-t71 * t85 + t56) - g(2) * (-t69 * t85 + t54) - g(3) * (t60 + t78) + t73) * MDP(25) + t75; (-g(1) * (-t71 * t84 + t56) - g(2) * (-t69 * t84 + t54) - g(3) * t81) * MDP(21) + (-g(1) * t56 - g(2) * t54 - g(3) * (t60 + t81) + t73) * MDP(25) + t75; (-MDP(21) - MDP(25)) * t48; t52 * MDP(25);];
taug = t1;
