% Calculate Gravitation load on the joints for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:52
% EndTime: 2019-12-31 17:30:53
% DurationCPUTime: 0.30s
% Computational Cost: add. (138->62), mult. (354->115), div. (0->0), fcn. (424->10), ass. (0->34)
t61 = sin(qJ(2));
t62 = sin(qJ(1));
t65 = cos(qJ(2));
t66 = cos(qJ(1));
t72 = cos(pkin(4));
t70 = t66 * t72;
t53 = t61 * t70 + t62 * t65;
t60 = sin(qJ(3));
t64 = cos(qJ(3));
t58 = sin(pkin(4));
t76 = t58 * t66;
t46 = t53 * t64 - t60 * t76;
t52 = t62 * t61 - t65 * t70;
t59 = sin(qJ(4));
t63 = cos(qJ(4));
t84 = t46 * t59 - t52 * t63;
t83 = t46 * t63 + t52 * t59;
t82 = g(3) * t58;
t79 = t58 * t61;
t78 = t58 * t64;
t77 = t58 * t65;
t75 = t59 * t64;
t74 = t63 * t64;
t73 = t64 * t65;
t71 = t62 * t72;
t69 = t53 * t60 + t64 * t76;
t55 = -t61 * t71 + t66 * t65;
t54 = t66 * t61 + t65 * t71;
t51 = t72 * t60 + t61 * t78;
t49 = t62 * t58 * t60 + t55 * t64;
t48 = -t55 * t60 + t62 * t78;
t44 = t49 * t63 + t54 * t59;
t43 = -t49 * t59 + t54 * t63;
t1 = [(g(1) * t62 - g(2) * t66) * MDP(2) + (g(1) * t66 + g(2) * t62) * MDP(3) + (g(1) * t53 - g(2) * t55) * MDP(9) + (-g(1) * t52 + g(2) * t54) * MDP(10) + (g(1) * t46 - g(2) * t49) * MDP(16) + (-g(1) * t69 - g(2) * t48) * MDP(17) + (g(1) * t83 - g(2) * t44) * MDP(23) + (-g(1) * t84 - g(2) * t43) * MDP(24); (g(1) * t55 + g(2) * t53 + g(3) * t79) * MDP(10) + (-g(1) * (-t54 * t74 + t55 * t59) - g(2) * (-t52 * t74 + t53 * t59) - (t59 * t61 + t63 * t73) * t82) * MDP(23) + (-g(1) * (t54 * t75 + t55 * t63) - g(2) * (t52 * t75 + t53 * t63) - (-t59 * t73 + t61 * t63) * t82) * MDP(24) + (-t64 * MDP(16) + MDP(17) * t60 - MDP(9)) * (-g(1) * t54 - g(2) * t52 + g(3) * t77); (g(1) * t49 + g(2) * t46 + g(3) * t51) * MDP(17) + (-MDP(23) * t63 + MDP(24) * t59 - MDP(16)) * (g(1) * t48 - g(2) * t69 + g(3) * (-t60 * t79 + t72 * t64)); (-g(1) * t43 + g(2) * t84 - g(3) * (-t51 * t59 - t63 * t77)) * MDP(23) + (g(1) * t44 + g(2) * t83 - g(3) * (-t51 * t63 + t59 * t77)) * MDP(24);];
taug = t1;
