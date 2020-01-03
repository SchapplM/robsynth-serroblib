% Calculate Gravitation load on the joints for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:48:00
% EndTime: 2019-12-31 19:48:01
% DurationCPUTime: 0.40s
% Computational Cost: add. (132->63), mult. (245->90), div. (0->0), fcn. (217->8), ass. (0->35)
t76 = sin(qJ(2));
t78 = cos(qJ(2));
t86 = t78 * pkin(2) + t76 * qJ(3);
t81 = -pkin(1) - t86;
t100 = MDP(10) - MDP(13);
t99 = MDP(9) - MDP(12) + MDP(17);
t77 = sin(qJ(1));
t79 = cos(qJ(1));
t60 = g(1) * t79 + g(2) * t77;
t98 = t60 * t76;
t57 = g(3) * t76 + t60 * t78;
t97 = g(1) * t77;
t93 = g(3) * t78;
t92 = t76 * t79;
t73 = pkin(8) + qJ(5);
t65 = sin(t73);
t91 = t77 * t65;
t66 = cos(t73);
t90 = t77 * t66;
t74 = sin(pkin(8));
t89 = t77 * t74;
t75 = cos(pkin(8));
t88 = t77 * t75;
t85 = qJ(3) * t78;
t84 = qJ(4) * t78;
t83 = t77 * pkin(6) - t81 * t79;
t70 = t79 * pkin(6);
t63 = t79 * t85;
t61 = t77 * t85;
t56 = -t93 + t98;
t55 = t66 * t79 - t76 * t91;
t54 = t65 * t79 + t76 * t90;
t53 = t65 * t92 + t90;
t52 = t66 * t92 - t91;
t1 = [(-g(1) * t70 - g(2) * t83 - t81 * t97) * MDP(14) + (-g(1) * (t75 * t79 - t76 * t89) - g(2) * (t74 * t92 + t88)) * MDP(15) + (-g(1) * (-t74 * t79 - t76 * t88) - g(2) * (t75 * t92 - t89)) * MDP(16) + (-g(1) * (pkin(3) * t79 + t70) - g(2) * (t79 * t84 + t83) + (-g(1) * (t81 - t84) - g(2) * pkin(3)) * t77) * MDP(18) + (-g(1) * t55 - g(2) * t53) * MDP(24) + (g(1) * t54 - g(2) * t52) * MDP(25) + (MDP(3) - MDP(11)) * t60 + (-t100 * t76 + t99 * t78 + MDP(2)) * (-g(2) * t79 + t97); (-g(1) * (-pkin(2) * t92 + t63) - g(2) * (-pkin(2) * t76 * t77 + t61) - g(3) * t86) * MDP(14) + (-g(1) * t63 - g(2) * t61 - g(3) * (t84 + t86) + (pkin(2) + qJ(4)) * t98) * MDP(18) + t99 * t56 + (-MDP(15) * t74 - MDP(16) * t75 - MDP(24) * t65 - MDP(25) * t66 + t100) * t57; (-MDP(14) - MDP(18)) * t56; -t57 * MDP(18); (-g(1) * t52 - g(2) * t54 + t66 * t93) * MDP(24) + (g(1) * t53 - g(2) * t55 - t65 * t93) * MDP(25);];
taug = t1;
