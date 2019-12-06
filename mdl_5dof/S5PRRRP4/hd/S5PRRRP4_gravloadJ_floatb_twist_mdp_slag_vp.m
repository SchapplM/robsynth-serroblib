% Calculate Gravitation load on the joints for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:21
% EndTime: 2019-12-05 16:46:21
% DurationCPUTime: 0.15s
% Computational Cost: add. (194->45), mult. (266->64), div. (0->0), fcn. (245->8), ass. (0->32)
t100 = MDP(13) + MDP(15);
t99 = MDP(14) - MDP(17);
t75 = sin(qJ(4));
t77 = cos(qJ(4));
t98 = pkin(4) * t77 + qJ(5) * t75 + pkin(3);
t72 = qJ(2) + qJ(3);
t70 = sin(t72);
t73 = sin(pkin(8));
t74 = cos(pkin(8));
t84 = g(1) * t74 + g(2) * t73;
t97 = t84 * t70;
t76 = sin(qJ(2));
t96 = pkin(2) * t76;
t71 = cos(t72);
t94 = pkin(7) * t71;
t91 = g(3) * t70;
t90 = t73 * t75;
t89 = t73 * t77;
t88 = t74 * t75;
t87 = t74 * t77;
t85 = t70 * pkin(7) + t98 * t71;
t82 = (-MDP(16) + MDP(7)) * (t84 * t71 + t91) + (t100 * t77 - t99 * t75 + MDP(6)) * (-g(3) * t71 + t97);
t58 = t71 * t90 + t87;
t60 = t71 * t88 - t89;
t47 = g(1) * t60 + g(2) * t58 + t75 * t91;
t79 = t98 * t97;
t78 = cos(qJ(2));
t64 = t74 * t94;
t62 = t73 * t94;
t61 = t71 * t87 + t90;
t59 = t71 * t89 - t88;
t1 = [(-MDP(1) - MDP(18)) * g(3); (-g(3) * t78 + t84 * t76) * MDP(3) + (g(3) * t76 + t84 * t78) * MDP(4) + (-g(1) * (-t74 * t96 + t64) - g(2) * (-t73 * t96 + t62) - g(3) * (t78 * pkin(2) + t85) + t79) * MDP(18) + t82; (-g(1) * t64 - g(2) * t62 - g(3) * t85 + t79) * MDP(18) + t82; (-g(1) * (-t60 * pkin(4) + t61 * qJ(5)) - g(2) * (-t58 * pkin(4) + t59 * qJ(5)) - (-pkin(4) * t75 + qJ(5) * t77) * t91) * MDP(18) + t99 * (g(1) * t61 + g(2) * t59 + t77 * t91) + t100 * t47; -t47 * MDP(18);];
taug = t1;
