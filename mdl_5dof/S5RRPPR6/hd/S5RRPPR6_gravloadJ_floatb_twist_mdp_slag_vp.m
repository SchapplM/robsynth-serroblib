% Calculate Gravitation load on the joints for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:47:33
% EndTime: 2021-01-15 19:47:35
% DurationCPUTime: 0.26s
% Computational Cost: add. (181->67), mult. (232->94), div. (0->0), fcn. (206->12), ass. (0->40)
t95 = MDP(12) - MDP(17);
t69 = qJ(2) + pkin(8);
t63 = sin(t69);
t93 = g(3) * t63;
t65 = cos(t69);
t92 = g(3) * t65;
t68 = pkin(9) + qJ(5);
t62 = sin(t68);
t76 = sin(qJ(1));
t91 = t76 * t62;
t64 = cos(t68);
t90 = t76 * t64;
t70 = sin(pkin(9));
t89 = t76 * t70;
t72 = cos(pkin(9));
t88 = t76 * t72;
t74 = -qJ(3) - pkin(6);
t87 = t76 * t74;
t78 = cos(qJ(1));
t86 = t78 * t62;
t85 = t78 * t64;
t84 = t78 * t70;
t83 = t78 * t72;
t82 = g(1) * t78 + g(2) * t76;
t58 = g(1) * t76 - g(2) * t78;
t77 = cos(qJ(2));
t75 = sin(qJ(2));
t73 = cos(pkin(8));
t71 = sin(pkin(8));
t66 = t77 * pkin(2);
t61 = t66 + pkin(1);
t60 = t74 * t78;
t57 = -t71 * pkin(3) + qJ(4) * t73;
t56 = pkin(3) * t73 + qJ(4) * t71 + pkin(2);
t53 = t65 * t85 + t91;
t52 = -t65 * t86 + t90;
t51 = -t65 * t90 + t86;
t50 = t65 * t91 + t85;
t47 = t56 * t77 + t57 * t75 + pkin(1);
t1 = [(-g(1) * (-t76 * t61 - t60) - g(2) * (t78 * t61 - t87)) * MDP(14) + (-g(1) * (-t65 * t88 + t84) - g(2) * (t65 * t83 + t89)) * MDP(15) + (-g(1) * (t65 * t89 + t83) - g(2) * (-t65 * t84 + t88)) * MDP(16) + (-g(1) * (-t47 * t76 - t60) - g(2) * (t47 * t78 - t87)) * MDP(18) + (-g(1) * t51 - g(2) * t53) * MDP(24) + (-g(1) * t50 - g(2) * t52) * MDP(25) + (MDP(3) - MDP(13)) * t82 + (-t75 * MDP(10) + t65 * MDP(11) + t77 * MDP(9) - t95 * t63 + MDP(2)) * t58; (g(3) * t75 + t82 * t77) * MDP(10) + (-g(3) * (t65 * pkin(3) + t63 * qJ(4) + t66) - t82 * (-t56 * t75 + t57 * t77)) * MDP(18) + (pkin(2) * MDP(14) + MDP(9)) * (-g(3) * t77 + t82 * t75) + t95 * (t82 * t65 + t93) + (t72 * MDP(15) - t70 * MDP(16) + t64 * MDP(24) - t62 * MDP(25) + MDP(11)) * (t82 * t63 - t92); (-MDP(14) - MDP(18)) * t58; (t92 - t82 * (t71 * t77 + t73 * t75)) * MDP(18); (-g(1) * t52 + g(2) * t50 + t62 * t93) * MDP(24) + (g(1) * t53 - g(2) * t51 + t64 * t93) * MDP(25);];
taug = t1;
