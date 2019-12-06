% Calculate Gravitation load on the joints for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(1,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [1x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RPRRR1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:09:52
% EndTime: 2019-12-05 18:09:53
% DurationCPUTime: 0.21s
% Computational Cost: add. (91->50), mult. (230->81), div. (0->0), fcn. (240->8), ass. (0->31)
t50 = sin(qJ(1));
t54 = cos(qJ(1));
t75 = -g(1) * t54 - g(2) * t50;
t49 = sin(qJ(3));
t72 = g(3) * t49;
t47 = sin(qJ(5));
t70 = t49 * t47;
t51 = cos(qJ(5));
t69 = t49 * t51;
t52 = cos(qJ(4));
t68 = t49 * t52;
t67 = t49 * t54;
t53 = cos(qJ(3));
t66 = t50 * t53;
t65 = t53 * t47;
t64 = t53 * t51;
t48 = sin(qJ(4));
t63 = t54 * t48;
t62 = t54 * t52;
t44 = g(1) * t50 - g(2) * t54;
t41 = t52 * t66 - t63;
t60 = -t41 * t47 + t50 * t69;
t59 = t41 * t51 + t50 * t70;
t58 = t51 * t68 - t65;
t57 = t47 * t68 + t64;
t43 = t50 * t48 + t53 * t62;
t42 = t50 * t52 - t53 * t63;
t40 = t48 * t66 + t62;
t39 = t43 * t51 + t47 * t67;
t38 = -t43 * t47 + t51 * t67;
t1 = [(g(1) * t41 - g(2) * t43) * MDP(19) + (-g(1) * t40 - g(2) * t42) * MDP(20) + (g(1) * t59 - g(2) * t39) * MDP(26) + (g(1) * t60 - g(2) * t38) * MDP(27) - (-MDP(6) * qJ(2) + MDP(3) - MDP(5)) * t75 + (t53 * MDP(12) - t49 * MDP(13) + MDP(2) + MDP(4)) * t44; -t44 * MDP(6); (-t53 * t75 + t72) * MDP(13) + (-g(3) * (t52 * t64 + t70) - t75 * t58) * MDP(26) + (-g(3) * (-t52 * t65 + t69) + t75 * t57) * MDP(27) + (t52 * MDP(19) - MDP(20) * t48 + MDP(12)) * (-g(3) * t53 - t49 * t75); (g(1) * t43 + g(2) * t41 + g(3) * t68) * MDP(20) + (MDP(26) * t51 - MDP(27) * t47 + MDP(19)) * (-g(1) * t42 + g(2) * t40 + t48 * t72); (-g(1) * t38 - g(2) * t60 + g(3) * t57) * MDP(26) + (g(1) * t39 + g(2) * t59 + g(3) * t58) * MDP(27);];
taug = t1;
