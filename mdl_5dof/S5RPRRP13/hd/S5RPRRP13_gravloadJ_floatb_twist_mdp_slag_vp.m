% Calculate Gravitation load on the joints for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP13_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP13_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP13_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:38
% EndTime: 2019-12-31 18:59:39
% DurationCPUTime: 0.37s
% Computational Cost: add. (119->59), mult. (272->79), div. (0->0), fcn. (259->6), ass. (0->26)
t63 = sin(qJ(1));
t66 = cos(qJ(1));
t89 = -g(1) * t63 + g(2) * t66;
t62 = sin(qJ(3));
t65 = cos(qJ(3));
t76 = t65 * pkin(7);
t88 = t62 * pkin(3) - t76;
t87 = MDP(13) - MDP(22);
t86 = MDP(19) + MDP(21);
t85 = MDP(20) - MDP(23);
t78 = g(3) * t65;
t61 = sin(qJ(4));
t75 = t63 * t61;
t64 = cos(qJ(4));
t74 = t63 * t64;
t73 = t66 * t61;
t72 = t66 * t64;
t71 = t66 * pkin(1) + t63 * qJ(2);
t69 = pkin(4) * t64 + qJ(5) * t61 + pkin(3);
t49 = t62 * t75 - t72;
t51 = t62 * t73 + t74;
t44 = g(1) * t49 - g(2) * t51 + t61 * t78;
t58 = t66 * qJ(2);
t52 = t62 * t72 - t75;
t50 = t62 * t74 + t73;
t1 = [(-g(1) * (-t63 * pkin(1) + t58) - g(2) * t71) * MDP(6) + (-g(1) * (t52 * pkin(4) + t51 * qJ(5) + t88 * t66 + t58) - g(2) * (t50 * pkin(4) + t66 * pkin(6) + t49 * qJ(5) + t71) + (-g(1) * (-pkin(1) - pkin(6)) - g(2) * t88) * t63) * MDP(24) + t85 * (g(1) * t51 + g(2) * t49) - (MDP(2) - MDP(4)) * t89 + t86 * (-g(1) * t52 - g(2) * t50) + (-t62 * MDP(12) - t87 * t65 + MDP(3) - MDP(5)) * (g(1) * t66 + g(2) * t63); -(-MDP(24) - MDP(6)) * t89; (-g(3) * (-t62 * t69 + t76) + t89 * (pkin(7) * t62 + t65 * t69)) * MDP(24) + t87 * (-t62 * t89 + t78) + (t85 * t61 - t86 * t64 - MDP(12)) * (-g(3) * t62 - t89 * t65); (-g(1) * (-t49 * pkin(4) + t50 * qJ(5)) - g(2) * (t51 * pkin(4) - t52 * qJ(5)) - (-pkin(4) * t61 + qJ(5) * t64) * t78) * MDP(24) + t85 * (g(1) * t50 - g(2) * t52 + t64 * t78) + t86 * t44; -t44 * MDP(24);];
taug = t1;
