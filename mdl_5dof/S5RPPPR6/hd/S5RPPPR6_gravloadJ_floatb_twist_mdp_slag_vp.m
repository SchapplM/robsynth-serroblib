% Calculate Gravitation load on the joints for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:51
% EndTime: 2019-12-31 17:47:52
% DurationCPUTime: 0.30s
% Computational Cost: add. (88->55), mult. (193->82), div. (0->0), fcn. (185->8), ass. (0->32)
t61 = sin(qJ(1));
t79 = g(1) * t61;
t59 = cos(pkin(7));
t60 = sin(qJ(5));
t78 = t59 * t60;
t62 = cos(qJ(5));
t77 = t59 * t62;
t63 = cos(qJ(1));
t76 = t59 * t63;
t56 = sin(pkin(8));
t75 = t61 * t56;
t58 = cos(pkin(8));
t74 = t61 * t58;
t73 = t63 * t56;
t72 = t63 * t58;
t71 = t63 * pkin(1) + t61 * qJ(2);
t57 = sin(pkin(7));
t70 = qJ(3) * t57;
t69 = qJ(4) * t59;
t68 = MDP(11) + MDP(15);
t67 = pkin(2) * t76 + t63 * t70 + t71;
t49 = g(1) * t63 + g(2) * t61;
t48 = -g(2) * t63 + t79;
t66 = -pkin(2) * t59 - pkin(1) - t70;
t44 = -t57 * t75 + t72;
t65 = t44 * t60 + t61 * t77;
t64 = -t44 * t62 + t61 * t78;
t53 = t63 * qJ(2);
t43 = t57 * t73 + t74;
t41 = t43 * t62 + t60 * t76;
t40 = -t43 * t60 + t62 * t76;
t1 = [(-g(1) * (-t61 * pkin(1) + t53) - g(2) * t71) * MDP(7) + (-g(1) * t53 - g(2) * t67 - t66 * t79) * MDP(11) + (-g(1) * t44 - g(2) * t43) * MDP(12) + (-g(1) * (-t57 * t74 - t73) - g(2) * (t57 * t72 - t75)) * MDP(13) + (-g(1) * (t63 * pkin(3) + t53) - g(2) * (t63 * t69 + t67) + (-g(1) * (t66 - t69) - g(2) * pkin(3)) * t61) * MDP(15) + (g(1) * t64 - g(2) * t41) * MDP(21) + (g(1) * t65 - g(2) * t40) * MDP(22) + (MDP(3) - MDP(6) - MDP(8)) * t49 + (MDP(2) + (-MDP(5) + MDP(10)) * t57 + (MDP(4) - MDP(9) + MDP(14)) * t59) * t48; (-MDP(7) - t68) * t48; t68 * (g(3) * t59 - t49 * t57); (-g(3) * t57 - t49 * t59) * MDP(15); (-g(1) * t40 - g(2) * t65 - g(3) * (t56 * t78 + t57 * t62)) * MDP(21) + (g(1) * t41 + g(2) * t64 - g(3) * (t56 * t77 - t57 * t60)) * MDP(22);];
taug = t1;
