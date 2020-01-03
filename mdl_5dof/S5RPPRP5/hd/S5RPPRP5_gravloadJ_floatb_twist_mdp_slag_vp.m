% Calculate Gravitation load on the joints for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPPRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:43
% EndTime: 2019-12-31 17:53:44
% DurationCPUTime: 0.26s
% Computational Cost: add. (115->56), mult. (263->75), div. (0->0), fcn. (260->6), ass. (0->27)
t62 = cos(pkin(7));
t63 = sin(qJ(4));
t61 = sin(pkin(7));
t77 = cos(qJ(4));
t71 = t61 * t77;
t49 = -t62 * t63 + t71;
t80 = MDP(17) + MDP(19);
t79 = MDP(18) - MDP(21);
t64 = sin(qJ(1));
t78 = g(1) * t64;
t65 = cos(qJ(1));
t75 = t62 * t65;
t74 = t65 * pkin(1) + t64 * qJ(2);
t73 = qJ(3) * t61;
t72 = MDP(11) + MDP(22);
t70 = pkin(2) * t75 + t65 * t73 + t74;
t51 = g(1) * t65 + g(2) * t64;
t50 = -g(2) * t65 + t78;
t68 = -pkin(2) * t62 - pkin(1) - t73;
t48 = t61 * t63 + t62 * t77;
t42 = t49 * t64;
t44 = t63 * t75 - t65 * t71;
t67 = g(1) * t44 - g(2) * t42 + g(3) * t48;
t58 = t65 * qJ(2);
t45 = t48 * t65;
t43 = t48 * t64;
t1 = [(-g(1) * (-t64 * pkin(1) + t58) - g(2) * t74) * MDP(7) + (-g(1) * t58 - g(2) * t70 - t68 * t78) * MDP(11) + (-g(1) * (-t43 * pkin(4) - pkin(6) * t65 + t42 * qJ(5) + t58) - g(2) * (pkin(3) * t75 + t45 * pkin(4) + t44 * qJ(5) + t70) + (-g(1) * (-pkin(3) * t62 + t68) + g(2) * pkin(6)) * t64) * MDP(22) - t79 * (-g(1) * t42 - g(2) * t44) + t80 * (g(1) * t43 - g(2) * t45) + (MDP(3) - MDP(6) - MDP(9) + MDP(20)) * t51 + (MDP(2) + (MDP(4) + MDP(8)) * t62 + (-MDP(5) + MDP(10)) * t61) * t50; (-MDP(7) - t72) * t50; t72 * (g(3) * t62 - t51 * t61); (-g(1) * (-pkin(4) * t44 + qJ(5) * t45) - g(2) * (pkin(4) * t42 + qJ(5) * t43) - g(3) * (-pkin(4) * t48 + qJ(5) * t49)) * MDP(22) + t80 * t67 + t79 * (g(1) * t45 + g(2) * t43 + g(3) * t49); -t67 * MDP(22);];
taug = t1;
