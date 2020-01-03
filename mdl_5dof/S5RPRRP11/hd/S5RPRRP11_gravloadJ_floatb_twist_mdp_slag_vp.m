% Calculate Gravitation load on the joints for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRP11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:33
% EndTime: 2019-12-31 18:54:34
% DurationCPUTime: 0.29s
% Computational Cost: add. (199->60), mult. (279->84), div. (0->0), fcn. (265->8), ass. (0->27)
t62 = sin(qJ(1));
t64 = cos(qJ(1));
t50 = g(1) * t64 + g(2) * t62;
t79 = MDP(14) - MDP(23);
t78 = MDP(20) + MDP(22);
t77 = MDP(21) - MDP(24);
t57 = pkin(8) + qJ(3);
t54 = sin(t57);
t74 = g(3) * t54;
t61 = sin(qJ(4));
t73 = t62 * t61;
t63 = cos(qJ(4));
t72 = t62 * t63;
t71 = t63 * t64;
t70 = t64 * t61;
t49 = g(1) * t62 - g(2) * t64;
t55 = cos(t57);
t59 = cos(pkin(8));
t68 = pkin(2) * t59 + pkin(3) * t55 + pkin(7) * t54 + pkin(1);
t67 = pkin(4) * t63 + qJ(5) * t61 + pkin(3);
t45 = t55 * t73 + t71;
t47 = t55 * t70 - t72;
t39 = g(1) * t47 + g(2) * t45 + t61 * t74;
t60 = -pkin(6) - qJ(2);
t48 = t55 * t71 + t73;
t46 = t55 * t72 - t70;
t1 = [(-g(1) * (-t62 * pkin(1) + qJ(2) * t64) - g(2) * (pkin(1) * t64 + t62 * qJ(2))) * MDP(7) + (-g(1) * (-t46 * pkin(4) - t45 * qJ(5)) - g(2) * (t48 * pkin(4) + t47 * qJ(5)) + (g(1) * t60 - g(2) * t68) * t64 + (g(1) * t68 + g(2) * t60) * t62) * MDP(25) - t77 * (g(1) * t45 - g(2) * t47) + (MDP(3) - MDP(6)) * t50 + t78 * (g(1) * t46 - g(2) * t48) + (t55 * MDP(13) + MDP(4) * t59 - MDP(5) * sin(pkin(8)) - t79 * t54 + MDP(2)) * t49; (-MDP(25) - MDP(7)) * t49; ((-t50 * pkin(7) - g(3) * t67) * t55 + (-g(3) * pkin(7) + t50 * t67) * t54) * MDP(25) + t79 * (t50 * t55 + t74) + (-t77 * t61 + t78 * t63 + MDP(13)) * (-g(3) * t55 + t50 * t54); (-g(1) * (-pkin(4) * t47 + qJ(5) * t48) - g(2) * (-pkin(4) * t45 + qJ(5) * t46) - (-pkin(4) * t61 + qJ(5) * t63) * t74) * MDP(25) + t77 * (g(1) * t48 + g(2) * t46 + t63 * t74) + t78 * t39; -t39 * MDP(25);];
taug = t1;
