% Calculate Gravitation load on the joints for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:57
% EndTime: 2019-12-31 19:54:57
% DurationCPUTime: 0.21s
% Computational Cost: add. (204->50), mult. (187->61), div. (0->0), fcn. (146->8), ass. (0->27)
t63 = cos(qJ(2));
t56 = t63 * pkin(2);
t59 = qJ(2) + pkin(8);
t55 = qJ(4) + t59;
t52 = sin(t55);
t53 = cos(t55);
t71 = t53 * pkin(4) + t52 * qJ(5);
t75 = pkin(3) * cos(t59) + t56 + t71;
t74 = MDP(18) + MDP(20);
t73 = MDP(19) - MDP(22);
t72 = pkin(4) * t52;
t60 = -qJ(3) - pkin(6);
t69 = qJ(5) * t53;
t61 = sin(qJ(2));
t68 = -pkin(3) * sin(t59) - t61 * pkin(2) - t72;
t62 = sin(qJ(1));
t64 = cos(qJ(1));
t48 = g(1) * t64 + g(2) * t62;
t39 = -g(3) * t53 + t48 * t52;
t67 = t73 * (g(3) * t52 + t48 * t53) + t74 * t39;
t47 = g(1) * t62 - g(2) * t64;
t66 = pkin(1) + t75;
t58 = -pkin(7) + t60;
t54 = t56 + pkin(1);
t46 = t64 * t69;
t45 = t62 * t69;
t1 = [(-g(1) * (-t62 * t54 - t64 * t60) - g(2) * (t64 * t54 - t62 * t60)) * MDP(12) + ((g(1) * t58 - g(2) * t66) * t64 + (g(1) * t66 + g(2) * t58) * t62) * MDP(23) + (MDP(3) - MDP(11) - MDP(21)) * t48 + (-t61 * MDP(10) + t63 * MDP(9) - t73 * t52 + t74 * t53 + MDP(2)) * t47; (g(3) * t61 + t48 * t63) * MDP(10) + (-g(1) * (t68 * t64 + t46) - g(2) * (t68 * t62 + t45) - g(3) * t75) * MDP(23) + t67 + (MDP(12) * pkin(2) + MDP(9)) * (-g(3) * t63 + t48 * t61); (-MDP(12) - MDP(23)) * t47; (-g(1) * (-t64 * t72 + t46) - g(2) * (-t62 * t72 + t45) - g(3) * t71) * MDP(23) + t67; -t39 * MDP(23);];
taug = t1;
