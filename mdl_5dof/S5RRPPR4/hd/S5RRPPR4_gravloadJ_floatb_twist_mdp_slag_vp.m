% Calculate Gravitation load on the joints for
% S5RRPPR4
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
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5RRPPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:49
% EndTime: 2019-12-31 19:27:49
% DurationCPUTime: 0.08s
% Computational Cost: add. (166->39), mult. (157->52), div. (0->0), fcn. (142->8), ass. (0->22)
t56 = sin(qJ(1));
t67 = t56 * pkin(1);
t53 = qJ(1) + qJ(2);
t50 = sin(t53);
t51 = cos(t53);
t66 = t51 * pkin(2) + t50 * qJ(3);
t65 = cos(pkin(8));
t58 = cos(qJ(1));
t64 = t58 * pkin(1) + t66;
t46 = t51 * qJ(3);
t63 = -t50 * pkin(2) + t46;
t62 = t46 + (-pkin(2) - pkin(3)) * t50;
t54 = sin(pkin(8));
t36 = -t50 * t65 + t51 * t54;
t37 = t50 * t54 + t51 * t65;
t61 = g(1) * t37 - g(2) * t36;
t42 = g(1) * t50 - g(2) * t51;
t55 = sin(qJ(5));
t57 = cos(qJ(5));
t59 = -t61 * MDP(11) + (MDP(6) - MDP(8)) * (g(1) * t51 + g(2) * t50) + (MDP(5) + MDP(7)) * t42 + (-t57 * MDP(18) + t55 * MDP(19) - MDP(10)) * (g(1) * t36 + g(2) * t37);
t47 = t51 * pkin(3);
t1 = [(g(1) * t56 - g(2) * t58) * MDP(2) + (g(1) * t58 + g(2) * t56) * MDP(3) + (-g(1) * (t63 - t67) - g(2) * t64) * MDP(9) + (-g(1) * (t62 - t67) - g(2) * (t47 + t64)) * MDP(12) + t59; (-g(1) * t63 - g(2) * t66) * MDP(9) + (-g(1) * t62 - g(2) * (t47 + t66)) * MDP(12) + t59; (-MDP(12) - MDP(9)) * t42; g(3) * MDP(12); (g(3) * t57 + t61 * t55) * MDP(18) + (-g(3) * t55 + t61 * t57) * MDP(19);];
taug = t1;
