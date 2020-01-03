% Calculate Gravitation load on the joints for
% S5RPRRP9
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
%   see S5RPRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP9_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRP9_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:11
% EndTime: 2019-12-31 18:49:12
% DurationCPUTime: 0.21s
% Computational Cost: add. (208->48), mult. (183->61), div. (0->0), fcn. (144->8), ass. (0->25)
t54 = pkin(8) + qJ(3);
t50 = cos(t54);
t51 = qJ(4) + t54;
t47 = sin(t51);
t48 = cos(t51);
t63 = t48 * pkin(4) + t47 * qJ(5);
t67 = pkin(3) * t50 + t63;
t66 = MDP(20) + MDP(22);
t65 = MDP(21) - MDP(24);
t64 = pkin(4) * t47;
t62 = qJ(5) * t48;
t57 = sin(qJ(1));
t58 = cos(qJ(1));
t43 = g(1) * t58 + g(2) * t57;
t35 = -g(3) * t48 + t43 * t47;
t61 = t65 * (g(3) * t47 + t43 * t48) + t66 * t35;
t49 = sin(t54);
t60 = -pkin(3) * t49 - t64;
t42 = g(1) * t57 - g(2) * t58;
t56 = cos(pkin(8));
t59 = t56 * pkin(2) + pkin(1) + t67;
t53 = -pkin(7) - pkin(6) - qJ(2);
t41 = t58 * t62;
t40 = t57 * t62;
t1 = [(-g(1) * (-t57 * pkin(1) + t58 * qJ(2)) - g(2) * (t58 * pkin(1) + t57 * qJ(2))) * MDP(7) + ((g(1) * t53 - g(2) * t59) * t58 + (g(1) * t59 + g(2) * t53) * t57) * MDP(25) + (MDP(3) - MDP(6) - MDP(23)) * t43 + (t66 * t48 - t65 * t47 + t50 * MDP(13) - t49 * MDP(14) + MDP(4) * t56 - MDP(5) * sin(pkin(8)) + MDP(2)) * t42; (-MDP(25) - MDP(7)) * t42; (-g(3) * t50 + t43 * t49) * MDP(13) + (g(3) * t49 + t43 * t50) * MDP(14) + (-g(1) * (t60 * t58 + t41) - g(2) * (t60 * t57 + t40) - g(3) * t67) * MDP(25) + t61; (-g(1) * (-t58 * t64 + t41) - g(2) * (-t57 * t64 + t40) - g(3) * t63) * MDP(25) + t61; -t35 * MDP(25);];
taug = t1;
