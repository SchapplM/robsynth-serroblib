% Calculate Gravitation load on the joints for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRRP6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:11
% EndTime: 2019-12-31 17:19:11
% DurationCPUTime: 0.15s
% Computational Cost: add. (73->38), mult. (157->57), div. (0->0), fcn. (139->6), ass. (0->23)
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t52 = g(1) * t47 + g(2) * t44;
t65 = MDP(10) - MDP(18);
t43 = sin(qJ(2));
t46 = cos(qJ(2));
t33 = -g(3) * t46 + t43 * t52;
t59 = g(3) * t43;
t57 = t44 * t46;
t42 = sin(qJ(3));
t56 = t47 * t42;
t45 = cos(qJ(3));
t55 = t47 * t45;
t53 = pkin(3) * t42 + pkin(5);
t40 = pkin(3) * t45 + pkin(2);
t41 = -qJ(4) - pkin(6);
t50 = t40 * t46 - t41 * t43;
t48 = pkin(1) + t50;
t37 = t44 * t45 - t46 * t56;
t35 = t42 * t57 + t55;
t38 = t42 * t44 + t46 * t55;
t36 = -t45 * t57 + t56;
t1 = [t52 * MDP(3) + (-g(1) * t36 - g(2) * t38) * MDP(16) + (-g(1) * t35 - g(2) * t37) * MDP(17) + ((-g(1) * t53 - g(2) * t48) * t47 + (g(1) * t48 - g(2) * t53) * t44) * MDP(19) + (t46 * MDP(9) - t43 * t65 + MDP(2)) * (g(1) * t44 - g(2) * t47); (-g(3) * t50 + t52 * (t40 * t43 + t41 * t46)) * MDP(19) + t65 * (t46 * t52 + t59) + (MDP(16) * t45 - MDP(17) * t42 + MDP(9)) * t33; (g(1) * t38 - g(2) * t36 + t45 * t59) * MDP(17) + (MDP(19) * pkin(3) + MDP(16)) * (-g(1) * t37 + g(2) * t35 + t42 * t59); -t33 * MDP(19);];
taug = t1;
