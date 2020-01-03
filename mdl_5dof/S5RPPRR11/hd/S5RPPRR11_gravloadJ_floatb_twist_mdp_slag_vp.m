% Calculate Gravitation load on the joints for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR11_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR11_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPPRR11_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:54
% EndTime: 2019-12-31 18:05:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (57->38), mult. (123->52), div. (0->0), fcn. (106->6), ass. (0->20)
t40 = cos(qJ(4));
t48 = g(3) * t40;
t36 = sin(qJ(5));
t38 = sin(qJ(1));
t47 = t38 * t36;
t39 = cos(qJ(5));
t46 = t38 * t39;
t41 = cos(qJ(1));
t45 = t41 * t36;
t44 = t41 * t39;
t43 = t41 * pkin(1) + t38 * qJ(2);
t31 = g(1) * t41 + g(2) * t38;
t30 = g(1) * t38 - g(2) * t41;
t37 = sin(qJ(4));
t33 = t41 * qJ(2);
t29 = t37 * t44 - t47;
t28 = -t37 * t45 - t46;
t27 = -t37 * t46 - t45;
t26 = t37 * t47 - t44;
t1 = [(-g(1) * (-t38 * pkin(1) + t33) - g(2) * t43) * MDP(6) + (-g(1) * (t33 + (-pkin(1) - qJ(3)) * t38) - g(2) * (t41 * qJ(3) + t43)) * MDP(9) + (-g(1) * t27 - g(2) * t29) * MDP(22) + (-g(1) * t26 - g(2) * t28) * MDP(23) + (MDP(3) - MDP(5) - MDP(7)) * t31 + (t37 * MDP(15) + t40 * MDP(16) + MDP(2) - MDP(4) + MDP(8)) * t30; (-MDP(6) - MDP(9)) * t30; -t31 * MDP(9); (t31 * t37 + t48) * MDP(16) + (-MDP(22) * t39 + MDP(23) * t36 - MDP(15)) * (-g(3) * t37 + t31 * t40); (-g(1) * t28 + g(2) * t26 + t36 * t48) * MDP(22) + (g(1) * t29 - g(2) * t27 + t39 * t48) * MDP(23);];
taug = t1;
