% Calculate Gravitation load on the joints for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:13
% EndTime: 2019-12-31 17:28:14
% DurationCPUTime: 0.13s
% Computational Cost: add. (102->40), mult. (166->63), div. (0->0), fcn. (166->8), ass. (0->28)
t45 = sin(qJ(2));
t64 = t45 * MDP(10);
t43 = qJ(3) + qJ(4);
t41 = sin(t43);
t42 = cos(t43);
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t63 = t47 * MDP(16) - t44 * MDP(17) + t42 * MDP(23) - t41 * MDP(24) + MDP(9);
t62 = g(3) * t45;
t46 = sin(qJ(1));
t48 = cos(qJ(2));
t61 = t46 * t48;
t49 = cos(qJ(1));
t60 = t49 * t41;
t59 = t49 * t42;
t58 = t49 * t44;
t57 = t49 * t47;
t33 = t41 * t61 + t59;
t34 = -t42 * t61 + t60;
t35 = t46 * t42 - t48 * t60;
t36 = t46 * t41 + t48 * t59;
t56 = (-g(1) * t35 + g(2) * t33 + t41 * t62) * MDP(23) + (g(1) * t36 - g(2) * t34 + t42 * t62) * MDP(24);
t51 = g(1) * t49 + g(2) * t46;
t40 = t46 * t44 + t48 * t57;
t39 = t46 * t47 - t48 * t58;
t38 = -t47 * t61 + t58;
t37 = t44 * t61 + t57;
t1 = [t51 * MDP(3) + (-g(1) * t38 - g(2) * t40) * MDP(16) + (-g(1) * t37 - g(2) * t39) * MDP(17) + (-g(1) * t34 - g(2) * t36) * MDP(23) + (-g(1) * t33 - g(2) * t35) * MDP(24) + (t48 * MDP(9) + MDP(2) - t64) * (g(1) * t46 - g(2) * t49); (-t63 * t48 + t64) * g(3) + (MDP(10) * t48 + t63 * t45) * t51; (-g(1) * t39 + g(2) * t37 + t44 * t62) * MDP(16) + (g(1) * t40 - g(2) * t38 + t47 * t62) * MDP(17) + t56; t56;];
taug = t1;
