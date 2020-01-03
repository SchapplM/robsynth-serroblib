% Calculate Gravitation load on the joints for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:44
% EndTime: 2019-12-31 16:40:44
% DurationCPUTime: 0.10s
% Computational Cost: add. (51->32), mult. (116->48), div. (0->0), fcn. (105->6), ass. (0->19)
t39 = sin(qJ(1));
t46 = g(1) * t39;
t41 = cos(qJ(1));
t45 = t41 * pkin(1) + t39 * qJ(2);
t31 = g(1) * t41 + g(2) * t39;
t30 = -g(2) * t41 + t46;
t36 = sin(pkin(6));
t37 = cos(pkin(6));
t44 = pkin(2) * t37 + qJ(3) * t36;
t38 = sin(qJ(4));
t40 = cos(qJ(4));
t43 = t36 * t40 - t37 * t38;
t42 = t36 * t38 + t37 * t40;
t33 = t41 * qJ(2);
t27 = t42 * t41;
t26 = t43 * t41;
t25 = t42 * t39;
t24 = t43 * t39;
t1 = [(-g(1) * (-t39 * pkin(1) + t33) - g(2) * t45) * MDP(7) + (-g(1) * t33 - g(2) * (t44 * t41 + t45) - (-pkin(1) - t44) * t46) * MDP(11) + (g(1) * t25 - g(2) * t27) * MDP(17) + (g(1) * t24 - g(2) * t26) * MDP(18) + (MDP(3) - MDP(6) - MDP(9)) * t31 + (MDP(2) + (MDP(4) + MDP(8)) * t37 + (-MDP(5) + MDP(10)) * t36) * t30; (-MDP(11) - MDP(7)) * t30; (g(3) * t37 - t31 * t36) * MDP(11); (-g(1) * t26 - g(2) * t24 + g(3) * t42) * MDP(17) + (g(1) * t27 + g(2) * t25 + g(3) * t43) * MDP(18);];
taug = t1;
