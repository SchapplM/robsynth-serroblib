% Calculate Gravitation load on the joints for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPRP4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:53
% EndTime: 2019-12-31 16:43:54
% DurationCPUTime: 0.09s
% Computational Cost: add. (77->27), mult. (97->38), div. (0->0), fcn. (73->6), ass. (0->16)
t27 = qJ(1) + pkin(6);
t25 = sin(t27);
t26 = cos(t27);
t38 = g(1) * t26 + g(2) * t25;
t42 = MDP(10) + MDP(12);
t41 = MDP(11) - MDP(14);
t29 = sin(qJ(1));
t31 = cos(qJ(1));
t36 = g(1) * t29 - g(2) * t31;
t28 = sin(qJ(3));
t30 = cos(qJ(3));
t35 = t30 * pkin(3) + t28 * qJ(4);
t33 = pkin(2) + t35;
t32 = t36 * pkin(1);
t21 = -g(3) * t30 + t38 * t28;
t1 = [t36 * MDP(2) + (g(1) * t31 + g(2) * t29) * MDP(3) + MDP(4) * t32 - t38 * MDP(13) + (t32 + (-g(1) * pkin(5) - g(2) * t33) * t26 + (-g(2) * pkin(5) + g(1) * t33) * t25) * MDP(15) + (-t41 * t28 + t42 * t30) * (g(1) * t25 - g(2) * t26); (-MDP(15) - MDP(4)) * g(3); (-g(3) * t35 + t38 * (pkin(3) * t28 - qJ(4) * t30)) * MDP(15) + t41 * (g(3) * t28 + t38 * t30) + t42 * t21; -t21 * MDP(15);];
taug = t1;
