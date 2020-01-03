% Calculate Gravitation load on the joints for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR6_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRPR6_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:42
% EndTime: 2019-12-31 16:24:42
% DurationCPUTime: 0.10s
% Computational Cost: add. (54->25), mult. (98->39), div. (0->0), fcn. (86->8), ass. (0->13)
t27 = sin(pkin(6));
t29 = cos(pkin(6));
t34 = g(1) * t29 + g(2) * t27;
t30 = sin(qJ(2));
t31 = cos(qJ(2));
t21 = -g(3) * t31 + t34 * t30;
t38 = g(3) * t30;
t36 = t27 * t31;
t35 = t29 * t31;
t25 = pkin(7) + qJ(4);
t24 = cos(t25);
t23 = sin(t25);
t1 = [(-MDP(1) - MDP(8)) * g(3); (-g(3) * (t31 * pkin(2) + t30 * qJ(3)) + t34 * (pkin(2) * t30 - qJ(3) * t31)) * MDP(8) + (MDP(4) - MDP(7)) * (t34 * t31 + t38) + (MDP(14) * t24 - MDP(15) * t23 + MDP(5) * cos(pkin(7)) - MDP(6) * sin(pkin(7)) + MDP(3)) * t21; -t21 * MDP(8); (-g(1) * (-t23 * t35 + t27 * t24) - g(2) * (-t23 * t36 - t29 * t24) + t23 * t38) * MDP(14) + (-g(1) * (-t27 * t23 - t24 * t35) - g(2) * (t29 * t23 - t24 * t36) + t24 * t38) * MDP(15);];
taug = t1;
