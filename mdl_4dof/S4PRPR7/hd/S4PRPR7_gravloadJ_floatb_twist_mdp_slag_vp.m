% Calculate Gravitation load on the joints for
% S4PRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRPR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:55
% EndTime: 2019-12-31 16:25:55
% DurationCPUTime: 0.08s
% Computational Cost: add. (35->23), mult. (86->37), div. (0->0), fcn. (75->6), ass. (0->12)
t22 = sin(pkin(6));
t23 = cos(pkin(6));
t30 = g(1) * t23 + g(2) * t22;
t27 = cos(qJ(2));
t33 = g(3) * t27;
t24 = sin(qJ(4));
t25 = sin(qJ(2));
t32 = t24 * t25;
t26 = cos(qJ(4));
t31 = t25 * t26;
t19 = t30 * t25 - t33;
t1 = [(-MDP(1) - MDP(7)) * g(3); (-g(3) * (t27 * pkin(2) + t25 * qJ(3)) + t30 * (pkin(2) * t25 - qJ(3) * t27)) * MDP(7) + (MDP(3) - MDP(5)) * t19 + (-MDP(13) * t24 - MDP(14) * t26 + MDP(4) - MDP(6)) * (g(3) * t25 + t30 * t27); -t19 * MDP(7); (-g(1) * (-t22 * t24 + t23 * t31) - g(2) * (t22 * t31 + t23 * t24) + t26 * t33) * MDP(13) + (-g(1) * (-t22 * t26 - t23 * t32) - g(2) * (-t22 * t32 + t23 * t26) - t24 * t33) * MDP(14);];
taug = t1;
