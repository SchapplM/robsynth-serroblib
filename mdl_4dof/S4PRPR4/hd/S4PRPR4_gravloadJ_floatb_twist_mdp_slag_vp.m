% Calculate Gravitation load on the joints for
% S4PRPR4
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
%   see S4PRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR4_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4PRPR4_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:00
% EndTime: 2019-12-31 16:22:01
% DurationCPUTime: 0.04s
% Computational Cost: add. (44->17), mult. (46->23), div. (0->0), fcn. (32->4), ass. (0->7)
t14 = pkin(6) + qJ(2);
t12 = sin(t14);
t13 = cos(t14);
t9 = g(1) * t12 - g(2) * t13;
t16 = cos(qJ(4));
t15 = sin(qJ(4));
t1 = [(-MDP(1) - MDP(7)) * g(3); (-g(1) * (-t12 * pkin(2) + t13 * qJ(3)) - g(2) * (t13 * pkin(2) + t12 * qJ(3))) * MDP(7) + (MDP(3) - MDP(5)) * t9 + (-t15 * MDP(13) - t16 * MDP(14) + MDP(4) - MDP(6)) * (g(1) * t13 + g(2) * t12); -t9 * MDP(7); (g(3) * t15 - t9 * t16) * MDP(13) + (g(3) * t16 + t9 * t15) * MDP(14);];
taug = t1;
