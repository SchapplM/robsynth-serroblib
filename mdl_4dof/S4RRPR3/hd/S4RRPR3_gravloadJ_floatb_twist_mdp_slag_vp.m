% Calculate Gravitation load on the joints for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4RRPR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:34
% EndTime: 2019-12-31 17:01:34
% DurationCPUTime: 0.04s
% Computational Cost: add. (66->22), mult. (63->35), div. (0->0), fcn. (44->8), ass. (0->14)
t26 = qJ(1) + qJ(2);
t23 = pkin(7) + t26;
t21 = sin(t23);
t22 = cos(t23);
t24 = sin(t26);
t25 = cos(t26);
t27 = sin(qJ(4));
t29 = cos(qJ(4));
t31 = g(1) * t24 - g(2) * t25;
t34 = t31 * MDP(5) + (g(1) * t25 + g(2) * t24) * MDP(6) + (t29 * MDP(13) - t27 * MDP(14)) * (g(1) * t21 - g(2) * t22);
t33 = g(1) * t22 + g(2) * t21;
t30 = cos(qJ(1));
t28 = sin(qJ(1));
t1 = [(g(1) * t28 - g(2) * t30) * MDP(2) + (g(1) * t30 + g(2) * t28) * MDP(3) + (-g(1) * (-t28 * pkin(1) - pkin(2) * t24) - g(2) * (t30 * pkin(1) + pkin(2) * t25)) * MDP(7) + t34; t31 * MDP(7) * pkin(2) + t34; -g(3) * MDP(7); (-g(3) * t29 + t33 * t27) * MDP(13) + (g(3) * t27 + t33 * t29) * MDP(14);];
taug = t1;
