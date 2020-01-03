% Calculate Gravitation load on the joints for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RPPRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:51:02
% EndTime: 2019-12-31 17:51:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (81->32), mult. (87->42), div. (0->0), fcn. (59->6), ass. (0->15)
t32 = sin(qJ(4));
t41 = pkin(4) * t32;
t40 = -MDP(16) - MDP(7);
t30 = qJ(1) + pkin(7);
t27 = sin(t30);
t28 = cos(t30);
t35 = cos(qJ(1));
t39 = t35 * pkin(1) + t28 * pkin(2) + t27 * qJ(3);
t33 = sin(qJ(1));
t38 = -t33 * pkin(1) + t28 * qJ(3);
t22 = -g(1) * t28 - g(2) * t27;
t21 = g(1) * t27 - g(2) * t28;
t34 = cos(qJ(4));
t31 = -qJ(5) - pkin(6);
t1 = [(g(1) * t35 + g(2) * t33) * MDP(3) + (-g(1) * (-t27 * pkin(2) + t38) - g(2) * t39) * MDP(7) + (-g(1) * (t28 * t41 + (-pkin(2) + t31) * t27 + t38) - g(2) * (t27 * t41 - t28 * t31 + t39)) * MDP(16) + (-MDP(5) + MDP(15)) * t21 + (pkin(1) * MDP(4) + MDP(2)) * (g(1) * t33 - g(2) * t35) + (t32 * MDP(13) + t34 * MDP(14) + MDP(6)) * t22; (-MDP(4) + t40) * g(3); t40 * t21; (g(3) * t34 + t21 * t32) * MDP(14) + (MDP(16) * pkin(4) + MDP(13)) * (g(3) * t32 - t21 * t34); t22 * MDP(16);];
taug = t1;
