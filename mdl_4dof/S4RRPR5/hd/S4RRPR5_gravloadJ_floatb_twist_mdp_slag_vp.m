% Calculate Gravitation load on the joints for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RRPR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:33
% EndTime: 2019-12-31 17:03:33
% DurationCPUTime: 0.09s
% Computational Cost: add. (80->24), mult. (81->33), div. (0->0), fcn. (58->6), ass. (0->12)
t32 = qJ(1) + qJ(2);
t30 = sin(t32);
t31 = cos(t32);
t39 = t31 * pkin(2) + t30 * qJ(3);
t38 = -t30 * pkin(2) + t31 * qJ(3);
t24 = g(1) * t30 - g(2) * t31;
t33 = sin(qJ(4));
t35 = cos(qJ(4));
t37 = (MDP(5) - MDP(7)) * t24 + (-t33 * MDP(15) - t35 * MDP(16) + MDP(6) - MDP(8)) * (g(1) * t31 + g(2) * t30);
t36 = cos(qJ(1));
t34 = sin(qJ(1));
t1 = [(g(1) * t34 - g(2) * t36) * MDP(2) + (g(1) * t36 + g(2) * t34) * MDP(3) + (-g(1) * (-t34 * pkin(1) + t38) - g(2) * (t36 * pkin(1) + t39)) * MDP(9) + t37; (-g(1) * t38 - g(2) * t39) * MDP(9) + t37; -t24 * MDP(9); (g(3) * t33 - t24 * t35) * MDP(15) + (g(3) * t35 + t24 * t33) * MDP(16);];
taug = t1;
