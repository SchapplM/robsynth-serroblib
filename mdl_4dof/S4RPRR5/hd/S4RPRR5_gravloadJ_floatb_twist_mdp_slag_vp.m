% Calculate Gravitation load on the joints for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RPRR5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:40
% EndTime: 2019-12-31 16:51:40
% DurationCPUTime: 0.10s
% Computational Cost: add. (52->21), mult. (106->32), div. (0->0), fcn. (108->6), ass. (0->12)
t31 = sin(qJ(3));
t32 = sin(qJ(1));
t33 = cos(qJ(3));
t34 = cos(qJ(1));
t18 = -t32 * t31 - t34 * t33;
t19 = t34 * t31 - t32 * t33;
t25 = sin(qJ(4));
t26 = cos(qJ(4));
t27 = g(1) * t18 + g(2) * t19;
t37 = (t26 * MDP(15) - t25 * MDP(16) + MDP(8)) * (g(1) * t19 - g(2) * t18) - t27 * MDP(9);
t20 = g(1) * t32 - g(2) * t34;
t1 = [(-g(1) * (-t32 * pkin(1) + t34 * qJ(2)) - g(2) * (t34 * pkin(1) + t32 * qJ(2))) * MDP(6) + (MDP(3) - MDP(5)) * (g(1) * t34 + g(2) * t32) + (MDP(2) + MDP(4)) * t20 - t37; -t20 * MDP(6); t37; (g(3) * t26 - t27 * t25) * MDP(15) + (-g(3) * t25 - t27 * t26) * MDP(16);];
taug = t1;
