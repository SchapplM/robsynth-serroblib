% Calculate Gravitation load on the joints for
% S4RPRP5
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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRP5_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S4RPRP5_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:01
% EndTime: 2019-12-31 16:45:02
% DurationCPUTime: 0.14s
% Computational Cost: add. (90->32), mult. (116->40), div. (0->0), fcn. (89->6), ass. (0->15)
t36 = sin(qJ(1));
t37 = cos(qJ(1));
t27 = g(1) * t37 + g(2) * t36;
t44 = MDP(13) + MDP(15);
t43 = MDP(14) - MDP(17);
t26 = g(1) * t36 - g(2) * t37;
t32 = pkin(6) + qJ(3);
t29 = sin(t32);
t30 = cos(t32);
t40 = t30 * pkin(3) + t29 * qJ(4);
t34 = cos(pkin(6));
t38 = t34 * pkin(2) + pkin(1) + t40;
t35 = -pkin(5) - qJ(2);
t22 = -g(3) * t30 + t27 * t29;
t1 = [(-g(1) * (-t36 * pkin(1) + t37 * qJ(2)) - g(2) * (t37 * pkin(1) + t36 * qJ(2))) * MDP(7) + ((g(1) * t35 - g(2) * t38) * t37 + (g(1) * t38 + g(2) * t35) * t36) * MDP(18) + (MDP(3) - MDP(6) - MDP(16)) * t27 + (MDP(4) * t34 - MDP(5) * sin(pkin(6)) - t43 * t29 + t44 * t30 + MDP(2)) * t26; (-MDP(18) - MDP(7)) * t26; (-g(3) * t40 + t27 * (pkin(3) * t29 - qJ(4) * t30)) * MDP(18) + t43 * (g(3) * t29 + t27 * t30) + t44 * t22; -t22 * MDP(18);];
taug = t1;
