% Calculate Gravitation load on the joints for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PPRR3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PPRR3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:26
% EndTime: 2019-12-31 16:17:27
% DurationCPUTime: 0.04s
% Computational Cost: add. (26->12), mult. (57->22), div. (0->0), fcn. (60->6), ass. (0->10)
t26 = cos(qJ(3));
t25 = sin(qJ(3));
t24 = sin(pkin(6));
t19 = cos(pkin(6));
t14 = -t19 * t26 - t24 * t25;
t15 = t19 * t25 - t24 * t26;
t22 = -g(1) * t14 - g(2) * t15;
t21 = cos(qJ(4));
t20 = sin(qJ(4));
t1 = [(-MDP(1) - MDP(2)) * g(3); (-g(1) * t24 + g(2) * t19) * MDP(2); t22 * MDP(5) + (t21 * MDP(11) - t20 * MDP(12) + MDP(4)) * (g(1) * t15 - g(2) * t14); (g(3) * t21 + t22 * t20) * MDP(11) + (-g(3) * t20 + t22 * t21) * MDP(12);];
taug = t1;
