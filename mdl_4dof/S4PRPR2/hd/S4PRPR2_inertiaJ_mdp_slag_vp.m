% Calculate joint inertia matrix for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PRPR2_inertiaJ_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:49
% EndTime: 2019-03-08 18:21:50
% DurationCPUTime: 0.04s
% Computational Cost: add. (39->21), mult. (77->34), div. (0->0), fcn. (78->6), ass. (0->14)
t23 = sin(pkin(6));
t33 = pkin(2) * t23;
t24 = cos(pkin(6));
t26 = sin(qJ(2));
t28 = cos(qJ(2));
t20 = -t23 * t26 + t24 * t28;
t21 = t23 * t28 + t24 * t26;
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t32 = (t27 * t20 - t25 * t21) * MDP(7) + (-t25 * t20 - t27 * t21) * MDP(8);
t22 = t24 * pkin(2) + pkin(3);
t31 = (t27 * t22 - t25 * t33) * MDP(7);
t30 = (-t25 * t22 - t27 * t33) * MDP(8);
t1 = [MDP(1) + (t20 ^ 2 + t21 ^ 2) * MDP(5); t28 * MDP(3) - t26 * MDP(4) + (t20 * t24 + t21 * t23) * MDP(5) * pkin(2) + t32; MDP(2) + MDP(6) + (t23 ^ 2 + t24 ^ 2) * MDP(5) * pkin(2) ^ 2 + 0.2e1 * t31 + 0.2e1 * t30; 0; 0; MDP(5); t32; MDP(6) + t30 + t31; 0; MDP(6);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
