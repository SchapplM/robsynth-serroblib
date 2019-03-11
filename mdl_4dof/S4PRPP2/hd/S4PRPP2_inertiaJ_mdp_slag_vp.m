% Calculate joint inertia matrix for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PRPP2_inertiaJ_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:00
% EndTime: 2019-03-08 18:19:00
% DurationCPUTime: 0.04s
% Computational Cost: add. (37->21), mult. (65->34), div. (0->0), fcn. (57->4), ass. (0->10)
t29 = MDP(5) + MDP(8);
t28 = cos(qJ(2));
t25 = sin(qJ(2));
t24 = cos(pkin(5));
t23 = sin(pkin(5));
t21 = t24 * pkin(2) + pkin(3);
t19 = t23 * pkin(2) + qJ(4);
t18 = t23 * t28 + t24 * t25;
t16 = t23 * t25 - t24 * t28;
t1 = [MDP(1) + t29 * (t16 ^ 2 + t18 ^ 2); t28 * MDP(3) - t25 * MDP(4) - t16 * MDP(6) + t18 * MDP(7) + (-t16 * t21 + t18 * t19) * MDP(8) + (-t16 * t24 + t18 * t23) * MDP(5) * pkin(2); MDP(2) + (t19 ^ 2 + t21 ^ 2) * MDP(8) + (t23 ^ 2 + t24 ^ 2) * MDP(5) * pkin(2) ^ 2 + 0.2e1 * t21 * MDP(6) + 0.2e1 * t19 * MDP(7); 0; 0; t29; t16 * MDP(8); -t21 * MDP(8) - MDP(6); 0; MDP(8);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
