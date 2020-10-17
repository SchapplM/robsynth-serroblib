% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:58
% EndTime: 2019-12-05 16:41:59
% DurationCPUTime: 0.19s
% Computational Cost: add. (81->25), mult. (159->49), div. (0->0), fcn. (125->4), ass. (0->25)
t14 = sin(qJ(4));
t12 = t14 ^ 2;
t16 = cos(qJ(4));
t20 = t16 ^ 2 + t12;
t25 = sin(qJ(3)) * pkin(2);
t8 = pkin(7) + t25;
t28 = t20 * t8;
t33 = -0.2e1 * t14;
t32 = -0.2e1 * t16;
t31 = 0.2e1 * t16;
t22 = cos(qJ(3)) * pkin(2);
t9 = -pkin(3) - t22;
t30 = pkin(3) - t9;
t19 = t16 * pkin(4) + t14 * qJ(5);
t2 = -pkin(3) - t19;
t1 = t2 - t22;
t29 = -t1 - t2;
t27 = t14 * pkin(7);
t26 = t14 * t8;
t24 = t16 * pkin(7);
t23 = t16 * t8;
t21 = t20 * pkin(7);
t3 = -t14 * pkin(4) + t16 * qJ(5);
t6 = t14 * t31;
t4 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 1, 0.2e1 * t22, -0.2e1 * t25, t12, t6, 0, 0, 0, t9 * t32, 0.2e1 * t9 * t14, t1 * t32, 0.2e1 * t28, t1 * t33, t20 * t8 ^ 2 + t1 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t22, -t25, t12, t6, 0, 0, 0, t30 * t16, -t30 * t14, t29 * t16, t21 + t28, t29 * t14, pkin(7) * t28 + t1 * t2; 0, 0, 0, 0, 1, 0, 0, t12, t6, 0, 0, 0, pkin(3) * t31, pkin(3) * t33, t2 * t32, 0.2e1 * t21, t2 * t33, t20 * pkin(7) ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, t16, 0, t14, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, -t26, -t23, -t26, t3, t23, t3 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, -t27, -t24, -t27, t3, t24, t3 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
