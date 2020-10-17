% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:15
% EndTime: 2019-12-31 17:14:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (77->23), mult. (155->49), div. (0->0), fcn. (116->4), ass. (0->24)
t14 = sin(qJ(3));
t12 = t14 ^ 2;
t16 = cos(qJ(3));
t19 = t16 ^ 2 + t12;
t24 = sin(qJ(2)) * pkin(1);
t8 = pkin(6) + t24;
t27 = t19 * t8;
t32 = -0.2e1 * t14;
t31 = -0.2e1 * t16;
t30 = 0.2e1 * t16;
t21 = cos(qJ(2)) * pkin(1);
t9 = -pkin(2) - t21;
t29 = pkin(2) - t9;
t2 = -t16 * pkin(3) - t14 * qJ(4) - pkin(2);
t1 = t2 - t21;
t28 = -t1 - t2;
t26 = t14 * pkin(6);
t25 = t14 * t8;
t23 = t16 * pkin(6);
t22 = t16 * t8;
t20 = t19 * pkin(6);
t3 = -t14 * pkin(3) + t16 * qJ(4);
t6 = t14 * t30;
t4 = [1, 0, 0, 1, 0.2e1 * t21, -0.2e1 * t24, t12, t6, 0, 0, 0, t9 * t31, 0.2e1 * t9 * t14, t1 * t31, 0.2e1 * t27, t1 * t32, t19 * t8 ^ 2 + t1 ^ 2; 0, 0, 0, 1, t21, -t24, t12, t6, 0, 0, 0, t29 * t16, -t29 * t14, t28 * t16, t20 + t27, t28 * t14, pkin(6) * t27 + t1 * t2; 0, 0, 0, 1, 0, 0, t12, t6, 0, 0, 0, pkin(2) * t30, pkin(2) * t32, t2 * t31, 0.2e1 * t20, t2 * t32, t19 * pkin(6) ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, -t25, -t22, -t25, t3, t22, t3 * t8; 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, -t26, -t23, -t26, t3, t23, t3 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
