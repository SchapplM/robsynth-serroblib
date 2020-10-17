% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:12
% EndTime: 2019-12-31 18:11:12
% DurationCPUTime: 0.16s
% Computational Cost: add. (92->31), mult. (146->53), div. (0->0), fcn. (123->4), ass. (0->25)
t16 = sin(qJ(3));
t29 = 0.2e1 * t16;
t17 = cos(qJ(3));
t28 = -0.2e1 * t17;
t27 = 0.2e1 * t17;
t14 = sin(pkin(7));
t8 = t14 * pkin(1) + pkin(6);
t26 = t16 * t8;
t10 = t16 * qJ(4);
t25 = t17 * pkin(3) + t10;
t12 = t16 ^ 2;
t6 = t17 ^ 2 + t12;
t24 = t17 * qJ(4);
t15 = cos(pkin(7));
t9 = -t15 * pkin(1) - pkin(2);
t23 = -t16 * pkin(3) + t24;
t2 = t9 - t25;
t21 = qJ(4) ^ 2;
t20 = 0.2e1 * qJ(4);
t18 = pkin(3) + pkin(4);
t5 = t17 * t8;
t4 = -t17 * qJ(5) + t5;
t3 = (-qJ(5) + t8) * t16;
t1 = t17 * pkin(4) - t2;
t7 = [1, 0, 0, (t14 ^ 2 + t15 ^ 2) * pkin(1) ^ 2, t12, t16 * t27, 0, 0, 0, t9 * t28, t9 * t29, t2 * t28, 0.2e1 * t6 * t8, -0.2e1 * t2 * t16, t6 * t8 ^ 2 + t2 ^ 2, t1 * t27, t1 * t29, -0.2e1 * t3 * t16 - 0.2e1 * t4 * t17, t1 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t16 - t3 * t17; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, t16, t17, 0, -t26, -t5, -t26, t23, t5, t23 * t8, -t3, t4, t18 * t16 - t24, t4 * qJ(4) - t3 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t17, 0, t16, t25, t17, t16, 0, t17 * t18 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, t20, pkin(3) ^ 2 + t21, 0.2e1 * t18, t20, 0, t18 ^ 2 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, t26, 0, 0, -t16, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), -1, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t16, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t7;
