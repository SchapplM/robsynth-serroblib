% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:33
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:32:32
% EndTime: 2021-01-15 15:32:34
% DurationCPUTime: 0.22s
% Computational Cost: add. (186->47), mult. (400->80), div. (0->0), fcn. (452->6), ass. (0->31)
t27 = cos(qJ(3));
t21 = -t27 * pkin(3) - pkin(2);
t42 = 0.2e1 * t21;
t41 = 0.2e1 * t27;
t23 = sin(pkin(8));
t40 = t23 * pkin(3);
t24 = cos(pkin(8));
t39 = t24 * pkin(3);
t25 = sin(qJ(3));
t13 = t23 * t25 - t24 * t27;
t28 = cos(qJ(2));
t38 = t28 * t13;
t14 = t23 * t27 + t24 * t25;
t37 = t28 * t14;
t36 = qJ(4) + pkin(6);
t16 = t36 * t27;
t31 = t36 * t25;
t6 = t23 * t16 + t24 * t31;
t8 = t24 * t16 - t23 * t31;
t35 = t6 ^ 2 + t8 ^ 2;
t26 = sin(qJ(2));
t10 = t14 * t26;
t12 = t13 * t26;
t34 = t10 * t6 - t12 * t8;
t33 = t10 * t14 + t12 * t13;
t32 = t10 ^ 2 + t12 ^ 2 + t28 ^ 2;
t30 = -0.2e1 * t8 * t13 + 0.2e1 * t6 * t14;
t19 = pkin(4) + t39;
t17 = qJ(5) + t40;
t3 = t13 * pkin(4) - t14 * qJ(5) + t21;
t1 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, t32; 0, 0, t28, -t26, 0, 0, 0, 0, 0, t28 * t27, -t28 * t25, -t38, -t37, t33, -t28 * t21 + t34, -t38, t33, t37, -t28 * t3 + t34; 0, 1, 0, 0, t25 ^ 2, t25 * t41, 0, 0, 0, pkin(2) * t41, -0.2e1 * pkin(2) * t25, t13 * t42, t14 * t42, t30, t21 ^ 2 + t35, 0.2e1 * t3 * t13, t30, -0.2e1 * t3 * t14, t3 ^ 2 + t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25 * t26, -t27 * t26, -t10, t12, 0, (-t10 * t24 - t12 * t23) * pkin(3), -t10, 0, -t12, -t10 * t19 - t12 * t17; 0, 0, 0, 0, 0, 0, t25, t27, 0, -t25 * pkin(6), -t27 * pkin(6), -t6, -t8, (-t13 * t23 - t14 * t24) * pkin(3), (t23 * t8 - t24 * t6) * pkin(3), -t6, -t17 * t13 - t19 * t14, t8, t8 * t17 - t6 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t39, -0.2e1 * t40, 0, (t23 ^ 2 + t24 ^ 2) * pkin(3) ^ 2, 0.2e1 * t19, 0, 0.2e1 * t17, t17 ^ 2 + t19 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, 0, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, t21, t13, 0, -t14, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
