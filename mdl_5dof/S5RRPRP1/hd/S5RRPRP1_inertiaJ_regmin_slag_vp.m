% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:59
% EndTime: 2022-01-20 10:20:00
% DurationCPUTime: 0.22s
% Computational Cost: add. (184->49), mult. (301->72), div. (0->0), fcn. (278->6), ass. (0->32)
t24 = sin(qJ(4));
t37 = 0.2e1 * t24;
t26 = cos(qJ(4));
t36 = -0.2e1 * t26;
t35 = t24 * pkin(4);
t34 = sin(qJ(2)) * pkin(1);
t33 = t26 * pkin(4);
t23 = cos(pkin(8));
t18 = -t23 * pkin(2) - pkin(3);
t12 = t18 - t33;
t20 = cos(qJ(2)) * pkin(1);
t19 = t20 + pkin(2);
t22 = sin(pkin(8));
t8 = t23 * t19 - t22 * t34;
t5 = -pkin(3) - t8;
t4 = t5 - t33;
t32 = t12 + t4;
t31 = t18 + t5;
t9 = t22 * t19 + t23 * t34;
t6 = pkin(7) + t9;
t30 = qJ(5) + t6;
t17 = t22 * pkin(2) + pkin(7);
t29 = qJ(5) + t17;
t21 = t24 ^ 2;
t16 = t26 * t37;
t11 = t29 * t26;
t10 = t29 * t24;
t7 = t11 * t26;
t3 = t30 * t26;
t2 = t30 * t24;
t1 = t3 * t26;
t13 = [1, 0, 0, 1, 0.2e1 * t20, -0.2e1 * t34, t8 ^ 2 + t9 ^ 2, t21, t16, 0, 0, 0, t5 * t36, t5 * t37, t4 * t36, t4 * t37, 0.2e1 * t2 * t24 + 0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 1, t20, -t34, (t22 * t9 + t23 * t8) * pkin(2), t21, t16, 0, 0, 0, -t31 * t26, t31 * t24, -t32 * t26, t32 * t24, t1 + t7 + (t10 + t2) * t24, t2 * t10 + t3 * t11 + t4 * t12; 0, 0, 0, 1, 0, 0, (t22 ^ 2 + t23 ^ 2) * pkin(2) ^ 2, t21, t16, 0, 0, 0, t18 * t36, t18 * t37, t12 * t36, t12 * t37, 0.2e1 * t10 * t24 + 0.2e1 * t7, t10 ^ 2 + t11 ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t26 + t3 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t26 + t11 * t24; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 ^ 2 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t26, 0, -t24 * t6, -t26 * t6, -t2, -t3, -t35, -t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t26, 0, -t24 * t17, -t26 * t17, -t10, -t11, -t35, -t10 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t24, t26, -t24, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t24, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t24, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t13;
