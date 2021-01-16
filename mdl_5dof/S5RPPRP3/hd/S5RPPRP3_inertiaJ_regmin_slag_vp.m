% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 17:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 17:04:41
% EndTime: 2021-01-15 17:04:42
% DurationCPUTime: 0.13s
% Computational Cost: add. (72->26), mult. (98->37), div. (0->0), fcn. (88->4), ass. (0->19)
t14 = sin(qJ(4));
t20 = 0.2e1 * t14;
t15 = cos(qJ(4));
t19 = 0.2e1 * t15;
t18 = t14 * pkin(4);
t17 = t15 * pkin(4);
t13 = cos(pkin(7));
t10 = -t13 * pkin(1) - pkin(2);
t12 = sin(pkin(7));
t8 = t12 * pkin(1) + qJ(3);
t11 = t15 ^ 2;
t7 = -pkin(6) + t10;
t6 = -t14 ^ 2 - t11;
t5 = t15 * t7;
t4 = t8 + t18;
t3 = -t15 * qJ(5) + t5;
t2 = (-qJ(5) + t7) * t14;
t1 = t2 * t14 + t3 * t15;
t9 = [1, 0, 0, (t12 ^ 2 + t13 ^ 2) * pkin(1) ^ 2, 0.2e1 * t10, 0.2e1 * t8, t10 ^ 2 + t8 ^ 2, t11, -0.2e1 * t15 * t14, 0, 0, 0, t8 * t20, t8 * t19, t4 * t20, t4 * t19, -0.2e1 * t1, t2 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t14 + t2 * t15; 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 1, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, t5, -t14 * t7, t3, -t2, -t17, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, -t14, -t15, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, t15, -t14, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t15, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t9;
