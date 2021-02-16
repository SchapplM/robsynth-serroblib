% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 14:56:28
% EndTime: 2021-01-15 14:56:29
% DurationCPUTime: 0.12s
% Computational Cost: add. (42->28), mult. (82->39), div. (0->0), fcn. (87->4), ass. (0->18)
t10 = cos(qJ(4));
t19 = 0.2e1 * t10;
t8 = sin(qJ(4));
t9 = sin(qJ(3));
t18 = t8 * t9;
t5 = t8 ^ 2;
t17 = t10 ^ 2 + t5;
t16 = t10 * pkin(4);
t15 = t10 * t9;
t11 = cos(qJ(3));
t14 = t11 * t8;
t13 = qJ(5) + pkin(6);
t1 = t13 * t8;
t2 = t13 * t10;
t12 = t1 * t8 + t2 * t10;
t4 = -pkin(3) - t16;
t3 = t11 * t10;
t6 = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t9 ^ 2 + t11 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t1 - t8 * t2; 0, 0, 0, t11, -t9, 0, 0, 0, 0, 0, t3, -t14, t3, -t14, t17 * t9, -t11 * t4 + t12 * t9; 0, 0, 1, 0, 0, t5, t8 * t19, 0, 0, 0, pkin(3) * t19, -0.2e1 * pkin(3) * t8, -0.2e1 * t4 * t10, 0.2e1 * t4 * t8, 0.2e1 * t12, t1 ^ 2 + t2 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t8, -t10, t8, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t15, -t18, -t15, 0, -pkin(4) * t18; 0, 0, 0, 0, 0, 0, 0, t8, t10, 0, -t8 * pkin(6), -t10 * pkin(6), -t1, -t2, -t8 * pkin(4), -t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t8, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
