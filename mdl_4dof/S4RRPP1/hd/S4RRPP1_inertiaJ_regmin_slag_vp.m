% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t13 = cos(qJ(2)) * pkin(1);
t11 = t13 + pkin(2);
t14 = sin(pkin(6));
t15 = cos(pkin(6));
t20 = sin(qJ(2)) * pkin(1);
t4 = t14 * t11 + t15 * t20;
t21 = t15 * pkin(2);
t19 = -t15 * t11 + t14 * t20;
t12 = t14 * pkin(2);
t9 = pkin(3) + t21;
t8 = t12 + qJ(4);
t2 = -pkin(3) + t19;
t1 = qJ(4) + t4;
t3 = [1, 0, 0, 1, 0.2e1 * t13, -0.2e1 * t20, t19 ^ 2 + t4 ^ 2, -0.2e1 * t2, 0.2e1 * t1, t1 ^ 2 + t2 ^ 2; 0, 0, 0, 1, t13, -t20 (t14 * t4 - t15 * t19) * pkin(2), 0.2e1 * pkin(3) - t19 + t21, t12 + 0.2e1 * qJ(4) + t4, t1 * t8 - t2 * t9; 0, 0, 0, 1, 0, 0 (t14 ^ 2 + t15 ^ 2) * pkin(2) ^ 2, 0.2e1 * t9, 0.2e1 * t8, t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, -1, 0, t2; 0, 0, 0, 0, 0, 0, 0, -1, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
