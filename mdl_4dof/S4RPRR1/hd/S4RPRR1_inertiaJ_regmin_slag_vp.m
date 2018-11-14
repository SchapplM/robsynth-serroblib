% Calculate minimal parameter regressor of joint inertia matrix for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MM_reg = S4RPRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t10 = sin(pkin(7));
t19 = pkin(1) * t10;
t12 = sin(qJ(4));
t18 = t12 * pkin(3);
t14 = cos(qJ(4));
t13 = sin(qJ(3));
t15 = cos(qJ(3));
t11 = cos(pkin(7));
t8 = t11 * pkin(1) + pkin(2);
t6 = t13 * t8 + t15 * t19;
t17 = t14 * t6;
t5 = -t13 * t19 + t15 * t8;
t4 = pkin(3) + t5;
t1 = -t12 * t6 + t14 * t4;
t9 = t14 * pkin(3);
t2 = -t12 * t4 - t17;
t3 = [1, 0, 0 (t10 ^ 2 + t11 ^ 2) * pkin(1) ^ 2, 1, 0.2e1 * t5, -0.2e1 * t6, 1, 0.2e1 * t1, 0.2e1 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t5, -t6, 1, t1 + t9, -t17 + (-pkin(3) - t4) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t9, -0.2e1 * t18; 0, 0, 0, 0, 0, 0, 0, 1, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, t9, -t18; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
