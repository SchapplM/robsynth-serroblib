% Calculate inertial parameters regressor of gravitation load for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% 
% Output:
% taug_reg [3x(3*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug_reg = S3RPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_gravloadJ_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_gravloadJ_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From joint_gravload_fixb_regressor_matlab.m
t12 = sin(qJ(1));
t13 = cos(qJ(1));
t18 = t13 * pkin(1) + t12 * qJ(2);
t17 = cos(qJ(3));
t16 = sin(qJ(3));
t1 = -t12 * t16 - t13 * t17;
t2 = -t12 * t17 + t13 * t16;
t15 = g(1) * t2 - g(2) * t1;
t14 = g(1) * t1 + g(2) * t2;
t9 = t13 * qJ(2);
t4 = g(1) * t13 + g(2) * t12;
t3 = g(1) * t12 - g(2) * t13;
t5 = [0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, -t4, -g(1) * (-t12 * pkin(1) + t9) - g(2) * t18, 0, 0, 0, 0, 0, 0, -t15, t14, 0, -g(1) * (t9 + (-pkin(1) - pkin(2)) * t12) - g(2) * (t13 * pkin(2) + t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, 0;];
taug_reg  = t5;
