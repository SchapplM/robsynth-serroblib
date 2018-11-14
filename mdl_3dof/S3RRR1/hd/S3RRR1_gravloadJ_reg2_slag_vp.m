% Calculate inertial parameters regressor of gravitation load for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% taug_reg [3x(3*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug_reg = S3RRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_gravloadJ_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_gravloadJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From joint_gravload_fixb_regressor_matlab.m
t10 = qJ(1) + qJ(2);
t7 = sin(t10);
t8 = cos(t10);
t3 = g(1) * t7 - g(2) * t8;
t11 = sin(qJ(1));
t12 = cos(qJ(1));
t13 = g(1) * t11 - g(2) * t12;
t9 = qJ(3) + t10;
t6 = cos(t9);
t5 = sin(t9);
t4 = g(1) * t8 + g(2) * t7;
t2 = g(1) * t6 + g(2) * t5;
t1 = g(1) * t5 - g(2) * t6;
t14 = [0, 0, 0, 0, 0, 0, t13, g(1) * t12 + g(2) * t11, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t13 * pkin(1), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t11 * pkin(1) - pkin(2) * t7) - g(2) * (t12 * pkin(1) + pkin(2) * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t14;
