% Calculate inertial parameters regressor of gravitation load for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
% 
% Output:
% taug_reg [3x(3*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug_reg = S3RPP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_gravloadJ_reg2_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_gravloadJ_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From joint_gravload_fixb_regressor_matlab.m
t7 = sin(qJ(1));
t8 = cos(qJ(1));
t9 = t8 * pkin(1) + t7 * qJ(2);
t4 = t8 * qJ(2);
t2 = g(1) * t8 + g(2) * t7;
t1 = g(1) * t7 - g(2) * t8;
t3 = [0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (-t7 * pkin(1) + t4) - g(2) * t9, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -g(1) * (t4 + (-pkin(1) - qJ(3)) * t7) - g(2) * (t8 * qJ(3) + t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2;];
taug_reg  = t3;
