% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:36
% EndTime: 2019-12-31 19:05:36
% DurationCPUTime: 0.08s
% Computational Cost: add. (83->37), mult. (89->28), div. (0->0), fcn. (141->8), ass. (0->20)
t27 = cos(qJ(3));
t26 = sin(qJ(1));
t25 = sin(qJ(3));
t15 = pkin(5) + 0;
t10 = -pkin(6) + t15;
t19 = cos(qJ(1));
t24 = t19 * pkin(1) + t26 * qJ(2) + 0;
t23 = t19 * pkin(2) + t24;
t22 = t26 * pkin(1) - t19 * qJ(2) + 0;
t21 = t26 * pkin(2) + t22;
t20 = -pkin(8) - pkin(7);
t18 = cos(qJ(4));
t17 = sin(qJ(4));
t16 = qJ(4) + qJ(5);
t8 = cos(t16);
t7 = sin(t16);
t6 = t18 * pkin(4) + pkin(3);
t2 = t19 * t25 - t26 * t27;
t1 = -t19 * t27 - t26 * t25;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t26, 0, 0; t26, t19, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t19, 0, t26, t24; t26, 0, -t19, t22; 0, 1, 0, t15; 0, 0, 0, 1; -t1, -t2, 0, t23; -t2, t1, 0, t21; 0, 0, -1, t10; 0, 0, 0, 1; -t1 * t18, t1 * t17, t2, -t1 * pkin(3) + t2 * pkin(7) + t23; -t2 * t18, t2 * t17, -t1, -t2 * pkin(3) - t1 * pkin(7) + t21; -t17, -t18, 0, t10; 0, 0, 0, 1; -t1 * t8, t1 * t7, t2, -t1 * t6 - t2 * t20 + t23; -t2 * t8, t2 * t7, -t1, t1 * t20 - t2 * t6 + t21; -t7, -t8, 0, -t17 * pkin(4) + t10; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
