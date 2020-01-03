% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:25
% EndTime: 2019-12-31 17:59:25
% DurationCPUTime: 0.08s
% Computational Cost: add. (99->33), mult. (54->30), div. (0->0), fcn. (87->8), ass. (0->25)
t16 = cos(qJ(4));
t11 = qJ(1) + pkin(8);
t6 = sin(t11);
t30 = t6 * t16;
t7 = cos(t11);
t29 = t7 * t16;
t12 = sin(qJ(5));
t13 = sin(qJ(4));
t28 = t12 * t13;
t15 = cos(qJ(5));
t27 = t13 * t15;
t26 = pkin(5) + 0;
t14 = sin(qJ(1));
t25 = t14 * pkin(1) + 0;
t17 = cos(qJ(1));
t24 = t17 * pkin(1) + 0;
t23 = t6 * pkin(2) + t25;
t8 = qJ(2) + t26;
t22 = t7 * pkin(2) + t6 * qJ(3) + t24;
t21 = pkin(3) + t8;
t20 = pkin(4) * t13 - pkin(7) * t16;
t19 = t7 * pkin(6) + t22;
t18 = -t7 * qJ(3) + t23;
t2 = t6 * pkin(6);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t14, 0, 0; t14, t17, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t7, -t6, 0, t24; t6, t7, 0, t25; 0, 0, 1, t8; 0, 0, 0, 1; 0, -t7, t6, t22; 0, -t6, -t7, t18; 1, 0, 0, t8; 0, 0, 0, 1; t6 * t13, t30, t7, t19; -t7 * t13, -t29, t6, t18 + t2; t16, -t13, 0, t21; 0, 0, 0, 1; t7 * t12 + t6 * t27, t7 * t15 - t6 * t28, -t30, t20 * t6 + t19; t6 * t12 - t7 * t27, t6 * t15 + t7 * t28, t29, t2 + (-qJ(3) - t20) * t7 + t23; t16 * t15, -t16 * t12, t13, t16 * pkin(4) + t13 * pkin(7) + t21; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
