% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:43
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPPRPR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:43:26
% EndTime: 2018-11-23 15:43:26
% DurationCPUTime: 0.10s
% Computational Cost: add. (124->52), mult. (98->45), div. (0->0), fcn. (143->8), ass. (0->35)
t19 = sin(qJ(1));
t13 = pkin(9) + qJ(4);
t7 = sin(t13);
t3 = t19 * t7;
t8 = cos(t13);
t41 = t19 * t8;
t21 = cos(qJ(1));
t40 = t21 * t7;
t39 = t21 * t8;
t17 = -pkin(7) - qJ(3);
t38 = pkin(5) - t17;
t37 = qJ(5) * t8;
t15 = sin(pkin(9));
t36 = t19 * t15;
t18 = sin(qJ(6));
t35 = t19 * t18;
t20 = cos(qJ(6));
t34 = t19 * t20;
t33 = t21 * t18;
t32 = t21 * t20;
t14 = pkin(6) + 0;
t31 = t19 * pkin(1) + 0;
t30 = pkin(2) + t14;
t29 = t21 * pkin(1) + t19 * qJ(2) + 0;
t28 = -pkin(3) * t15 - qJ(2);
t16 = cos(pkin(9));
t27 = t16 * pkin(3) + t30;
t26 = pkin(3) * t36 + t29;
t25 = -t19 * t17 + t31;
t24 = -t21 * qJ(2) + t31;
t23 = t8 * pkin(4) + t7 * qJ(5) + t27;
t22 = -t21 * t17 + t26;
t2 = pkin(4) * t3;
t1 = t21 * t37;
t4 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t19, 0, 0; t19, t21, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; 0, -t21, t19, t29; 0, -t19, -t21, t24; 1, 0, 0, t14; 0, 0, 0, 1; t36, t19 * t16, t21, t21 * qJ(3) + t29; -t21 * t15, -t21 * t16, t19, t19 * qJ(3) + t24; t16, -t15, 0, t30; 0, 0, 0, 1; t3, t41, t21, t22; -t40, -t39, t19, t28 * t21 + t25; t8, -t7, 0, t27; 0, 0, 0, 1; t21, -t3, -t41, -t19 * t37 + t2 + t22; t19, t40, t39, t1 + (-pkin(4) * t7 + t28) * t21 + t25; 0, -t8, t7, t23; 0, 0, 0, 1; -t8 * t35 + t32, -t8 * t34 - t33, t3, t2 + t38 * t21 + (pkin(8) * t7 - t37) * t19 + t26; t8 * t33 + t34, t8 * t32 - t35, -t40, t1 + t38 * t19 + ((-pkin(4) - pkin(8)) * t7 + t28) * t21 + t31; t7 * t18, t7 * t20, t8, t8 * pkin(8) + t23; 0, 0, 0, 1;];
T_ges = t4;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
