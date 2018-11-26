% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2018-11-23 15:50
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPPRRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:50:12
% EndTime: 2018-11-23 15:50:12
% DurationCPUTime: 0.10s
% Computational Cost: add. (98->51), mult. (97->48), div. (0->0), fcn. (143->8), ass. (0->33)
t12 = sin(qJ(5));
t36 = pkin(5) * t12;
t13 = sin(qJ(4));
t14 = sin(qJ(1));
t35 = t14 * t13;
t15 = cos(qJ(5));
t34 = t14 * t15;
t16 = cos(qJ(4));
t33 = t14 * t16;
t17 = cos(qJ(1));
t32 = t17 * t13;
t31 = t17 * t15;
t30 = t17 * t16;
t10 = pkin(6) + 0;
t29 = t14 * pkin(1) + 0;
t28 = pkin(2) + t10;
t27 = t17 * pkin(1) + t14 * qJ(2) + 0;
t26 = pkin(3) + t28;
t25 = t17 * qJ(3) + t27;
t24 = pkin(4) * t13 - pkin(8) * t16;
t1 = t15 * pkin(5) + pkin(4);
t18 = -pkin(9) - pkin(8);
t23 = t1 * t13 + t16 * t18;
t22 = -t17 * qJ(2) + t29;
t4 = t14 * qJ(3);
t21 = t22 + t4;
t20 = -t14 * pkin(7) + t25;
t8 = t17 * pkin(7);
t19 = t21 + t8;
t11 = qJ(5) + qJ(6);
t3 = cos(t11);
t2 = sin(t11);
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t14, 0, 0; t14, t17, 0, 0; 0, 0, 1, t10; 0, 0, 0, 1; 0, -t17, t14, t27; 0, -t14, -t17, t22; 1, 0, 0, t10; 0, 0, 0, 1; 0, t14, t17, t25; 0, -t17, t14, t21; 1, 0, 0, t28; 0, 0, 0, 1; t32, t30, -t14, t20; t35, t33, t17, t19; t16, -t13, 0, t26; 0, 0, 0, 1; -t14 * t12 + t13 * t31, -t12 * t32 - t34, -t30, t24 * t17 + t20; t17 * t12 + t13 * t34, -t12 * t35 + t31, -t33, t24 * t14 + t19; t16 * t15, -t16 * t12, t13, t16 * pkin(4) + t13 * pkin(8) + t26; 0, 0, 0, 1; -t14 * t2 + t3 * t32, -t14 * t3 - t2 * t32, -t30, t23 * t17 + (-pkin(7) - t36) * t14 + t25; t17 * t2 + t3 * t35, t17 * t3 - t2 * t35, -t33, t4 + t8 + (-qJ(2) + t36) * t17 + t23 * t14 + t29; t16 * t3, -t16 * t2, t13, t16 * t1 - t13 * t18 + t26; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
