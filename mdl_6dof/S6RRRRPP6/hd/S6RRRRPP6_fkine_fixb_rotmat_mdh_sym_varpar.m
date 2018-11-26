% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2018-11-23 18:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRPP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:07:50
% EndTime: 2018-11-23 18:07:50
% DurationCPUTime: 0.15s
% Computational Cost: add. (184->60), mult. (204->57), div. (0->0), fcn. (278->8), ass. (0->37)
t27 = qJ(3) + qJ(4);
t21 = sin(t27);
t22 = cos(t27);
t33 = cos(qJ(1));
t30 = sin(qJ(1));
t32 = cos(qJ(2));
t48 = t30 * t32;
t3 = t21 * t48 + t33 * t22;
t4 = -t33 * t21 + t22 * t48;
t52 = t4 * pkin(4) + t3 * qJ(5);
t47 = t33 * t32;
t5 = t21 * t47 - t30 * t22;
t6 = t30 * t21 + t22 * t47;
t51 = t6 * pkin(4) + t5 * qJ(5);
t29 = sin(qJ(2));
t12 = t29 * t21;
t13 = t29 * t22;
t28 = sin(qJ(3));
t49 = t30 * t28;
t17 = t30 * t29;
t18 = t33 * t29;
t26 = pkin(6) + 0;
t45 = t30 * pkin(1) + 0;
t34 = -pkin(9) - pkin(8);
t44 = t29 * (pkin(5) - t34);
t43 = t33 * pkin(1) + t30 * pkin(7) + 0;
t31 = cos(qJ(3));
t19 = t31 * pkin(3) + pkin(2);
t42 = t29 * t19 + t32 * t34 + t26;
t41 = pkin(2) * t32 + pkin(8) * t29;
t40 = -t33 * pkin(7) + t45;
t39 = pkin(3) * t49 + t19 * t47 + t43;
t38 = pkin(4) * t13 + qJ(5) * t12 + t42;
t37 = t19 * t48 + (-pkin(3) * t28 - pkin(7)) * t33 + t45;
t36 = -t34 * t18 + t39;
t35 = -t34 * t17 + t37;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t33, -t30, 0, 0; t30, t33, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t47, -t18, t30, t43; t48, -t17, -t33, t40; t29, t32, 0, t26; 0, 0, 0, 1; t31 * t47 + t49, -t28 * t47 + t30 * t31, t18, t41 * t33 + t43; -t33 * t28 + t31 * t48, -t28 * t48 - t33 * t31, t17, t41 * t30 + t40; t29 * t31, -t29 * t28, -t32, t29 * pkin(2) - t32 * pkin(8) + t26; 0, 0, 0, 1; t6, -t5, t18, t36; t4, -t3, t17, t35; t13, -t12, -t32, t42; 0, 0, 0, 1; t18, -t6, t5, t36 + t51; t17, -t4, t3, t35 + t52; -t32, -t13, t12, t38; 0, 0, 0, 1; t18, t5, t6, t6 * qJ(6) + t33 * t44 + t39 + t51; t17, t3, t4, t4 * qJ(6) + t30 * t44 + t37 + t52; -t32, t12, t13, -t32 * pkin(5) + qJ(6) * t13 + t38; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
