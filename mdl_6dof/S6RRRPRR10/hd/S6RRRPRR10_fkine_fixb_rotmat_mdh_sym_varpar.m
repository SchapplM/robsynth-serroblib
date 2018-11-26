% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2018-11-23 17:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRPRR10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:58:15
% EndTime: 2018-11-23 17:58:15
% DurationCPUTime: 0.14s
% Computational Cost: add. (166->68), mult. (286->76), div. (0->0), fcn. (388->10), ass. (0->39)
t28 = sin(qJ(2));
t31 = cos(qJ(3));
t13 = t28 * t31;
t27 = sin(qJ(3));
t49 = t28 * t27;
t51 = pkin(3) * t13 + qJ(4) * t49;
t26 = sin(qJ(5));
t50 = t26 * t27;
t29 = sin(qJ(1));
t14 = t29 * t28;
t32 = cos(qJ(2));
t48 = t29 * t32;
t33 = cos(qJ(1));
t16 = t33 * t28;
t47 = t33 * t32;
t24 = pkin(6) + 0;
t46 = t28 * pkin(2) + t24;
t45 = pkin(5) * t26 + qJ(4);
t44 = t33 * pkin(1) + t29 * pkin(7) + 0;
t43 = t46 + t51;
t42 = t29 * pkin(1) - t33 * pkin(7) + 0;
t41 = pkin(2) * t47 + pkin(8) * t16 + t44;
t40 = -t32 * pkin(8) + t46;
t6 = t29 * t27 + t31 * t47;
t39 = t6 * pkin(3) + t41;
t38 = pkin(2) * t48 + pkin(8) * t14 + t42;
t4 = -t33 * t27 + t31 * t48;
t37 = t4 * pkin(3) + t38;
t5 = t27 * t47 - t29 * t31;
t36 = t5 * qJ(4) + t39;
t3 = t27 * t48 + t33 * t31;
t35 = t3 * qJ(4) + t37;
t34 = -pkin(10) - pkin(9);
t30 = cos(qJ(5));
t25 = qJ(5) + qJ(6);
t19 = cos(t25);
t18 = sin(t25);
t17 = t30 * pkin(5) + pkin(4);
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t33, -t29, 0, 0; t29, t33, 0, 0; 0, 0, 1, t24; 0, 0, 0, 1; t47, -t16, t29, t44; t48, -t14, -t33, t42; t28, t32, 0, t24; 0, 0, 0, 1; t6, -t5, t16, t41; t4, -t3, t14, t38; t13, -t49, -t32, t40; 0, 0, 0, 1; t6, t16, t5, t36; t4, t14, t3, t35; t13, -t32, t49, t40 + t51; 0, 0, 0, 1; t5 * t26 + t6 * t30, -t6 * t26 + t5 * t30, -t16, t6 * pkin(4) - pkin(9) * t16 + t36; t3 * t26 + t4 * t30, -t4 * t26 + t3 * t30, -t14, t4 * pkin(4) - pkin(9) * t14 + t35; (t30 * t31 + t50) * t28 (-t26 * t31 + t27 * t30) * t28, t32, pkin(4) * t13 + (-pkin(8) + pkin(9)) * t32 + t43; 0, 0, 0, 1; t5 * t18 + t6 * t19, -t6 * t18 + t5 * t19, -t16, t34 * t16 + t6 * t17 + t45 * t5 + t39; t3 * t18 + t4 * t19, -t4 * t18 + t3 * t19, -t14, t34 * t14 + t4 * t17 + t45 * t3 + t37; (t18 * t27 + t19 * t31) * t28 (-t18 * t31 + t19 * t27) * t28, t32 (-pkin(8) - t34) * t32 + (pkin(5) * t50 + t17 * t31) * t28 + t43; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
