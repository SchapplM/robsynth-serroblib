% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2018-11-23 16:44
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPPPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:44:30
% EndTime: 2018-11-23 16:44:30
% DurationCPUTime: 0.15s
% Computational Cost: add. (141->60), mult. (259->57), div. (0->0), fcn. (351->8), ass. (0->37)
t27 = sin(pkin(9));
t30 = sin(qJ(2));
t16 = t30 * t27;
t28 = cos(pkin(9));
t17 = t30 * t28;
t53 = pkin(3) * t17 + qJ(4) * t16;
t34 = cos(qJ(1));
t31 = sin(qJ(1));
t33 = cos(qJ(2));
t51 = t31 * t33;
t3 = t27 * t51 + t34 * t28;
t52 = t3 * qJ(4);
t20 = t31 * t30;
t21 = t34 * t30;
t50 = t34 * t33;
t5 = t27 * t50 - t31 * t28;
t49 = t5 * qJ(4);
t48 = -pkin(4) - qJ(3);
t47 = pkin(5) + qJ(4);
t46 = qJ(3) * t30;
t26 = pkin(6) + 0;
t45 = t30 * pkin(2) + t26;
t44 = t34 * pkin(1) + t31 * pkin(7) + 0;
t43 = t31 * pkin(1) - t34 * pkin(7) + 0;
t42 = pkin(2) * t50 + t34 * t46 + t44;
t41 = qJ(5) * t17 + t45 + t53;
t40 = -t33 * qJ(3) + t45;
t6 = t31 * t27 + t28 * t50;
t39 = t6 * pkin(3) + t42;
t38 = pkin(2) * t51 + t31 * t46 + t43;
t4 = -t34 * t27 + t28 * t51;
t37 = t4 * pkin(3) + t38;
t36 = pkin(4) * t21 + t6 * qJ(5) + t39;
t35 = pkin(4) * t20 + t4 * qJ(5) + t37;
t32 = cos(qJ(6));
t29 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t31, 0, 0; t31, t34, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t50, -t21, t31, t44; t51, -t20, -t34, t43; t30, t33, 0, t26; 0, 0, 0, 1; t6, -t5, t21, t42; t4, -t3, t20, t38; t17, -t16, -t33, t40; 0, 0, 0, 1; t21, -t6, t5, t39 + t49; t20, -t4, t3, t37 + t52; -t33, -t17, t16, t40 + t53; 0, 0, 0, 1; t5, -t21, t6, t36 + t49; t3, -t20, t4, t35 + t52; t16, t33, t17, t48 * t33 + t41; 0, 0, 0, 1; t6 * t29 + t5 * t32, -t5 * t29 + t6 * t32, t21, pkin(8) * t21 + t47 * t5 + t36; t4 * t29 + t3 * t32, -t3 * t29 + t4 * t32, t20, pkin(8) * t20 + t47 * t3 + t35; (t27 * t32 + t28 * t29) * t30 (-t27 * t29 + t28 * t32) * t30, -t33, pkin(5) * t16 + (-pkin(8) + t48) * t33 + t41; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
