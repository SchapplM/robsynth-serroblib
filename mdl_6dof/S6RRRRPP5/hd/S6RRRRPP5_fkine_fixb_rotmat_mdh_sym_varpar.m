% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPP5
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
% Datum: 2018-11-23 18:07
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRRPP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:07:02
% EndTime: 2018-11-23 18:07:02
% DurationCPUTime: 0.15s
% Computational Cost: add. (184->58), mult. (204->57), div. (0->0), fcn. (278->8), ass. (0->37)
t26 = qJ(3) + qJ(4);
t20 = sin(t26);
t21 = cos(t26);
t32 = cos(qJ(1));
t29 = sin(qJ(1));
t31 = cos(qJ(2));
t48 = t29 * t31;
t3 = t20 * t48 + t32 * t21;
t4 = -t32 * t20 + t21 * t48;
t51 = t4 * pkin(4) + t3 * qJ(5);
t46 = t32 * t31;
t5 = t20 * t46 - t29 * t21;
t6 = t29 * t20 + t21 * t46;
t50 = t6 * pkin(4) + t5 * qJ(5);
t28 = sin(qJ(2));
t12 = t28 * t20;
t13 = t28 * t21;
t27 = sin(qJ(3));
t49 = t29 * t27;
t16 = t29 * t28;
t17 = t32 * t28;
t25 = pkin(6) + 0;
t44 = t29 * pkin(1) + 0;
t33 = -pkin(9) - pkin(8);
t43 = t28 * (-qJ(6) - t33);
t42 = t32 * pkin(1) + t29 * pkin(7) + 0;
t30 = cos(qJ(3));
t18 = t30 * pkin(3) + pkin(2);
t41 = t28 * t18 + t31 * t33 + t25;
t40 = pkin(2) * t31 + pkin(8) * t28;
t39 = -t32 * pkin(7) + t44;
t38 = pkin(3) * t49 + t18 * t46 + t42;
t37 = pkin(4) * t13 + qJ(5) * t12 + t41;
t36 = t18 * t48 + (-pkin(3) * t27 - pkin(7)) * t32 + t44;
t35 = -t33 * t17 + t38;
t34 = -t33 * t16 + t36;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t32, -t29, 0, 0; t29, t32, 0, 0; 0, 0, 1, t25; 0, 0, 0, 1; t46, -t17, t29, t42; t48, -t16, -t32, t39; t28, t31, 0, t25; 0, 0, 0, 1; t30 * t46 + t49, -t27 * t46 + t29 * t30, t17, t40 * t32 + t42; -t32 * t27 + t30 * t48, -t27 * t48 - t32 * t30, t16, t40 * t29 + t39; t28 * t30, -t28 * t27, -t31, t28 * pkin(2) - t31 * pkin(8) + t25; 0, 0, 0, 1; t6, -t5, t17, t35; t4, -t3, t16, t34; t13, -t12, -t31, t41; 0, 0, 0, 1; t6, t17, t5, t35 + t50; t4, t16, t3, t34 + t51; t13, -t31, t12, t37; 0, 0, 0, 1; t6, t5, -t17, t6 * pkin(5) + t32 * t43 + t38 + t50; t4, t3, -t16, t4 * pkin(5) + t29 * t43 + t36 + t51; t13, t12, t31, pkin(5) * t13 + t31 * qJ(6) + t37; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
