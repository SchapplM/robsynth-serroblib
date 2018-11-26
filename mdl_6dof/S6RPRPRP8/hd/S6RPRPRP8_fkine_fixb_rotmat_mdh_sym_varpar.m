% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2018-11-23 16:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRPRP8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:01:16
% EndTime: 2018-11-23 16:01:16
% DurationCPUTime: 0.12s
% Computational Cost: add. (138->48), mult. (123->44), div. (0->0), fcn. (177->8), ass. (0->37)
t19 = qJ(3) + pkin(9);
t14 = cos(t19);
t22 = sin(qJ(5));
t46 = t14 * t22;
t13 = sin(t19);
t24 = sin(qJ(1));
t45 = t24 * t13;
t44 = t24 * t14;
t43 = t24 * t22;
t23 = sin(qJ(3));
t42 = t24 * t23;
t25 = cos(qJ(5));
t41 = t24 * t25;
t27 = cos(qJ(1));
t8 = t27 * t14;
t40 = t27 * t22;
t39 = t27 * t25;
t20 = pkin(6) + 0;
t38 = t24 * pkin(1) + 0;
t37 = pkin(2) + t20;
t36 = -pkin(3) * t23 - qJ(2);
t35 = t27 * pkin(1) + t24 * qJ(2) + 0;
t26 = cos(qJ(3));
t34 = t26 * pkin(3) + t37;
t21 = -qJ(4) - pkin(7);
t33 = -t24 * t21 + t38;
t32 = -t27 * qJ(2) + t38;
t31 = t14 * pkin(4) + t13 * pkin(8) + t34;
t30 = pkin(3) * t42 - t27 * t21 + t35;
t29 = pkin(4) * t45 - pkin(8) * t44 + t30;
t28 = pkin(8) * t8 + (-pkin(4) * t13 + t36) * t27 + t33;
t7 = t14 * t25;
t4 = -t13 * t39 + t43;
t3 = t13 * t40 + t41;
t2 = t13 * t41 + t40;
t1 = t13 * t43 - t39;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t27, -t24, 0, 0; t24, t27, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; 0, -t27, t24, t35; 0, -t24, -t27, t32; 1, 0, 0, t20; 0, 0, 0, 1; t42, t24 * t26, t27, t27 * pkin(7) + t35; -t27 * t23, -t27 * t26, t24, t24 * pkin(7) + t32; t26, -t23, 0, t37; 0, 0, 0, 1; t45, t44, t27, t30; -t27 * t13, -t8, t24, t27 * t36 + t33; t14, -t13, 0, t34; 0, 0, 0, 1; t2, -t1, -t44, t29; t4, t3, t8, t28; t7, -t46, t13, t31; 0, 0, 0, 1; t2, -t44, t1, t2 * pkin(5) + t1 * qJ(6) + t29; t4, t8, -t3, t4 * pkin(5) - t3 * qJ(6) + t28; t7, t13, t46 (pkin(5) * t25 + qJ(6) * t22) * t14 + t31; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
