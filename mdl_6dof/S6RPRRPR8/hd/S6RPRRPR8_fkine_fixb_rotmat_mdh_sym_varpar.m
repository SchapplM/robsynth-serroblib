% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 16:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRPR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:20:09
% EndTime: 2018-11-23 16:20:10
% DurationCPUTime: 0.14s
% Computational Cost: add. (140->63), mult. (129->66), div. (0->0), fcn. (184->10), ass. (0->37)
t23 = cos(qJ(4));
t6 = t23 * pkin(4) + pkin(3);
t20 = sin(qJ(4));
t22 = sin(qJ(1));
t41 = t22 * t20;
t21 = sin(qJ(3));
t40 = t22 * t21;
t39 = t22 * t23;
t24 = cos(qJ(3));
t38 = t22 * t24;
t25 = cos(qJ(1));
t37 = t25 * t20;
t36 = t25 * t21;
t35 = t25 * t23;
t19 = -qJ(5) - pkin(8);
t18 = pkin(6) + 0;
t17 = qJ(4) + pkin(10);
t34 = t22 * pkin(1) + 0;
t33 = pkin(2) + t18;
t11 = t22 * pkin(7);
t32 = t11 + t34;
t31 = t25 * pkin(1) + t22 * qJ(2) + 0;
t30 = t25 * pkin(7) + t31;
t29 = pkin(3) * t21 - pkin(8) * t24;
t8 = cos(t17);
t1 = pkin(5) * t8 + t6;
t16 = -pkin(9) + t19;
t28 = t1 * t21 + t16 * t24;
t27 = t19 * t24 + t21 * t6;
t26 = -t25 * qJ(2) + t34;
t9 = qJ(6) + t17;
t7 = sin(t17);
t5 = t25 * t24;
t4 = cos(t9);
t3 = sin(t9);
t2 = t20 * pkin(4) + pkin(5) * t7;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t25, -t22, 0, 0; t22, t25, 0, 0; 0, 0, 1, t18; 0, 0, 0, 1; 0, -t25, t22, t31; 0, -t22, -t25, t26; 1, 0, 0, t18; 0, 0, 0, 1; t40, t38, t25, t30; -t36, -t5, t22, t11 + t26; t24, -t21, 0, t33; 0, 0, 0, 1; t21 * t39 + t37, -t20 * t40 + t35, -t38, t29 * t22 + t30; -t21 * t35 + t41, t20 * t36 + t39, t5 (-qJ(2) - t29) * t25 + t32; t24 * t23, -t24 * t20, t21, t24 * pkin(3) + t21 * pkin(8) + t33; 0, 0, 0, 1; t25 * t7 + t8 * t40, t25 * t8 - t7 * t40, -t38, pkin(4) * t37 + t27 * t22 + t30; t22 * t7 - t8 * t36, t22 * t8 + t7 * t36, t5, pkin(4) * t41 + (-qJ(2) - t27) * t25 + t32; t24 * t8, -t24 * t7, t21, -t21 * t19 + t24 * t6 + t33; 0, 0, 0, 1; t25 * t3 + t4 * t40, t25 * t4 - t3 * t40, -t38, t25 * t2 + t28 * t22 + t30; t22 * t3 - t4 * t36, t22 * t4 + t3 * t36, t5, t22 * t2 + (-qJ(2) - t28) * t25 + t32; t24 * t4, -t24 * t3, t21, t24 * t1 - t21 * t16 + t33; 0, 0, 0, 1;];
T_ges = t10;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
