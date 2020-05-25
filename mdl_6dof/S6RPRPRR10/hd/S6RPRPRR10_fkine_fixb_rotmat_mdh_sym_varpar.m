% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 16:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRPRR10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:08:49
% EndTime: 2018-11-23 16:08:49
% DurationCPUTime: 0.13s
% Computational Cost: add. (140->63), mult. (129->66), div. (0->0), fcn. (184->10), ass. (0->35)
t20 = cos(pkin(10));
t5 = t20 * pkin(4) + pkin(3);
t19 = sin(pkin(10));
t23 = sin(qJ(1));
t39 = t23 * t19;
t22 = sin(qJ(3));
t38 = t23 * t22;
t24 = cos(qJ(3));
t37 = t23 * t24;
t25 = cos(qJ(1));
t36 = t25 * t19;
t35 = t25 * t22;
t21 = -pkin(8) - qJ(4);
t18 = pkin(6) + 0;
t17 = pkin(10) + qJ(5);
t34 = t23 * pkin(1) + 0;
t33 = pkin(2) + t18;
t12 = t23 * pkin(7);
t32 = t12 + t34;
t31 = t25 * pkin(1) + t23 * qJ(2) + 0;
t30 = t25 * pkin(7) + t31;
t8 = cos(t17);
t1 = pkin(5) * t8 + t5;
t16 = -pkin(9) + t21;
t29 = t1 * t22 + t16 * t24;
t28 = t21 * t24 + t22 * t5;
t27 = pkin(3) * t22 - qJ(4) * t24;
t26 = -t25 * qJ(2) + t34;
t9 = qJ(6) + t17;
t7 = sin(t17);
t6 = t25 * t24;
t4 = cos(t9);
t3 = sin(t9);
t2 = t19 * pkin(4) + pkin(5) * t7;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t25, -t23, 0, 0; t23, t25, 0, 0; 0, 0, 1, t18; 0, 0, 0, 1; 0, -t25, t23, t31; 0, -t23, -t25, t26; 1, 0, 0, t18; 0, 0, 0, 1; t38, t37, t25, t30; -t35, -t6, t23, t12 + t26; t24, -t22, 0, t33; 0, 0, 0, 1; t20 * t38 + t36, -t19 * t38 + t25 * t20, -t37, t27 * t23 + t30; -t20 * t35 + t39, t19 * t35 + t23 * t20, t6 (-qJ(2) - t27) * t25 + t32; t24 * t20, -t24 * t19, t22, t24 * pkin(3) + t22 * qJ(4) + t33; 0, 0, 0, 1; t25 * t7 + t8 * t38, t25 * t8 - t7 * t38, -t37, pkin(4) * t36 + t28 * t23 + t30; t23 * t7 - t8 * t35, t23 * t8 + t7 * t35, t6, pkin(4) * t39 + (-qJ(2) - t28) * t25 + t32; t24 * t8, -t24 * t7, t22, -t22 * t21 + t24 * t5 + t33; 0, 0, 0, 1; t25 * t3 + t4 * t38, t25 * t4 - t3 * t38, -t37, t25 * t2 + t29 * t23 + t30; t23 * t3 - t4 * t35, t23 * t4 + t3 * t35, t6, t23 * t2 + (-qJ(2) - t29) * t25 + t32; t24 * t4, -t24 * t3, t22, t24 * t1 - t22 * t16 + t33; 0, 0, 0, 1;];
T_ges = t10;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
