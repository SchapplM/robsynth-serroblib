% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRRPP3
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
% Datum: 2018-11-23 18:06
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRRPP3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:05:33
% EndTime: 2018-11-23 18:05:33
% DurationCPUTime: 0.13s
% Computational Cost: add. (185->53), mult. (173->46), div. (0->0), fcn. (242->8), ass. (0->37)
t28 = qJ(2) + qJ(3);
t24 = sin(t28);
t29 = sin(qJ(4));
t15 = t24 * t29;
t32 = cos(qJ(4));
t16 = t24 * t32;
t51 = pkin(4) * t16 + qJ(5) * t15;
t31 = sin(qJ(1));
t17 = t31 * t24;
t25 = cos(t28);
t50 = t31 * t25;
t49 = t31 * t29;
t48 = t31 * t32;
t34 = cos(qJ(1));
t18 = t34 * t24;
t47 = t34 * t25;
t46 = t34 * t29;
t45 = t34 * t32;
t27 = pkin(6) + 0;
t30 = sin(qJ(2));
t44 = t30 * pkin(2) + t27;
t33 = cos(qJ(2));
t22 = t33 * pkin(2) + pkin(1);
t35 = -pkin(8) - pkin(7);
t43 = t31 * t22 + t34 * t35 + 0;
t42 = t24 * pkin(3) + t44;
t41 = t34 * t22 - t31 * t35 + 0;
t40 = pkin(3) * t50 + pkin(9) * t17 + t43;
t39 = -t25 * pkin(9) + t42;
t38 = pkin(3) * t47 + pkin(9) * t18 + t41;
t3 = t25 * t49 + t45;
t4 = t25 * t48 - t46;
t37 = t4 * pkin(4) + t3 * qJ(5) + t40;
t5 = t25 * t46 - t48;
t6 = t25 * t45 + t49;
t36 = t6 * pkin(4) + t5 * qJ(5) + t38;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t31, 0, 0; t31, t34, 0, 0; 0, 0, 1, t27; 0, 0, 0, 1; t34 * t33, -t34 * t30, t31, t34 * pkin(1) + t31 * pkin(7) + 0; t31 * t33, -t31 * t30, -t34, t31 * pkin(1) - t34 * pkin(7) + 0; t30, t33, 0, t27; 0, 0, 0, 1; t47, -t18, t31, t41; t50, -t17, -t34, t43; t24, t25, 0, t44; 0, 0, 0, 1; t6, -t5, t18, t38; t4, -t3, t17, t40; t16, -t15, -t25, t39; 0, 0, 0, 1; t18, -t6, t5, t36; t17, -t4, t3, t37; -t25, -t16, t15, t39 + t51; 0, 0, 0, 1; t18, t5, t6, pkin(5) * t18 + t6 * qJ(6) + t36; t17, t3, t4, pkin(5) * t17 + t4 * qJ(6) + t37; -t25, t15, t16, qJ(6) * t16 + (-pkin(5) - pkin(9)) * t25 + t42 + t51; 0, 0, 0, 1;];
T_ges = t1;
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
