% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2018-11-23 15:39
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPPPRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:38:54
% EndTime: 2018-11-23 15:38:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (106->41), mult. (136->38), div. (0->0), fcn. (205->8), ass. (0->28)
t16 = sin(qJ(5));
t14 = cos(pkin(9));
t19 = cos(qJ(1));
t30 = sin(pkin(9));
t33 = sin(qJ(1));
t3 = t19 * t14 - t33 * t30;
t35 = t3 * t16;
t4 = t33 * t14 + t19 * t30;
t34 = t4 * t16;
t15 = sin(qJ(6));
t18 = cos(qJ(5));
t32 = t15 * t18;
t17 = cos(qJ(6));
t31 = t17 * t18;
t13 = pkin(6) + 0;
t29 = t33 * pkin(1) + 0;
t28 = pkin(2) + t13;
t27 = t19 * pkin(1) + t33 * qJ(2) + 0;
t26 = t19 * qJ(3) + t27;
t6 = qJ(4) + t28;
t25 = pkin(5) * t18 + pkin(8) * t16;
t24 = t33 * pkin(3) + t26;
t23 = -t19 * qJ(2) + t29;
t7 = t33 * qJ(3);
t22 = t7 + (-pkin(3) - qJ(2)) * t19 + t29;
t21 = t4 * pkin(4) - t3 * pkin(7) + t24;
t20 = -t3 * pkin(4) - t4 * pkin(7) + t22;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t33, 0, 0; t33, t19, 0, 0; 0, 0, 1, t13; 0, 0, 0, 1; 0, -t19, t33, t27; 0, -t33, -t19, t23; 1, 0, 0, t13; 0, 0, 0, 1; t33, 0, t19, t26; -t19, 0, t33, t23 + t7; 0, -1, 0, t28; 0, 0, 0, 1; t4, t3, 0, t24; -t3, t4, 0, t22; 0, 0, 1, t6; 0, 0, 0, 1; t4 * t18, -t34, -t3, t21; -t3 * t18, t35, -t4, t20; t16, t18, 0, t6; 0, 0, 0, 1; -t3 * t15 + t4 * t31, -t3 * t17 - t4 * t32, t34, t25 * t4 + t21; -t4 * t15 - t3 * t31, -t4 * t17 + t3 * t32, -t35, -t25 * t3 + t20; t16 * t17, -t16 * t15, -t18, t16 * pkin(5) - t18 * pkin(8) + t6; 0, 0, 0, 1;];
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
