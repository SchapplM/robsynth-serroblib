% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2018-11-23 15:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRPPR8_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:55:56
% EndTime: 2018-11-23 15:55:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (92->53), mult. (110->42), div. (0->0), fcn. (155->6), ass. (0->27)
t20 = cos(qJ(1));
t16 = sin(qJ(3));
t17 = sin(qJ(1));
t4 = t17 * t16;
t37 = pkin(4) * t4 - t20 * qJ(5);
t36 = -pkin(3) - pkin(4);
t19 = cos(qJ(3));
t35 = t17 * t19;
t34 = t20 * t16;
t5 = t20 * t19;
t33 = qJ(4) * t19;
t14 = pkin(6) + 0;
t31 = t17 * pkin(1) + 0;
t30 = pkin(2) + t14;
t29 = t20 * pkin(1) + t17 * qJ(2) + 0;
t8 = t17 * pkin(7);
t28 = t20 * t33 + t31 + t8;
t27 = t20 * pkin(7) + t29;
t26 = t19 * pkin(3) + t16 * qJ(4) + t30;
t25 = pkin(3) * t4 + t27;
t24 = -t20 * qJ(2) + t31;
t23 = t19 * pkin(4) + t26;
t22 = -t17 * qJ(5) + t28;
t21 = -t17 * t33 + t25;
t18 = cos(qJ(6));
t15 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t20, -t17, 0, 0; t17, t20, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; 0, -t20, t17, t29; 0, -t17, -t20, t24; 1, 0, 0, t14; 0, 0, 0, 1; t4, t35, t20, t27; -t34, -t5, t17, t24 + t8; t19, -t16, 0, t30; 0, 0, 0, 1; t4, t20, -t35, t21; -t34, t17, t5 (-pkin(3) * t16 - qJ(2)) * t20 + t28; t19, 0, t16, t26; 0, 0, 0, 1; -t35, -t4, -t20, t21 + t37; t5, t34, -t17 (t36 * t16 - qJ(2)) * t20 + t22; t16, -t19, 0, t23; 0, 0, 0, 1; -t20 * t15 - t18 * t35, t15 * t35 - t20 * t18, t4 (pkin(8) * t16 + (-pkin(5) - qJ(4)) * t19) * t17 + t25 + t37; -t17 * t15 + t18 * t5, -t15 * t5 - t17 * t18, -t34 (pkin(5) * t19 - qJ(2) + (-pkin(8) + t36) * t16) * t20 + t22; t16 * t18, -t16 * t15, t19, t16 * pkin(5) + t19 * pkin(8) + t23; 0, 0, 0, 1;];
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
