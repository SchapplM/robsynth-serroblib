% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2018-11-23 17:43
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRRPRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:42:44
% EndTime: 2018-11-23 17:42:45
% DurationCPUTime: 0.11s
% Computational Cost: add. (167->53), mult. (140->47), div. (0->0), fcn. (198->8), ass. (0->38)
t25 = qJ(2) + qJ(3);
t21 = cos(t25);
t26 = sin(qJ(5));
t50 = t21 * t26;
t29 = cos(qJ(5));
t49 = t21 * t29;
t20 = sin(t25);
t28 = sin(qJ(1));
t48 = t28 * t20;
t14 = t28 * t21;
t47 = t28 * t26;
t46 = t28 * t29;
t31 = cos(qJ(1));
t45 = t31 * t20;
t15 = t31 * t21;
t44 = t31 * t26;
t43 = t31 * t29;
t42 = qJ(4) * t20;
t24 = pkin(6) + 0;
t27 = sin(qJ(2));
t41 = t27 * pkin(2) + t24;
t30 = cos(qJ(2));
t18 = t30 * pkin(2) + pkin(1);
t32 = -pkin(8) - pkin(7);
t40 = t28 * t18 + t31 * t32 + 0;
t39 = t20 * pkin(3) + t41;
t38 = t31 * t18 - t28 * t32 + 0;
t37 = pkin(3) * t14 + t28 * t42 + t40;
t36 = pkin(3) * t15 + t31 * t42 + t38;
t35 = -t21 * qJ(4) + t39;
t34 = -t31 * pkin(4) + pkin(9) * t14 + t37;
t33 = t28 * pkin(4) + pkin(9) * t15 + t36;
t16 = t20 * pkin(9);
t4 = t20 * t47 - t43;
t3 = t20 * t46 + t44;
t2 = t20 * t44 + t46;
t1 = -t20 * t43 + t47;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t31, -t28, 0, 0; t28, t31, 0, 0; 0, 0, 1, t24; 0, 0, 0, 1; t31 * t30, -t31 * t27, t28, t31 * pkin(1) + t28 * pkin(7) + 0; t28 * t30, -t28 * t27, -t31, t28 * pkin(1) - t31 * pkin(7) + 0; t27, t30, 0, t24; 0, 0, 0, 1; t15, -t45, t28, t38; t14, -t48, -t31, t40; t20, t21, 0, t41; 0, 0, 0, 1; t28, -t15, t45, t36; -t31, -t14, t48, t37; 0, -t20, -t21, t35; 0, 0, 0, 1; t2, -t1, t15, t33; t4, t3, t14, t34; -t50, -t49, t20, t16 + t35; 0, 0, 0, 1; t2, t15, t1, t2 * pkin(5) + t1 * qJ(6) + t33; t4, t14, -t3, t4 * pkin(5) - t3 * qJ(6) + t34; -t50, t20, t49, t16 + (-pkin(5) * t26 + qJ(6) * t29 - qJ(4)) * t21 + t39; 0, 0, 0, 1;];
T_ges = t5;
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
