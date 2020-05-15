% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2018-11-23 16:14
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RPRRPP7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:14:40
% EndTime: 2018-11-23 16:14:40
% DurationCPUTime: 0.14s
% Computational Cost: add. (108->51), mult. (165->40), div. (0->0), fcn. (230->6), ass. (0->35)
t24 = sin(qJ(4));
t25 = sin(qJ(3));
t29 = cos(qJ(1));
t45 = t29 * t25;
t26 = sin(qJ(1));
t27 = cos(qJ(4));
t46 = t26 * t27;
t5 = t24 * t45 + t46;
t44 = t29 * t27;
t6 = t26 * t24 - t25 * t44;
t48 = t6 * pkin(4) - t5 * qJ(5);
t47 = t26 * t25;
t28 = cos(qJ(3));
t11 = t26 * t28;
t12 = t28 * t24;
t13 = t28 * t27;
t15 = t29 * t28;
t23 = pkin(6) + 0;
t42 = pkin(8) * t11;
t41 = t26 * pkin(1) + 0;
t40 = pkin(2) + t23;
t39 = -pkin(3) * t25 - qJ(2);
t38 = t29 * pkin(1) + t26 * qJ(2) + 0;
t18 = t26 * pkin(7);
t37 = pkin(8) * t15 + t18 + t41;
t36 = t29 * pkin(7) + t38;
t35 = t28 * pkin(3) + t25 * pkin(8) + t40;
t34 = pkin(3) * t47 + t36;
t33 = -t29 * qJ(2) + t41;
t32 = pkin(4) * t13 + qJ(5) * t12 + t35;
t3 = t24 * t47 - t44;
t4 = t29 * t24 + t25 * t46;
t31 = t4 * pkin(4) + t3 * qJ(5) + t34;
t30 = t39 * t29 + t37;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t26, 0, 0; t26, t29, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; 0, -t29, t26, t38; 0, -t26, -t29, t33; 1, 0, 0, t23; 0, 0, 0, 1; t47, t11, t29, t36; -t45, -t15, t26, t18 + t33; t28, -t25, 0, t40; 0, 0, 0, 1; t4, -t3, -t11, t34 - t42; t6, t5, t15, t30; t13, -t12, t25, t35; 0, 0, 0, 1; t4, -t11, t3, t31 - t42; t6, t15, -t5, t30 + t48; t13, t25, t12, t32; 0, 0, 0, 1; t4, t3, t11, t4 * pkin(5) + (-pkin(8) + qJ(6)) * t11 + t31; t6, -t5, -t15, t6 * pkin(5) + (-qJ(6) * t28 + t39) * t29 + t37 + t48; t13, t12, -t25, pkin(5) * t13 - t25 * qJ(6) + t32; 0, 0, 0, 1;];
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
