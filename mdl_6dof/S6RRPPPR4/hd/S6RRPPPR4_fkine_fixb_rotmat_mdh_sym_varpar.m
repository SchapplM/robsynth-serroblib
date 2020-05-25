% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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

function [T_c_mdh, Tc_stack] = S6RRPPPR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:43:46
% EndTime: 2018-11-23 16:43:46
% DurationCPUTime: 0.16s
% Computational Cost: add. (130->61), mult. (218->59), div. (0->0), fcn. (297->8), ass. (0->35)
t28 = sin(qJ(1));
t30 = cos(qJ(2));
t15 = t28 * t30;
t27 = sin(qJ(2));
t44 = qJ(3) * t27;
t49 = pkin(2) * t15 + t28 * t44;
t48 = t28 * t27;
t24 = sin(pkin(9));
t47 = t30 * t24;
t25 = cos(pkin(9));
t46 = t30 * t25;
t31 = cos(qJ(1));
t45 = t31 * t27;
t16 = t31 * t30;
t43 = qJ(4) * t30;
t23 = pkin(6) + 0;
t42 = t28 * pkin(1) + 0;
t41 = t27 * pkin(2) + t23;
t40 = t31 * pkin(1) + t28 * pkin(7) + 0;
t39 = -t31 * pkin(7) + t42;
t17 = t27 * qJ(4);
t38 = qJ(5) * t46 + t17 + t41;
t37 = pkin(2) * t16 + t31 * t44 + t40;
t36 = -t30 * qJ(3) + t41;
t35 = t28 * pkin(3) + t31 * t43 + t37;
t34 = t28 * t43 + (-pkin(3) - pkin(7)) * t31 + t42 + t49;
t3 = t28 * t24 - t25 * t45;
t4 = t24 * t45 + t28 * t25;
t33 = t4 * pkin(4) + t3 * qJ(5) + t35;
t5 = t31 * t24 + t25 * t48;
t6 = t24 * t48 - t31 * t25;
t32 = t6 * pkin(4) - t5 * qJ(5) + t34;
t29 = cos(qJ(6));
t26 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t31, -t28, 0, 0; t28, t31, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t16, -t45, t28, t40; t15, -t48, -t31, t39; t27, t30, 0, t23; 0, 0, 0, 1; t28, -t16, t45, t37; -t31, -t15, t48, t39 + t49; 0, -t27, -t30, t36; 0, 0, 0, 1; t4, -t3, t16, t35; t6, t5, t15, t34; -t47, -t46, t27, t17 + t36; 0, 0, 0, 1; t4, t16, t3, t33; t6, t15, -t5, t32; -t47, t27, t46 (-pkin(4) * t24 - qJ(3)) * t30 + t38; 0, 0, 0, 1; t3 * t26 + t4 * t29, -t4 * t26 + t3 * t29, -t16, t4 * pkin(5) - pkin(8) * t16 + t33; -t5 * t26 + t6 * t29, -t6 * t26 - t5 * t29, -t15, t6 * pkin(5) - pkin(8) * t15 + t32; (-t24 * t29 + t25 * t26) * t30 (t24 * t26 + t25 * t29) * t30, -t27, -t27 * pkin(8) + (-qJ(3) + (-pkin(4) - pkin(5)) * t24) * t30 + t38; 0, 0, 0, 1;];
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
