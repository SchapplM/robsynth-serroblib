% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2018-11-23 16:59
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6RRPRPP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:59:01
% EndTime: 2018-11-23 16:59:02
% DurationCPUTime: 0.14s
% Computational Cost: add. (121->57), mult. (192->45), div. (0->0), fcn. (261->6), ass. (0->35)
t27 = sin(qJ(1));
t29 = cos(qJ(2));
t15 = t27 * t29;
t26 = sin(qJ(2));
t43 = qJ(3) * t26;
t49 = pkin(2) * t15 + t27 * t43;
t48 = t27 * t26;
t28 = cos(qJ(4));
t47 = t27 * t28;
t25 = sin(qJ(4));
t46 = t29 * t25;
t16 = t29 * t28;
t30 = cos(qJ(1));
t45 = t30 * t26;
t44 = t30 * t28;
t17 = t30 * t29;
t42 = qJ(6) * t29;
t24 = pkin(6) + 0;
t41 = t27 * pkin(1) + 0;
t40 = t26 * pkin(2) + t24;
t39 = t30 * pkin(1) + t27 * pkin(7) + 0;
t38 = -t30 * pkin(7) + t41;
t18 = t26 * pkin(8);
t37 = qJ(5) * t16 + t18 + t40;
t36 = pkin(2) * t17 + t30 * t43 + t39;
t35 = -t29 * qJ(3) + t40;
t34 = t27 * pkin(3) + pkin(8) * t17 + t36;
t33 = pkin(8) * t15 + (-pkin(3) - pkin(7)) * t30 + t41 + t49;
t3 = t27 * t25 - t26 * t44;
t4 = t25 * t45 + t47;
t32 = t4 * pkin(4) + t3 * qJ(5) + t34;
t5 = t30 * t25 + t26 * t47;
t6 = t25 * t48 - t44;
t31 = t6 * pkin(4) - t5 * qJ(5) + t33;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t27, 0, 0; t27, t30, 0, 0; 0, 0, 1, t24; 0, 0, 0, 1; t17, -t45, t27, t39; t15, -t48, -t30, t38; t26, t29, 0, t24; 0, 0, 0, 1; t27, -t17, t45, t36; -t30, -t15, t48, t38 + t49; 0, -t26, -t29, t35; 0, 0, 0, 1; t4, -t3, t17, t34; t6, t5, t15, t33; -t46, -t16, t26, t18 + t35; 0, 0, 0, 1; t4, t17, t3, t32; t6, t15, -t5, t31; -t46, t26, t16 (-pkin(4) * t25 - qJ(3)) * t29 + t37; 0, 0, 0, 1; t4, t3, -t17, t4 * pkin(5) - t30 * t42 + t32; t6, -t5, -t15, t6 * pkin(5) - t27 * t42 + t31; -t46, t16, -t26, -t26 * qJ(6) + (-qJ(3) + (-pkin(4) - pkin(5)) * t25) * t29 + t37; 0, 0, 0, 1;];
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
