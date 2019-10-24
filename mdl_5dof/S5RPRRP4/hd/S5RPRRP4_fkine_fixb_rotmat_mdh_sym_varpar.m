% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:44
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:44:22
% EndTime: 2019-10-24 10:44:23
% DurationCPUTime: 0.11s
% Computational Cost: add. (114->55), mult. (117->54), div. (0->0), fcn. (168->8), ass. (0->35)
t25 = -pkin(7) - pkin(6);
t23 = cos(qJ(3));
t9 = t23 * pkin(3) + pkin(2);
t18 = qJ(3) + qJ(4);
t10 = sin(t18);
t19 = sin(pkin(8));
t38 = t19 * t10;
t22 = sin(qJ(1));
t37 = t22 * t19;
t20 = cos(pkin(8));
t36 = t22 * t20;
t21 = sin(qJ(3));
t35 = t22 * t21;
t34 = t22 * t23;
t24 = cos(qJ(1));
t33 = t24 * t20;
t32 = t24 * t21;
t31 = t24 * t23;
t17 = pkin(5) + 0;
t30 = t24 * qJ(2) + 0;
t29 = t24 * pkin(1) + t22 * qJ(2) + 0;
t28 = pkin(2) * t20 + pkin(6) * t19;
t16 = -qJ(5) + t25;
t11 = cos(t18);
t5 = pkin(4) * t11 + t9;
t27 = -t16 * t19 + t20 * t5;
t26 = -t19 * t25 + t20 * t9;
t8 = t24 * t19;
t7 = t19 * t11;
t6 = t21 * pkin(3) + pkin(4) * t10;
t4 = t22 * t10 + t11 * t33;
t3 = -t10 * t33 + t22 * t11;
t2 = t24 * t10 - t11 * t36;
t1 = t10 * t36 + t24 * t11;
t12 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t17; -t22, -t24, 0, 0; t24, -t22, 0, 0; 0, 0, 0, 1; t19, t20, 0, t17; -t36, t37, t24, -t22 * pkin(1) + t30; t33, -t8, t22, t29; 0, 0, 0, 1; t19 * t23, -t19 * t21, -t20, t19 * pkin(2) - t20 * pkin(6) + t17; -t20 * t34 + t32, t20 * t35 + t31, -t37, (-pkin(1) - t28) * t22 + t30; t20 * t31 + t35, -t20 * t32 + t34, t8, t28 * t24 + t29; 0, 0, 0, 1; t7, -t38, -t20, t19 * t9 + t20 * t25 + t17; t2, t1, -t37, pkin(3) * t32 + (-pkin(1) - t26) * t22 + t30; t4, t3, t8, pkin(3) * t35 + t26 * t24 + t29; 0, 0, 0, 1; t7, -t38, -t20, t20 * t16 + t19 * t5 + t17; t2, t1, -t37, t24 * t6 + (-pkin(1) - t27) * t22 + t30; t4, t3, t8, t22 * t6 + t27 * t24 + t29; 0, 0, 0, 1;];
T_ges = t12;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
