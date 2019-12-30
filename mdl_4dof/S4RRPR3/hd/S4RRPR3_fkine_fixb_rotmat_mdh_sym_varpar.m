% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4RRPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 13:29:59
% EndTime: 2019-12-29 13:29:59
% DurationCPUTime: 0.13s
% Computational Cost: add. (65->19), mult. (18->12), div. (0->0), fcn. (38->8), ass. (0->18)
t11 = qJ(1) + qJ(2);
t21 = pkin(4) + 0;
t13 = sin(qJ(1));
t20 = t13 * pkin(1) + 0;
t15 = cos(qJ(1));
t19 = t15 * pkin(1) + 0;
t18 = pkin(5) + t21;
t7 = sin(t11);
t17 = pkin(2) * t7 + t20;
t8 = cos(t11);
t16 = pkin(2) * t8 + t19;
t14 = cos(qJ(4));
t12 = sin(qJ(4));
t6 = pkin(7) + t11;
t5 = qJ(3) + t18;
t2 = cos(t6);
t1 = sin(t6);
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t15, -t13, 0, 0; t13, t15, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t8, -t7, 0, t19; t7, t8, 0, t20; 0, 0, 1, t18; 0, 0, 0, 1; t2, -t1, 0, t16; t1, t2, 0, t17; 0, 0, 1, t5; 0, 0, 0, 1; t2 * t14, -t2 * t12, t1, t2 * pkin(3) + t1 * pkin(6) + t16; t1 * t14, -t1 * t12, -t2, t1 * pkin(3) - t2 * pkin(6) + t17; t12, t14, 0, t5; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
