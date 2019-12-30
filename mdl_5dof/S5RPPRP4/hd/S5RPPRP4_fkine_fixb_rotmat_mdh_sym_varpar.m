% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-29 16:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 16:00:52
% EndTime: 2019-12-29 16:00:52
% DurationCPUTime: 0.18s
% Computational Cost: add. (77->36), mult. (89->24), div. (0->0), fcn. (141->6), ass. (0->21)
t17 = sin(qJ(4));
t19 = cos(qJ(1));
t24 = sin(pkin(7));
t25 = cos(pkin(7));
t26 = sin(qJ(1));
t3 = -t19 * t25 - t26 * t24;
t28 = t3 * t17;
t4 = t19 * t24 - t26 * t25;
t27 = t4 * t17;
t15 = pkin(5) + 0;
t9 = -qJ(3) + t15;
t23 = t19 * pkin(1) + t26 * qJ(2) + 0;
t22 = t19 * pkin(2) + t23;
t21 = t26 * pkin(1) - t19 * qJ(2) + 0;
t20 = t26 * pkin(2) + t21;
t18 = cos(qJ(4));
t16 = -qJ(5) - pkin(6);
t8 = t18 * pkin(4) + pkin(3);
t2 = t3 * t18;
t1 = t4 * t18;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t26, 0, 0; t26, t19, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t19, 0, t26, t23; t26, 0, -t19, t21; 0, 1, 0, t15; 0, 0, 0, 1; -t3, -t4, 0, t22; -t4, t3, 0, t20; 0, 0, -1, t9; 0, 0, 0, 1; -t2, t28, t4, -t3 * pkin(3) + t4 * pkin(6) + t22; -t1, t27, -t3, -t4 * pkin(3) - t3 * pkin(6) + t20; -t17, -t18, 0, t9; 0, 0, 0, 1; -t2, t28, t4, -t4 * t16 - t3 * t8 + t22; -t1, t27, -t3, t3 * t16 - t4 * t8 + t20; -t17, -t18, 0, -t17 * pkin(4) + t9; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
