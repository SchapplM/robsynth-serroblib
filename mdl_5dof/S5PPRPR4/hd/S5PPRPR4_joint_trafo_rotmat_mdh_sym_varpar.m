% Calculate homogenous joint transformation matrices for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5PPRPR4_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:11
% EndTime: 2019-12-31 17:32:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (7->7), mult. (0->0), div. (0->0), fcn. (16->8), ass. (0->9)
t35 = cos(qJ(3));
t34 = cos(qJ(5));
t33 = sin(qJ(3));
t32 = sin(qJ(5));
t31 = cos(pkin(7));
t30 = cos(pkin(8));
t29 = sin(pkin(7));
t28 = sin(pkin(8));
t1 = [t31, -t29, 0, 0; t29, t31, 0, 0; 0, 0, 1, qJ(1); 0, 0, 0, 1; 1, 0, 0, pkin(1); 0, 0, -1, -qJ(2); 0, 1, 0, 0; 0, 0, 0, 1; t35, -t33, 0, pkin(2); 0, 0, -1, -pkin(5); t33, t35, 0, 0; 0, 0, 0, 1; t30, -t28, 0, pkin(3); 0, 0, -1, -qJ(4); t28, t30, 0, 0; 0, 0, 0, 1; t34, -t32, 0, pkin(4); t32, t34, 0, 0; 0, 0, 1, pkin(6); 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
