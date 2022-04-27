% Calculate homogenous joint transformation matrices for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-02-03 18:11
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5RRPRP8_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-02-03 18:11:23
% EndTime: 2022-02-03 18:11:23
% DurationCPUTime: 0.02s
% Computational Cost: add. (6->6), mult. (0->0), div. (0->0), fcn. (12->6), ass. (0->7)
t39 = cos(qJ(1));
t38 = cos(qJ(2));
t37 = cos(qJ(4));
t36 = sin(qJ(1));
t35 = sin(qJ(2));
t34 = sin(qJ(4));
t1 = [t39, -t36, 0, 0; t36, t39, 0, 0; 0, 0, 1, pkin(5); t38, -t35, 0, pkin(1); 0, 0, -1, -pkin(6); t35, t38, 0, 0; 1, 0, 0, pkin(2); 0, 0, -1, -qJ(3); 0, 1, 0, 0; t37, -t34, 0, pkin(3); 0, 0, -1, -pkin(7); t34, t37, 0, 0; 1, 0, 0, pkin(4); 0, 1, 0, 0; 0, 0, 1, qJ(5);];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
