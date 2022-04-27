% Calculate homogenous joint transformation matrices for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-02-03 15:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5RPRRR14_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-02-03 15:18:36
% EndTime: 2022-02-03 15:18:36
% DurationCPUTime: 0.09s
% Computational Cost: add. (11->11), mult. (12->12), div. (0->0), fcn. (40->14), ass. (0->15)
t150 = cos(qJ(1));
t149 = cos(qJ(3));
t148 = cos(qJ(4));
t147 = cos(qJ(5));
t146 = sin(qJ(1));
t145 = sin(qJ(3));
t144 = sin(qJ(4));
t143 = sin(qJ(5));
t142 = cos(pkin(5));
t141 = cos(pkin(6));
t140 = cos(pkin(11));
t139 = sin(pkin(5));
t138 = sin(pkin(6));
t137 = sin(pkin(11));
t1 = [t150, -t146, 0, 0; t146, t150, 0, 0; 0, 0, 1, pkin(7); t140, -t137, 0, pkin(1); t142 * t137, t142 * t140, -t139, -t139 * qJ(2); t139 * t137, t139 * t140, t142, t142 * qJ(2); t149, -t145, 0, pkin(2); t141 * t145, t141 * t149, -t138, -t138 * pkin(8); t138 * t145, t138 * t149, t141, t141 * pkin(8); t148, -t144, 0, pkin(3); 0, 0, -1, -pkin(9); t144, t148, 0, 0; t147, -t143, 0, pkin(4); 0, 0, -1, -pkin(10); t143, t147, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
