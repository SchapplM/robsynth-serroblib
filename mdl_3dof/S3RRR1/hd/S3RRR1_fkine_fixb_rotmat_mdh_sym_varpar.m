% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% Tc_mdh [4x4x(3+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   4:  mdh base (link 0) -> mdh frame (4-1), link (4-1)
%   ...
%   3+1:  mdh base (link 0) -> mdh frame (3)
% T_c_stack [(3+1)*3 x 4]
%   stacked matrices from Tc_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2026-07-26 17:44
% Revision: 05cfe754a3b58d1958e6f60508538bf5b2996c81 (2026-06-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Tc_mdh, Tc_stack] = S3RRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2026-07-26 17:44:32
% EndTime: 2026-07-26 17:44:32
% DurationCPUTime: 0.11s
% Computational Cost: add. (29->12), mult. (6->4), div. (0->0), fcn. (18->6), ass. (0->13)
t8 = qJ(1) + qJ(2);
t14 = pkin(3) + 0;
t9 = sin(qJ(1));
t13 = t9 * pkin(1) + 0;
t10 = cos(qJ(1));
t12 = t10 * pkin(1) + 0;
t11 = pkin(4) + t14;
t5 = qJ(3) + t8;
t4 = cos(t8);
t3 = sin(t8);
t2 = cos(t5);
t1 = sin(t5);
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; t10, -t9, 0, 0; t9, t10, 0, 0; 0, 0, 1, t14; t4, -t3, 0, t12; t3, t4, 0, t13; 0, 0, 1, t11; t2, -t1, 0, pkin(2) * t4 + t12; t1, t2, 0, pkin(2) * t3 + t13; 0, 0, 1, pkin(5) + t11;];
Tc_stack = t6;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), Tc_mdh = NaN(4,4,3+1);               % numerisch
else,                         Tc_mdh = sym('xx', [4,4,3+1]); end % symbolisch
for i = 1:3+1
  Tc_mdh(:,:,i) = [Tc_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
