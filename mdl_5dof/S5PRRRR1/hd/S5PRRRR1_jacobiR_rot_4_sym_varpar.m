% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRRR1_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobiR_rot_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobiR_rot_4_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:29:13
% EndTime: 2019-07-18 13:29:14
% DurationCPUTime: 0.02s
% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
t34 = qJ(3) + qJ(4);
t32 = sin(t34);
t35 = sin(qJ(2));
t40 = t35 * t32;
t33 = cos(t34);
t39 = t35 * t33;
t36 = cos(qJ(2));
t38 = t36 * t32;
t37 = t36 * t33;
t1 = [0, -t39, -t38, -t38, 0; 0, 0, -t33, -t33, 0; 0, t37, -t40, -t40, 0; 0, t40, -t37, -t37, 0; 0, 0, t32, t32, 0; 0, -t38, -t39, -t39, 0; 0, t36, 0, 0, 0; 0, 0, 0, 0, 0; 0, t35, 0, 0, 0;];
JR_rot  = t1;