% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRPRR1_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_jacobiR_rot_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_jacobiR_rot_4_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:52
% EndTime: 2019-07-18 17:22:52
% DurationCPUTime: 0.02s
% Computational Cost: add. (28->13), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->10)
t22 = qJ(2) + qJ(4);
t20 = sin(t22);
t23 = sin(qJ(1));
t28 = t23 * t20;
t21 = cos(t22);
t27 = t23 * t21;
t24 = cos(qJ(1));
t26 = t24 * t20;
t25 = t24 * t21;
t1 = [-t27, -t26, 0, -t26, 0; t25, -t28, 0, -t28, 0; 0, t21, 0, t21, 0; t28, -t25, 0, -t25, 0; -t26, -t27, 0, -t27, 0; 0, -t20, 0, -t20, 0; t24, 0, 0, 0, 0; t23, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
JR_rot  = t1;
