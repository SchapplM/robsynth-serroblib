% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
%
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-06 15:46
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRRR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_jacobiR_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR3_jacobiR_rot_5_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-06 15:46:08
% EndTime: 2019-06-06 15:46:08
% DurationCPUTime: 0.03s
% Computational Cost: add. (55->12), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
t29 = qJ(2) + qJ(3) + qJ(4);
t27 = sin(t29);
t31 = cos(qJ(5));
t33 = t27 * t31;
t28 = cos(t29);
t30 = sin(qJ(5));
t32 = t28 * t30;
t26 = t28 * t31;
t25 = t27 * t30;
t1 = [0, -t33, -t33, -t33, -t32; 0, t26, t26, t26, -t25; 0, 0, 0, 0, t31; 0, t25, t25, t25, -t26; 0, -t32, -t32, -t32, -t33; 0, 0, 0, 0, -t30; 0, t28, t28, t28, 0; 0, t27, t27, t27, 0; 0, 0, 0, 0, 0;];
JR_rot  = t1;
