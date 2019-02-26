% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPPRR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:24:17
% EndTime: 2019-02-26 20:24:17
% DurationCPUTime: 0.02s
% Computational Cost: add. (12->5), mult. (28->8), div. (0->0), fcn. (50->6), ass. (0->13)
t34 = cos(qJ(1));
t33 = sin(qJ(1));
t28 = sin(pkin(9));
t29 = cos(pkin(9));
t22 = -t33 * t28 - t34 * t29;
t26 = sin(qJ(5));
t32 = t22 * t26;
t27 = cos(qJ(5));
t31 = t22 * t27;
t23 = t34 * t28 - t33 * t29;
t30 = t23 * t26;
t21 = t23 * t27;
t1 = [t32, 0, 0, 0, t21, 0; t30, 0, 0, 0, -t31, 0; 0, 0, 0, 0, t26, 0; t31, 0, 0, 0, -t30, 0; t21, 0, 0, 0, t32, 0; 0, 0, 0, 0, t27, 0; t23, 0, 0, 0, 0, 0; -t22, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
