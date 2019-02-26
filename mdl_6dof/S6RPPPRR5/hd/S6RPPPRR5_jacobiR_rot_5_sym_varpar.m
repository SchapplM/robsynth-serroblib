% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPPRR5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:25:03
% EndTime: 2019-02-26 20:25:03
% DurationCPUTime: 0.02s
% Computational Cost: add. (17->8), mult. (28->8), div. (0->0), fcn. (50->6), ass. (0->13)
t32 = sin(qJ(1));
t23 = cos(pkin(9));
t26 = cos(qJ(1));
t27 = sin(pkin(9));
t20 = t26 * t23 - t32 * t27;
t24 = sin(qJ(5));
t31 = t20 * t24;
t25 = cos(qJ(5));
t30 = t20 * t25;
t21 = t32 * t23 + t26 * t27;
t29 = t21 * t24;
t28 = t21 * t25;
t1 = [t30, 0, 0, 0, -t29, 0; t28, 0, 0, 0, t31, 0; 0, 0, 0, 0, t25, 0; -t31, 0, 0, 0, -t28, 0; -t29, 0, 0, 0, t30, 0; 0, 0, 0, 0, -t24, 0; t21, 0, 0, 0, 0, 0; -t20, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
