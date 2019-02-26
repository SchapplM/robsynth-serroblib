% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR6_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:25
% EndTime: 2019-02-26 21:31:25
% DurationCPUTime: 0.03s
% Computational Cost: add. (12->7), mult. (36->8), div. (0->0), fcn. (58->6), ass. (0->13)
t22 = sin(pkin(10));
t23 = cos(pkin(10));
t24 = sin(qJ(2));
t26 = cos(qJ(2));
t29 = t26 * t22 - t24 * t23;
t28 = t24 * t22 + t26 * t23;
t27 = cos(qJ(1));
t25 = sin(qJ(1));
t21 = t28 * t27;
t20 = t29 * t27;
t19 = t28 * t25;
t18 = t29 * t25;
t1 = [-t19, t20, 0, 0, 0, 0; t21, t18, 0, 0, 0, 0; 0, t28, 0, 0, 0, 0; t18, t21, 0, 0, 0, 0; -t20, t19, 0, 0, 0, 0; 0, -t29, 0, 0, 0, 0; -t27, 0, 0, 0, 0, 0; -t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
