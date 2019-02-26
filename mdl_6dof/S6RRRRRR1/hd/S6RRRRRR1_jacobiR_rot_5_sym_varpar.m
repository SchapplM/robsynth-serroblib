% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:47:03
% EndTime: 2019-02-26 22:47:03
% DurationCPUTime: 0.03s
% Computational Cost: add. (106->23), mult. (20->4), div. (0->0), fcn. (50->4), ass. (0->10)
t34 = qJ(2) + qJ(3) + qJ(4) + qJ(5);
t32 = sin(t34);
t35 = sin(qJ(1));
t40 = t35 * t32;
t33 = cos(t34);
t39 = t35 * t33;
t36 = cos(qJ(1));
t38 = t36 * t32;
t37 = t36 * t33;
t1 = [-t39, -t38, -t38, -t38, -t38, 0; t37, -t40, -t40, -t40, -t40, 0; 0, t33, t33, t33, t33, 0; t40, -t37, -t37, -t37, -t37, 0; -t38, -t39, -t39, -t39, -t39, 0; 0, -t32, -t32, -t32, -t32, 0; t36, 0, 0, 0, 0, 0; t35, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
