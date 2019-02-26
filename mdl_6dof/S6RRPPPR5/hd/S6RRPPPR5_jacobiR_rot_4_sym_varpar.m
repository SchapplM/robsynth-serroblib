% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR5_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobiR_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:13
% EndTime: 2019-02-26 21:24:13
% DurationCPUTime: 0.03s
% Computational Cost: add. (8->8), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
t59 = sin(qJ(2));
t60 = sin(qJ(1));
t67 = t60 * t59;
t57 = sin(pkin(9));
t61 = cos(qJ(2));
t66 = t61 * t57;
t58 = cos(pkin(9));
t65 = t61 * t58;
t62 = cos(qJ(1));
t64 = t62 * t59;
t63 = t62 * t61;
t1 = [-t67, t63, 0, 0, 0, 0; t64, t60 * t61, 0, 0, 0, 0; 0, t59, 0, 0, 0, 0; -t62 * t57 + t60 * t65, t58 * t64, 0, 0, 0, 0; -t60 * t57 - t58 * t63, t58 * t67, 0, 0, 0, 0; 0, -t65, 0, 0, 0, 0; -t62 * t58 - t60 * t66, -t57 * t64, 0, 0, 0, 0; t57 * t63 - t60 * t58, -t57 * t67, 0, 0, 0, 0; 0, t66, 0, 0, 0, 0;];
JR_rot  = t1;
