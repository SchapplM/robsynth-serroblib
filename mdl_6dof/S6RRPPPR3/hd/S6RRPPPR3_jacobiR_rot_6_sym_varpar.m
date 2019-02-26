% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPPR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:08
% EndTime: 2019-02-26 21:23:08
% DurationCPUTime: 0.05s
% Computational Cost: add. (38->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t71 = sin(qJ(2));
t72 = sin(qJ(1));
t79 = t72 * t71;
t70 = pkin(9) + qJ(6);
t68 = sin(t70);
t73 = cos(qJ(2));
t78 = t73 * t68;
t69 = cos(t70);
t77 = t73 * t69;
t74 = cos(qJ(1));
t76 = t74 * t71;
t75 = t74 * t73;
t67 = -t72 * t68 + t69 * t76;
t66 = -t68 * t76 - t72 * t69;
t65 = -t74 * t68 - t69 * t79;
t64 = t68 * t79 - t74 * t69;
t1 = [t65, t69 * t75, 0, 0, 0, t66; t67, t72 * t77, 0, 0, 0, -t64; 0, t71 * t69, 0, 0, 0, t78; t64, -t68 * t75, 0, 0, 0, -t67; t66, -t72 * t78, 0, 0, 0, t65; 0, -t71 * t68, 0, 0, 0, t77; -t72 * t73, -t76, 0, 0, 0, 0; t75, -t79, 0, 0, 0, 0; 0, t73, 0, 0, 0, 0;];
JR_rot  = t1;
