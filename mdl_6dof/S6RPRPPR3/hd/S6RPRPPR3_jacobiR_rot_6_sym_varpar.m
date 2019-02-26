% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:24
% EndTime: 2019-02-26 20:40:24
% DurationCPUTime: 0.03s
% Computational Cost: add. (38->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t68 = sin(qJ(6));
t69 = sin(qJ(3));
t75 = t69 * t68;
t70 = cos(qJ(6));
t74 = t69 * t70;
t71 = cos(qJ(3));
t73 = t71 * t68;
t72 = t71 * t70;
t67 = qJ(1) + pkin(9);
t66 = cos(t67);
t65 = sin(t67);
t64 = -t65 * t68 + t66 * t74;
t63 = -t65 * t70 - t66 * t75;
t62 = -t65 * t74 - t66 * t68;
t61 = t65 * t75 - t66 * t70;
t1 = [t62, 0, t66 * t72, 0, 0, t63; t64, 0, t65 * t72, 0, 0, -t61; 0, 0, t74, 0, 0, t73; t61, 0, -t66 * t73, 0, 0, -t64; t63, 0, -t65 * t73, 0, 0, t62; 0, 0, -t75, 0, 0, t72; -t65 * t71, 0, -t66 * t69, 0, 0, 0; t66 * t71, 0, -t65 * t69, 0, 0, 0; 0, 0, t71, 0, 0, 0;];
JR_rot  = t1;
