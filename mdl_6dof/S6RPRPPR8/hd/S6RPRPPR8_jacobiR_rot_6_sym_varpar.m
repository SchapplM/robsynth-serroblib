% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPPR8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobiR_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:07
% EndTime: 2019-02-26 20:43:07
% DurationCPUTime: 0.04s
% Computational Cost: add. (15->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t64 = sin(qJ(3));
t65 = sin(qJ(1));
t75 = t65 * t64;
t66 = cos(qJ(6));
t74 = t65 * t66;
t63 = sin(qJ(6));
t67 = cos(qJ(3));
t73 = t67 * t63;
t72 = t67 * t66;
t68 = cos(qJ(1));
t71 = t68 * t64;
t70 = t68 * t66;
t69 = t68 * t67;
t62 = t65 * t63 - t66 * t69;
t61 = t63 * t69 + t74;
t60 = t68 * t63 + t65 * t72;
t59 = t65 * t73 - t70;
t1 = [t62, 0, t64 * t74, 0, 0, t59; -t60, 0, -t64 * t70, 0, 0, -t61; 0, 0, t72, 0, 0, -t64 * t63; t61, 0, -t63 * t75, 0, 0, t60; t59, 0, t63 * t71, 0, 0, t62; 0, 0, -t73, 0, 0, -t64 * t66; t71, 0, t65 * t67, 0, 0, 0; t75, 0, -t69, 0, 0, 0; 0, 0, -t64, 0, 0, 0;];
JR_rot  = t1;
