% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRP5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_jacobiR_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:32:33
% EndTime: 2019-02-26 20:32:33
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t65 = sin(qJ(4));
t66 = sin(qJ(1));
t76 = t66 * t65;
t67 = cos(qJ(5));
t75 = t66 * t67;
t64 = sin(qJ(5));
t68 = cos(qJ(4));
t74 = t68 * t64;
t73 = t68 * t67;
t69 = cos(qJ(1));
t72 = t69 * t65;
t71 = t69 * t67;
t70 = t69 * t68;
t63 = -t66 * t64 + t65 * t71;
t62 = -t64 * t72 - t75;
t61 = -t69 * t64 - t65 * t75;
t60 = t64 * t76 - t71;
t1 = [t61, 0, 0, t67 * t70, t62, 0; t63, 0, 0, t66 * t73, -t60, 0; 0, 0, 0, -t65 * t67, -t74, 0; t60, 0, 0, -t64 * t70, -t63, 0; t62, 0, 0, -t66 * t74, t61, 0; 0, 0, 0, t65 * t64, -t73, 0; t66 * t68, 0, 0, t72, 0, 0; -t70, 0, 0, t76, 0, 0; 0, 0, 0, t68, 0, 0;];
JR_rot  = t1;
