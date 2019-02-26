% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPP3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:56
% EndTime: 2019-02-26 21:35:56
% DurationCPUTime: 0.03s
% Computational Cost: add. (36->11), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t80 = sin(qJ(2));
t81 = sin(qJ(1));
t88 = t81 * t80;
t79 = pkin(9) + qJ(4);
t77 = sin(t79);
t82 = cos(qJ(2));
t87 = t82 * t77;
t78 = cos(t79);
t86 = t82 * t78;
t83 = cos(qJ(1));
t85 = t83 * t80;
t84 = t83 * t82;
t76 = t81 * t77 + t78 * t84;
t75 = t77 * t84 - t81 * t78;
t74 = -t83 * t77 + t81 * t86;
t73 = t83 * t78 + t81 * t87;
t1 = [-t88, t84, 0, 0, 0, 0; t85, t81 * t82, 0, 0, 0, 0; 0, t80, 0, 0, 0, 0; t74, t78 * t85, 0, t75, 0, 0; -t76, t78 * t88, 0, t73, 0, 0; 0, -t86, 0, t80 * t77, 0, 0; -t73, -t77 * t85, 0, t76, 0, 0; t75, -t77 * t88, 0, t74, 0, 0; 0, t87, 0, t80 * t78, 0, 0;];
JR_rot  = t1;
