% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JR_rot = S6RRPRPP3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:56
% EndTime: 2019-02-26 21:35:56
% DurationCPUTime: 0.03s
% Computational Cost: add. (38->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
t79 = sin(qJ(2));
t80 = sin(qJ(1));
t87 = t80 * t79;
t78 = pkin(9) + qJ(4);
t76 = sin(t78);
t81 = cos(qJ(2));
t86 = t81 * t76;
t77 = cos(t78);
t85 = t81 * t77;
t82 = cos(qJ(1));
t84 = t82 * t79;
t83 = t82 * t81;
t75 = t80 * t76 + t77 * t83;
t74 = t76 * t83 - t80 * t77;
t73 = -t82 * t76 + t80 * t85;
t72 = -t82 * t77 - t80 * t86;
t1 = [-t87, t83, 0, 0, 0, 0; t84, t80 * t81, 0, 0, 0, 0; 0, t79, 0, 0, 0, 0; t72, -t76 * t84, 0, t75, 0, 0; t74, -t76 * t87, 0, t73, 0, 0; 0, t86, 0, t79 * t77, 0, 0; -t73, -t77 * t84, 0, -t74, 0, 0; t75, -t77 * t87, 0, t72, 0, 0; 0, t85, 0, -t79 * t76, 0, 0;];
JR_rot  = t1;
