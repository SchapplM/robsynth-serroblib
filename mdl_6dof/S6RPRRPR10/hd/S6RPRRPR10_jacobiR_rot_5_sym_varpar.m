% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR10_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_jacobiR_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:04
% EndTime: 2019-02-26 21:06:04
% DurationCPUTime: 0.03s
% Computational Cost: add. (16->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
t77 = sin(qJ(3));
t78 = sin(qJ(1));
t88 = t78 * t77;
t79 = cos(qJ(4));
t87 = t78 * t79;
t76 = sin(qJ(4));
t80 = cos(qJ(3));
t86 = t80 * t76;
t85 = t80 * t79;
t81 = cos(qJ(1));
t84 = t81 * t77;
t83 = t81 * t79;
t82 = t81 * t80;
t75 = -t78 * t76 + t77 * t83;
t74 = t76 * t84 + t87;
t73 = t81 * t76 + t77 * t87;
t72 = t76 * t88 - t83;
t1 = [t75, 0, t78 * t85, -t72, 0, 0; t73, 0, -t79 * t82, t74, 0, 0; 0, 0, -t77 * t79, -t86, 0, 0; -t82, 0, t88, 0, 0, 0; -t78 * t80, 0, -t84, 0, 0, 0; 0, 0, t80, 0, 0, 0; t74, 0, t78 * t86, t73, 0, 0; t72, 0, -t76 * t82, -t75, 0, 0; 0, 0, -t77 * t76, t85, 0, 0;];
JR_rot  = t1;
