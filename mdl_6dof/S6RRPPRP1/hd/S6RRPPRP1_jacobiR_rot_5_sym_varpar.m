% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRP1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:45
% EndTime: 2019-02-26 21:24:45
% DurationCPUTime: 0.03s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t80 = qJ(2) + pkin(9);
t76 = sin(t80);
t81 = sin(qJ(1));
t88 = t81 * t76;
t79 = pkin(10) + qJ(5);
t77 = cos(t79);
t87 = t81 * t77;
t78 = cos(t80);
t86 = t81 * t78;
t82 = cos(qJ(1));
t85 = t82 * t76;
t84 = t82 * t77;
t83 = t82 * t78;
t75 = sin(t79);
t74 = t81 * t75 + t77 * t83;
t73 = -t75 * t83 + t87;
t72 = t82 * t75 - t77 * t86;
t71 = t75 * t86 + t84;
t1 = [t72, -t76 * t84, 0, 0, t73, 0; t74, -t76 * t87, 0, 0, -t71, 0; 0, t78 * t77, 0, 0, -t76 * t75, 0; t71, t75 * t85, 0, 0, -t74, 0; t73, t75 * t88, 0, 0, t72, 0; 0, -t78 * t75, 0, 0, -t76 * t77, 0; -t88, t83, 0, 0, 0, 0; t85, t86, 0, 0, 0, 0; 0, t76, 0, 0, 0, 0;];
JR_rot  = t1;
