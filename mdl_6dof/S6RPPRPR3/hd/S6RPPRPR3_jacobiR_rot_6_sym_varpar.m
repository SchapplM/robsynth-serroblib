% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRPR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:47
% EndTime: 2019-02-26 20:26:47
% DurationCPUTime: 0.03s
% Computational Cost: add. (61->16), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t80 = qJ(1) + pkin(9);
t76 = sin(t80);
t81 = sin(qJ(6));
t88 = t76 * t81;
t82 = cos(qJ(6));
t87 = t76 * t82;
t79 = qJ(4) + pkin(10);
t77 = cos(t79);
t86 = t77 * t81;
t85 = t77 * t82;
t78 = cos(t80);
t84 = t78 * t81;
t83 = t78 * t82;
t75 = sin(t79);
t74 = t75 * t83 - t88;
t73 = t75 * t84 + t87;
t72 = t75 * t87 + t84;
t71 = -t75 * t88 + t83;
t1 = [t74, 0, 0, t76 * t85, 0, t71; t72, 0, 0, -t77 * t83, 0, t73; 0, 0, 0, -t75 * t82, 0, -t86; -t73, 0, 0, -t76 * t86, 0, -t72; t71, 0, 0, t77 * t84, 0, t74; 0, 0, 0, t75 * t81, 0, -t85; -t78 * t77, 0, 0, t76 * t75, 0, 0; -t76 * t77, 0, 0, -t78 * t75, 0, 0; 0, 0, 0, t77, 0, 0;];
JR_rot  = t1;
