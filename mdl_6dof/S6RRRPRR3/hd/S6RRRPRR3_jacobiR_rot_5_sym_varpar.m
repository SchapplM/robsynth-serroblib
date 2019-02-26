% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:17
% EndTime: 2019-02-26 22:17:17
% DurationCPUTime: 0.04s
% Computational Cost: add. (68->14), mult. (76->8), div. (0->0), fcn. (122->6), ass. (0->14)
t60 = qJ(2) + qJ(3);
t58 = sin(t60);
t59 = cos(t60);
t61 = sin(qJ(5));
t63 = cos(qJ(5));
t68 = -t58 * t63 + t59 * t61;
t65 = t58 * t61 + t59 * t63;
t64 = cos(qJ(1));
t62 = sin(qJ(1));
t52 = t65 * t64;
t51 = t68 * t64;
t50 = t65 * t62;
t49 = t68 * t62;
t1 = [-t50, t51, t51, 0, -t51, 0; t52, t49, t49, 0, -t49, 0; 0, t65, t65, 0, -t65, 0; t49, t52, t52, 0, -t52, 0; -t51, t50, t50, 0, -t50, 0; 0, -t68, -t68, 0, t68, 0; -t64, 0, 0, 0, 0, 0; -t62, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
